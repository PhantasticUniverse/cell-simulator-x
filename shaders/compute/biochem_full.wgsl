// Full biochemistry on GPU — Phase 11.2.C.
//
// Extends the 18-species glycolysis+shunt kernel (11.2.A/B) to all 38 species
// of the `FullyIntegratedSolver` minus the post-RK4 hemoglobin and pH
// (those are 11.2.D). One workgroup-thread per cell; one dispatch advances
// every cell by `n_steps` RK4 ms-steps.
//
// State layout: cell-major, 38 f32 per cell, matching `FullyIntegratedIndices`.
//   0–16  Glycolysis (Glucose .. Pi)
//   17    2,3-BPG
//   18    6-Phosphogluconolactone (6-PGL)
//   19    6-Phosphogluconate (6-PG)
//   20    Ribulose-5-phosphate (Ru5P)
//   21    Ribose-5-phosphate (R5P)
//   22    Xylulose-5-phosphate (Xu5P)
//   23    Sedoheptulose-7-phosphate (S7P)
//   24    Erythrose-4-phosphate (E4P)
//   25    NADPH
//   26    NADP+
//   27    GSH
//   28    GSSG
//   29    H2O2
//   30    Cytosolic Ca²⁺ — STORED IN µM (not mM) for numerical headroom
//   31    Glutamate (GSH precursor)
//   32    Cysteine (GSH precursor)
//   33    Glycine (GSH precursor)
//   34    γ-Glutamylcysteine
//   35    Cytosolic Na⁺
//   36    Cytosolic K⁺
//   37    Cytosolic Cl⁻
//
// Mixed-precision: f32 throughout. Glycolytic / PPP / glutathione species span
// 0.001–5 mM; ions span 5–150 mM; Ca is in µM (0.1–~1.0). All within f32
// headroom. CPU-side keeps f64. Host wrapper does f64↔f32 at upload/download.
//
// Out of scope (still): post-RK4 hemoglobin Adair binding (11.2.D), pH buffer
// (11.2.D), and the three inline homeostasis hacks
// (basal NADPH consumption, basal GSH oxidation, ATP regen — 11.2.E).

const N_SPECIES: u32 = 38u;

// === Glycolysis indices (0–17) ===
const GLUCOSE: u32 = 0u;
const G6P:     u32 = 1u;
const F6P:     u32 = 2u;
const FBP:     u32 = 3u;
const DHAP:    u32 = 4u;
const GAP:     u32 = 5u;
const BPG13:   u32 = 6u;
const PG3:     u32 = 7u;
const PG2:     u32 = 8u;
const PEP:     u32 = 9u;
const PYR:     u32 = 10u;
const LAC:     u32 = 11u;
const ATP:     u32 = 12u;
const ADP:     u32 = 13u;
const NAD:     u32 = 14u;
const NADH:    u32 = 15u;
const PI:      u32 = 16u;
const BPG23:   u32 = 17u;

// === Redox indices (18–34) ===
const PGL6:    u32 = 18u;
const PG6:     u32 = 19u;
const RU5P:    u32 = 20u;
const R5P:     u32 = 21u;
const XU5P:    u32 = 22u;
const S7P:     u32 = 23u;
const E4P:     u32 = 24u;
const NADPH:   u32 = 25u;
const NADP:    u32 = 26u;
const GSH:     u32 = 27u;
const GSSG:    u32 = 28u;
const H2O2:    u32 = 29u;
const CA:      u32 = 30u;  // µM units!
const GLU:     u32 = 31u;
const CYS:     u32 = 32u;
const GLY:     u32 = 33u;
const GGC:     u32 = 34u;  // γ-Glu-Cys

// === Ion indices (35–37) ===
const NA:      u32 = 35u;
const K_ION:   u32 = 36u;
const CL:      u32 = 37u;

// === Uniforms — world-level constants and dispatch parameters ===
struct U {
    n_cells: u32,
    n_steps: u32,
    dt_sec: f32,
    external_glucose_mM: f32,
    atp_consumption_mM_per_sec: f32,
    enable_glucose_transport: u32,
    enable_lactate_export: u32,
    max_change_mM: f32,
    min_concentration_mM: f32,
    membrane_tension_pN_per_nm: f32,
    enable_piezo1: u32,
    enable_ion_homeostasis: u32,
    h2o2_production_rate_mM_per_sec: f32,
    na_external_mM: f32,
    k_external_mM: f32,
    g_na_per_sec: f32,
    g_k_per_sec: f32,
    _pad0: f32,
    _pad1: f32,
    _pad2: f32,
}

@group(0) @binding(0) var<storage, read_write> state: array<f32>;
@group(0) @binding(1) var<uniform>             u: U;

// === Kinetic helpers (mirror src/biochemistry/enzyme.rs) ============

fn michaelis_menten(vmax: f32, km: f32, s: f32) -> f32 {
    if (s <= 0.0) { return 0.0; }
    return vmax * s / (km + s);
}

fn hill_kinetics(vmax: f32, k_half: f32, s: f32, n: f32) -> f32 {
    if (s <= 0.0) { return 0.0; }
    let s_n = pow(s, n);
    let k_n = pow(k_half, n);
    return vmax * s_n / (k_n + s_n);
}

fn reversible_mm(vmax_f: f32, km_s: f32, s: f32, vmax_r: f32, km_p: f32, p: f32) -> f32 {
    let s_term = s / km_s;
    let p_term = p / km_p;
    let num = vmax_f * s_term - vmax_r * p_term;
    let den = 1.0 + s_term + p_term;
    if (den <= 0.0) { return 0.0; }
    return num / den;
}

fn ordered_bi_bi(vmax: f32, km_a: f32, km_b: f32, ki_a: f32, a: f32, b: f32) -> f32 {
    if (a <= 0.0 || b <= 0.0) { return 0.0; }
    let num = vmax * a * b;
    let den = ki_a * km_b + km_b * a + km_a * b + a * b;
    if (den <= 0.0) { return 0.0; }
    return num / den;
}

fn reversible_ordered_bi_bi(
    vmax_f: f32, vmax_r: f32,
    km_b: f32, km_p: f32, ki_a: f32, ki_q: f32,
    a: f32, b: f32, p: f32, q: f32,
) -> f32 {
    let num = vmax_f * a * b / (ki_a * km_b) - vmax_r * p * q / (ki_q * km_p);
    let den = 1.0 + a / ki_a + b / km_b + a * b / (ki_a * km_b)
            + p / km_p + q / ki_q + p * q / (ki_q * km_p);
    if (den <= 0.0) { return 0.0; }
    return num / den;
}

// === Glycolysis enzyme rates (Phase 11.2.A) =========================

fn rate_hexokinase(glucose: f32, atp: f32, g6p: f32) -> f32 {
    let base = ordered_bi_bi(0.035, 0.1, 1.0, 0.1, glucose, atp);
    return base / (1.0 + g6p / 0.5);
}

fn rate_gpi(g6p: f32, f6p: f32) -> f32 {
    return reversible_mm(1.0, 0.5, g6p, 2.5, 0.1, f6p);
}

fn rate_pfk(f6p: f32, atp: f32) -> f32 {
    let f6p_term = hill_kinetics(1.0, 0.04, f6p, 2.5);
    let atp_term = michaelis_menten(1.0, 0.1, atp);
    return 0.2 * f6p_term * atp_term;
}

fn rate_aldolase(fbp: f32, dhap: f32, gap: f32) -> f32 {
    let mass_action_ratio: f32 = select(1.0, (dhap * gap) / (fbp * 7.8e-2), fbp > 1e-9);
    let v_forward = 0.3 * fbp / (0.01 + fbp);
    let driving = clamp(1.0 - mass_action_ratio, -1.0, 1.0);
    return v_forward * driving;
}

fn rate_tpi(dhap: f32, gap: f32) -> f32 {
    return 10.0 * max(dhap, 0.0) - 222.222 * max(gap, 0.0);
}

fn rate_gapdh(gap: f32, nad: f32, pi: f32, bpg: f32, nadh: f32) -> f32 {
    let km_gap = 0.05;
    let km_nad = 0.05;
    let km_pi = 1.0;
    let km_bpg = 0.001;
    let km_nadh = 0.01;
    let forward = 1.5 * (gap / km_gap) * (nad / km_nad) * (pi / km_pi);
    let reverse = 0.5 * (bpg / km_bpg) * (nadh / km_nadh);
    let denom = 1.0 + gap / km_gap + nad / km_nad + pi / km_pi
              + bpg / km_bpg + nadh / km_nadh;
    return (forward - reverse) / denom;
}

fn rate_pgk(bpg: f32, adp: f32, pg3: f32, atp: f32) -> f32 {
    return reversible_ordered_bi_bi(
        5.0, 0.1,
        0.2, 0.5,
        0.002, 0.5,
        bpg, adp, pg3, atp
    );
}

fn rate_pgm(pg3: f32, pg2: f32) -> f32 {
    return reversible_mm(2.0, 0.2, pg3, 11.0, 0.03, pg2);
}

fn rate_enolase(pg2: f32, pep: f32) -> f32 {
    return reversible_mm(1.5, 0.05, pg2, 3.0, 0.1, pep);
}

fn rate_pk(pep: f32, adp: f32) -> f32 {
    let pep_term = hill_kinetics(1.0, 0.2, pep, 2.0);
    let adp_term = michaelis_menten(1.0, 0.3, adp);
    return 1.5 * pep_term * adp_term;
}

fn rate_ldh(pyr: f32, nadh: f32, lac: f32, nad: f32) -> f32 {
    return reversible_ordered_bi_bi(
        10.0, 0.01,
        0.01, 10.0,
        0.2, 0.1,
        pyr, nadh, lac, nad
    );
}

// 2,3-BPG shunt (Rapoport-Luebering) — Phase 11.2.B.

fn rate_bpgm(bpg13: f32, bpg23: f32) -> f32 {
    let base = michaelis_menten(0.045, 0.004, bpg13);
    return base / (1.0 + bpg23 / 2.0);
}

fn rate_bpgp(bpg23: f32, pi: f32, pg3: f32) -> f32 {
    let km_apparent = 2.0 * (1.0 + pi / 0.5);
    let base = michaelis_menten(0.0015, km_apparent, bpg23);
    return base / (1.0 + pg3 / 1.0);
}

// === PPP enzyme rates (Phase 11.2.C) ================================

// G6PDH: G6P + NADP+ → 6-PGL + NADPH (NADPH competitively inhibits NADP+
// binding; vmax tuned in Phase 10 from 0.06 to 0.004).
fn rate_g6pdh(g6p: f32, nadp: f32, nadph: f32) -> f32 {
    if (g6p <= 0.0 || nadp <= 0.0) { return 0.0; }
    let vmax: f32 = 0.004;
    let km_g6p: f32 = 0.039;
    let km_nadp: f32 = 0.005;
    let ki_nadph: f32 = 0.005;
    let km_nadp_app = km_nadp * (1.0 + nadph / ki_nadph);
    let g6p_term = g6p / (km_g6p + g6p);
    let nadp_term = nadp / (km_nadp_app + nadp);
    return vmax * g6p_term * nadp_term;
}

// 6-Phosphogluconolactonase: first-order k = 10/s.
fn rate_pgl_enzyme(pgl6: f32) -> f32 {
    return 10.0 * pgl6;
}

// 6PGDH: 6-PG + NADP+ → Ru5P + NADPH (CO2 not tracked).
fn rate_pgdh6(pg6: f32, nadp: f32) -> f32 {
    if (pg6 <= 0.0 || nadp <= 0.0) { return 0.0; }
    let vmax: f32 = 0.003;
    let km_pg6: f32 = 0.035;
    let km_nadp: f32 = 0.008;
    let pg6_term = pg6 / (km_pg6 + pg6);
    let nadp_term = nadp / (km_nadp + nadp);
    return vmax * pg6_term * nadp_term;
}

// RPE: Ru5P ⇌ Xu5P (Keq=1.4, vmax_r = vmax_f*km_xu5p/(Keq*km_ru5p) = 1.4286).
fn rate_rpe(ru5p: f32, xu5p: f32) -> f32 {
    return reversible_mm(2.0, 0.5, ru5p, 1.4285714, 0.5, xu5p);
}

// RPI: Ru5P ⇌ R5P (Keq=3.0, vmax_r = 2.0 * 1.5 / (3.0 * 0.8) = 1.25).
fn rate_rpi(ru5p: f32, r5p: f32) -> f32 {
    return reversible_mm(2.0, 0.8, ru5p, 1.25, 1.5, r5p);
}

// TK: Xu5P + R5P ⇌ G3P + S7P (bi-bi with mass-action equilibrium correction).
fn rate_tk(xu5p: f32, r5p: f32, g3p: f32, s7p: f32) -> f32 {
    let vmax: f32 = 1.0;
    let km_xu5p: f32 = 0.15;
    let km_r5p: f32 = 0.5;
    let keq: f32 = 1.2;
    let forward = vmax * xu5p * r5p / ((km_xu5p + xu5p) * (km_r5p + r5p));
    let substrate_product = xu5p * r5p;
    let mass_action = select(0.0, (g3p * s7p) / (substrate_product * keq), substrate_product > 1e-12);
    return forward * clamp(1.0 - mass_action, -1.0, 1.0);
}

// TA: S7P + G3P ⇌ E4P + F6P (same form as TK).
fn rate_ta(s7p: f32, g3p: f32, e4p: f32, f6p: f32) -> f32 {
    let vmax: f32 = 0.8;
    let km_s7p: f32 = 0.3;
    let km_g3p: f32 = 0.1;
    let keq: f32 = 1.05;
    let forward = vmax * s7p * g3p / ((km_s7p + s7p) * (km_g3p + g3p));
    let substrate_product = s7p * g3p;
    let mass_action = select(0.0, (e4p * f6p) / (substrate_product * keq), substrate_product > 1e-12);
    return forward * clamp(1.0 - mass_action, -1.0, 1.0);
}

// === Glutathione enzyme rates (Phase 11.2.C) ========================

// GPx: 2 GSH + H2O2 → GSSG + 2 H2O (ping-pong bi-bi).
fn rate_gpx(gsh: f32, h2o2: f32) -> f32 {
    if (gsh <= 0.0 || h2o2 <= 0.0) { return 0.0; }
    let vmax: f32 = 0.02;
    let km_gsh: f32 = 1.0;
    let km_h2o2: f32 = 0.002;
    let num = vmax * gsh * h2o2;
    let den = km_gsh * h2o2 + km_h2o2 * gsh + gsh * h2o2;
    if (den <= 0.0) { return 0.0; }
    return num / den;
}

// GR: GSSG + NADPH → 2 GSH + NADP+ (vmax tuned in Phase 10 from 0.15 to 0.025).
fn rate_gr(gssg: f32, nadph: f32) -> f32 {
    if (gssg <= 0.0 || nadph <= 0.0) { return 0.0; }
    let vmax: f32 = 0.025;
    let km_gssg: f32 = 0.015;
    let km_nadph: f32 = 0.015;
    let gssg_term = gssg / (km_gssg + gssg);
    let nadph_term = nadph / (km_nadph + nadph);
    return vmax * gssg_term * nadph_term;
}

// γ-GCS: Glu + Cys + ATP → γ-Glu-Cys + ADP (with non-competitive GSH feedback).
fn rate_gamma_gcs(glu: f32, cys: f32, atp: f32, gsh: f32) -> f32 {
    if (glu <= 0.0 || cys <= 0.0 || atp <= 0.0) { return 0.0; }
    let vmax: f32 = 0.01;
    let km_glu: f32 = 1.8;
    let km_cys: f32 = 0.3;
    let km_atp: f32 = 0.4;
    let ki_gsh: f32 = 2.3;
    let glu_term = glu / (km_glu + glu);
    let cys_term = cys / (km_cys + cys);
    let atp_term = atp / (km_atp + atp);
    let inhib = 1.0 / (1.0 + gsh / ki_gsh);
    return vmax * glu_term * cys_term * atp_term * inhib;
}

// GS: γ-Glu-Cys + Gly + ATP → GSH + ADP.
fn rate_gs(ggc: f32, gly: f32, atp: f32) -> f32 {
    if (ggc <= 0.0 || gly <= 0.0 || atp <= 0.0) { return 0.0; }
    let vmax: f32 = 0.02;
    let km_ggc: f32 = 0.15;
    let km_gly: f32 = 0.4;
    let km_atp: f32 = 0.3;
    let ggc_term = ggc / (km_ggc + ggc);
    let gly_term = gly / (km_gly + gly);
    let atp_term = atp / (km_atp + atp);
    return vmax * ggc_term * gly_term * atp_term;
}

// === Piezo1 mechanotransduction + PMCA (Phase 11.2.C) ===============

// Open probability via Hill kinetics on tension.
fn piezo1_open_prob(tension: f32) -> f32 {
    if (tension <= 0.0) { return 0.0; }
    let half_act: f32 = 1.5;  // pN/nm
    let hill_n: f32 = 2.0;
    let x = tension / half_act;
    let x_n = pow(x, hill_n);
    return x_n / (1.0 + x_n);
}

// Ca²⁺ influx in µM/s (driving-force × open-probability).
fn piezo1_ca_influx_uM_per_sec(tension: f32, ca_uM: f32) -> f32 {
    let p_open = piezo1_open_prob(tension);
    let max_perm: f32 = 0.5;
    let ca_ext_uM: f32 = 1800.0;  // 1.8 mM in µM
    let driving = max(ca_ext_uM - ca_uM, 0.0);
    let driving_term = driving / (500.0 + driving);
    return p_open * max_perm * driving_term * 1000.0;
}

// PMCA Ca²⁺ extrusion in µM/s (ATP-dependent).
fn pmca_extrusion_uM_per_sec(ca_uM: f32, atp_mM: f32) -> f32 {
    let basal: f32 = 0.1;  // µM
    let extrusion_rate: f32 = 0.5;  // s⁻¹
    let km_atp: f32 = 0.1;  // mM
    let elevation = max(ca_uM - basal, 0.0);
    let atp_term = atp_mM / (km_atp + atp_mM);
    return extrusion_rate * elevation * atp_term;
}

// ATP release triggered by elevated cytosolic Ca²⁺ (Pannexin-1, Hill n=2).
fn atp_release_mM_per_sec(ca_uM: f32) -> f32 {
    let basal: f32 = 0.1;
    let half: f32 = 0.5;
    let max_release: f32 = 0.05;
    let elevation = max(ca_uM - basal, 0.0);
    let x = elevation / half;
    let x2 = x * x;
    return max_release * x2 / (1.0 + x2);
}

// === Na⁺/K⁺-ATPase (Phase 11.2.C) ===================================

// Pump rate (mM ATP / s); Hill on Na, Hill on extracellular K, MM on ATP.
fn nakatpase_rate(na_cyt: f32, atp_mM: f32) -> f32 {
    if (na_cyt <= 0.0 || u.k_external_mM <= 0.0 || atp_mM <= 0.0) { return 0.0; }
    let vmax: f32 = 0.055;
    let km_na: f32 = 15.0;
    let km_k: f32 = 1.5;
    let km_atp: f32 = 0.2;
    let hill_na: f32 = 3.0;
    let hill_k: f32 = 2.0;
    let na_n = pow(na_cyt, hill_na);
    let km_na_n = pow(km_na, hill_na);
    let na_term = na_n / (km_na_n + na_n);
    let k_n = pow(u.k_external_mM, hill_k);
    let km_k_n = pow(km_k, hill_k);
    let k_term = k_n / (km_k_n + k_n);
    let atp_term = atp_mM / (km_atp + atp_mM);
    return vmax * na_term * k_term * atp_term;
}

// === Derivatives ====================================================

fn derivatives(y: ptr<function, array<f32, 38>>, dydt: ptr<function, array<f32, 38>>) {
    for (var i = 0u; i < N_SPECIES; i = i + 1u) {
        (*dydt)[i] = 0.0;
    }

    // Read state into locals.
    let glucose = (*y)[GLUCOSE];
    let g6p     = (*y)[G6P];
    let f6p     = (*y)[F6P];
    let fbp     = (*y)[FBP];
    let dhap    = (*y)[DHAP];
    let gap     = (*y)[GAP];
    let bpg13   = (*y)[BPG13];
    let pg3     = (*y)[PG3];
    let pg2     = (*y)[PG2];
    let pep     = (*y)[PEP];
    let pyr     = (*y)[PYR];
    let lac     = (*y)[LAC];
    let atp     = (*y)[ATP];
    let adp     = (*y)[ADP];
    let nad     = (*y)[NAD];
    let nadh    = (*y)[NADH];
    let pi      = (*y)[PI];
    let bpg23   = (*y)[BPG23];
    let pgl6    = (*y)[PGL6];
    let pg6     = (*y)[PG6];
    let ru5p    = (*y)[RU5P];
    let r5p     = (*y)[R5P];
    let xu5p    = (*y)[XU5P];
    let s7p     = (*y)[S7P];
    let e4p     = (*y)[E4P];
    let nadph   = (*y)[NADPH];
    let nadp    = (*y)[NADP];
    let gsh     = (*y)[GSH];
    let gssg    = (*y)[GSSG];
    let h2o2    = (*y)[H2O2];
    let ca      = (*y)[CA];
    let glu     = (*y)[GLU];
    let cys     = (*y)[CYS];
    let gly     = (*y)[GLY];
    let ggc     = (*y)[GGC];
    let na      = (*y)[NA];
    let k_ion   = (*y)[K_ION];

    // === Glycolysis (11 enzymes) ===
    let v_hk = rate_hexokinase(glucose, atp, g6p);
    (*dydt)[GLUCOSE] = (*dydt)[GLUCOSE] - v_hk;
    (*dydt)[ATP]     = (*dydt)[ATP]     - v_hk;
    (*dydt)[G6P]     = (*dydt)[G6P]     + v_hk;
    (*dydt)[ADP]     = (*dydt)[ADP]     + v_hk;

    let v_gpi = rate_gpi(g6p, f6p);
    (*dydt)[G6P] = (*dydt)[G6P] - v_gpi;
    (*dydt)[F6P] = (*dydt)[F6P] + v_gpi;

    let v_pfk = rate_pfk(f6p, atp);
    (*dydt)[F6P] = (*dydt)[F6P] - v_pfk;
    (*dydt)[ATP] = (*dydt)[ATP] - v_pfk;
    (*dydt)[FBP] = (*dydt)[FBP] + v_pfk;
    (*dydt)[ADP] = (*dydt)[ADP] + v_pfk;

    let v_ald = rate_aldolase(fbp, dhap, gap);
    (*dydt)[FBP]  = (*dydt)[FBP]  - v_ald;
    (*dydt)[DHAP] = (*dydt)[DHAP] + v_ald;
    (*dydt)[GAP]  = (*dydt)[GAP]  + v_ald;

    let v_tpi = rate_tpi(dhap, gap);
    (*dydt)[DHAP] = (*dydt)[DHAP] - v_tpi;
    (*dydt)[GAP]  = (*dydt)[GAP]  + v_tpi;

    let v_gapdh = rate_gapdh(gap, nad, pi, bpg13, nadh);
    (*dydt)[GAP]   = (*dydt)[GAP]   - v_gapdh;
    (*dydt)[NAD]   = (*dydt)[NAD]   - v_gapdh;
    (*dydt)[PI]    = (*dydt)[PI]    - v_gapdh;
    (*dydt)[BPG13] = (*dydt)[BPG13] + v_gapdh;
    (*dydt)[NADH]  = (*dydt)[NADH]  + v_gapdh;

    let v_pgk = rate_pgk(bpg13, adp, pg3, atp);
    (*dydt)[BPG13] = (*dydt)[BPG13] - v_pgk;
    (*dydt)[ADP]   = (*dydt)[ADP]   - v_pgk;
    (*dydt)[PG3]   = (*dydt)[PG3]   + v_pgk;
    (*dydt)[ATP]   = (*dydt)[ATP]   + v_pgk;

    let v_pgm = rate_pgm(pg3, pg2);
    (*dydt)[PG3] = (*dydt)[PG3] - v_pgm;
    (*dydt)[PG2] = (*dydt)[PG2] + v_pgm;

    let v_eno = rate_enolase(pg2, pep);
    (*dydt)[PG2] = (*dydt)[PG2] - v_eno;
    (*dydt)[PEP] = (*dydt)[PEP] + v_eno;

    let v_pk = rate_pk(pep, adp);
    (*dydt)[PEP] = (*dydt)[PEP] - v_pk;
    (*dydt)[ADP] = (*dydt)[ADP] - v_pk;
    (*dydt)[PYR] = (*dydt)[PYR] + v_pk;
    (*dydt)[ATP] = (*dydt)[ATP] + v_pk;

    let v_ldh = rate_ldh(pyr, nadh, lac, nad);
    (*dydt)[PYR]  = (*dydt)[PYR]  - v_ldh;
    (*dydt)[NADH] = (*dydt)[NADH] - v_ldh;
    (*dydt)[LAC]  = (*dydt)[LAC]  + v_ldh;
    (*dydt)[NAD]  = (*dydt)[NAD]  + v_ldh;

    // === 2,3-BPG shunt ===
    let v_bpgm = rate_bpgm(bpg13, bpg23);
    (*dydt)[BPG13] = (*dydt)[BPG13] - v_bpgm;
    (*dydt)[BPG23] = (*dydt)[BPG23] + v_bpgm;

    let v_bpgp = rate_bpgp(bpg23, pi, pg3);
    (*dydt)[BPG23] = (*dydt)[BPG23] - v_bpgp;
    (*dydt)[PG3]   = (*dydt)[PG3]   + v_bpgp;
    (*dydt)[PI]    = (*dydt)[PI]    + v_bpgp;

    // === PPP (7 enzymes) ===
    let v_g6pdh = rate_g6pdh(g6p, nadp, nadph);
    (*dydt)[G6P]   = (*dydt)[G6P]   - v_g6pdh;
    (*dydt)[NADP]  = (*dydt)[NADP]  - v_g6pdh;
    (*dydt)[PGL6]  = (*dydt)[PGL6]  + v_g6pdh;
    (*dydt)[NADPH] = (*dydt)[NADPH] + v_g6pdh;

    let v_pgl = rate_pgl_enzyme(pgl6);
    (*dydt)[PGL6] = (*dydt)[PGL6] - v_pgl;
    (*dydt)[PG6]  = (*dydt)[PG6]  + v_pgl;

    let v_pgdh6 = rate_pgdh6(pg6, nadp);
    (*dydt)[PG6]   = (*dydt)[PG6]   - v_pgdh6;
    (*dydt)[NADP]  = (*dydt)[NADP]  - v_pgdh6;
    (*dydt)[RU5P]  = (*dydt)[RU5P]  + v_pgdh6;
    (*dydt)[NADPH] = (*dydt)[NADPH] + v_pgdh6;

    let v_rpe = rate_rpe(ru5p, xu5p);
    (*dydt)[RU5P] = (*dydt)[RU5P] - v_rpe;
    (*dydt)[XU5P] = (*dydt)[XU5P] + v_rpe;

    let v_rpi = rate_rpi(ru5p, r5p);
    (*dydt)[RU5P] = (*dydt)[RU5P] - v_rpi;
    (*dydt)[R5P]  = (*dydt)[R5P]  + v_rpi;

    let v_tk = rate_tk(xu5p, r5p, gap, s7p);
    (*dydt)[XU5P] = (*dydt)[XU5P] - v_tk;
    (*dydt)[R5P]  = (*dydt)[R5P]  - v_tk;
    (*dydt)[GAP]  = (*dydt)[GAP]  + v_tk;
    (*dydt)[S7P]  = (*dydt)[S7P]  + v_tk;

    let v_ta = rate_ta(s7p, gap, e4p, f6p);
    (*dydt)[S7P]  = (*dydt)[S7P]  - v_ta;
    (*dydt)[GAP]  = (*dydt)[GAP]  - v_ta;
    (*dydt)[E4P]  = (*dydt)[E4P]  + v_ta;
    (*dydt)[F6P]  = (*dydt)[F6P]  + v_ta;

    // === Glutathione cycle (4 enzymes) ===
    let v_gpx = rate_gpx(gsh, h2o2);
    (*dydt)[GSH]  = (*dydt)[GSH]  - 2.0 * v_gpx;
    (*dydt)[H2O2] = (*dydt)[H2O2] - v_gpx;
    (*dydt)[GSSG] = (*dydt)[GSSG] + v_gpx;

    let v_gr = rate_gr(gssg, nadph);
    (*dydt)[GSSG]  = (*dydt)[GSSG]  - v_gr;
    (*dydt)[NADPH] = (*dydt)[NADPH] - v_gr;
    (*dydt)[GSH]   = (*dydt)[GSH]   + 2.0 * v_gr;
    (*dydt)[NADP]  = (*dydt)[NADP]  + v_gr;

    let v_gamma_gcs = rate_gamma_gcs(glu, cys, atp, gsh);
    (*dydt)[GLU] = (*dydt)[GLU] - v_gamma_gcs;
    (*dydt)[CYS] = (*dydt)[CYS] - v_gamma_gcs;
    (*dydt)[ATP] = (*dydt)[ATP] - v_gamma_gcs;
    (*dydt)[GGC] = (*dydt)[GGC] + v_gamma_gcs;
    (*dydt)[ADP] = (*dydt)[ADP] + v_gamma_gcs;

    let v_gs = rate_gs(ggc, gly, atp);
    (*dydt)[GGC] = (*dydt)[GGC] - v_gs;
    (*dydt)[GLY] = (*dydt)[GLY] - v_gs;
    (*dydt)[ATP] = (*dydt)[ATP] - v_gs;
    (*dydt)[GSH] = (*dydt)[GSH] + v_gs;
    (*dydt)[ADP] = (*dydt)[ADP] + v_gs;

    // Basal H2O2 production from Hb autoxidation.
    (*dydt)[H2O2] = (*dydt)[H2O2] + u.h2o2_production_rate_mM_per_sec;

    // === Piezo1 + PMCA + Ca²⁺-triggered ATP release ===
    if (u.enable_piezo1 != 0u) {
        let influx = piezo1_ca_influx_uM_per_sec(u.membrane_tension_pN_per_nm, ca);
        let extrusion = pmca_extrusion_uM_per_sec(ca, atp);
        (*dydt)[CA] = (*dydt)[CA] + influx - extrusion;
        // PMCA: 1 Ca²⁺ extruded = 1 ATP. Convert µM/s → mM/s.
        let pmca_atp = extrusion / 1000.0;
        (*dydt)[ATP] = (*dydt)[ATP] - pmca_atp;
        (*dydt)[ADP] = (*dydt)[ADP] + pmca_atp;
        // Pannexin-1 ATP release.
        let release = atp_release_mM_per_sec(ca);
        (*dydt)[ATP] = (*dydt)[ATP] - release;
    }

    // === Ion homeostasis (Na⁺/K⁺-ATPase + leaks) ===
    if (u.enable_ion_homeostasis != 0u) {
        let pump = nakatpase_rate(na, atp);
        (*dydt)[NA]    = (*dydt)[NA]    - 3.0 * pump;
        (*dydt)[K_ION] = (*dydt)[K_ION] + 2.0 * pump;
        (*dydt)[ATP]   = (*dydt)[ATP]   - pump;
        (*dydt)[ADP]   = (*dydt)[ADP]   + pump;
        // Passive leaks (Na⁺ inward, K⁺ outward).
        let na_leak = u.g_na_per_sec * (u.na_external_mM - na);
        (*dydt)[NA] = (*dydt)[NA] + na_leak;
        let k_leak = u.g_k_per_sec * (k_ion - u.k_external_mM);
        (*dydt)[K_ION] = (*dydt)[K_ION] - k_leak;
    }

    // === World-level dynamics (external ATP demand, GLUT1, MCT1) ===
    (*dydt)[ATP] = (*dydt)[ATP] - u.atp_consumption_mM_per_sec;
    (*dydt)[ADP] = (*dydt)[ADP] + u.atp_consumption_mM_per_sec;

    if (u.enable_glucose_transport != 0u) {
        let transport = 0.5 * (u.external_glucose_mM - glucose);
        (*dydt)[GLUCOSE] = (*dydt)[GLUCOSE] + transport;
    }

    if (u.enable_lactate_export != 0u && lac > 1.0) {
        let export_rate = 0.1 * (lac - 1.0);
        (*dydt)[LAC] = (*dydt)[LAC] - export_rate;
    }
}

// === RK4 step =======================================================

fn apply_clamp_and_floor(y_old: f32, dy: f32) -> f32 {
    let dy_clamped = clamp(dy, -u.max_change_mM, u.max_change_mM);
    return max(y_old + dy_clamped, u.min_concentration_mM);
}

@compute @workgroup_size(64)
fn rk4_step(@builtin(global_invocation_id) gid: vec3<u32>) {
    let cell = gid.x;
    if (cell >= u.n_cells) { return; }

    let base = cell * N_SPECIES;

    var y: array<f32, 38>;
    for (var i = 0u; i < N_SPECIES; i = i + 1u) {
        y[i] = state[base + i];
    }

    for (var step = 0u; step < u.n_steps; step = step + 1u) {
        var k1: array<f32, 38>;
        var k2: array<f32, 38>;
        var k3: array<f32, 38>;
        var k4: array<f32, 38>;
        var y_temp: array<f32, 38>;

        derivatives(&y, &k1);

        for (var i = 0u; i < N_SPECIES; i = i + 1u) {
            y_temp[i] = y[i] + 0.5 * u.dt_sec * k1[i];
        }
        derivatives(&y_temp, &k2);

        for (var i = 0u; i < N_SPECIES; i = i + 1u) {
            y_temp[i] = y[i] + 0.5 * u.dt_sec * k2[i];
        }
        derivatives(&y_temp, &k3);

        for (var i = 0u; i < N_SPECIES; i = i + 1u) {
            y_temp[i] = y[i] + u.dt_sec * k3[i];
        }
        derivatives(&y_temp, &k4);

        let dt6 = u.dt_sec / 6.0;
        for (var i = 0u; i < N_SPECIES; i = i + 1u) {
            let dy = dt6 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
            y[i] = apply_clamp_and_floor(y[i], dy);
        }
    }

    for (var i = 0u; i < N_SPECIES; i = i + 1u) {
        state[base + i] = y[i];
    }
}
