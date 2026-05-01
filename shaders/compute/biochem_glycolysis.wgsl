// Glycolysis on GPU — Phase 11.2.A.
//
// Ports the 11-enzyme glycolytic backbone (HK → LDH) plus glucose import
// and lactate export to a wgpu compute kernel. Each thread handles one
// cell; one dispatch advances every cell by `n_steps` RK4 ms-steps.
//
// State layout: cell-major, 18 f32 per cell (indices 0–17 of the
// `FullyIntegratedIndices` layout).
//   0  Glucose
//   1  Glucose-6-phosphate (G6P)
//   2  Fructose-6-phosphate (F6P)
//   3  Fructose-1,6-bisphosphate (F1,6BP)
//   4  Dihydroxyacetone phosphate (DHAP)
//   5  Glyceraldehyde-3-phosphate (GAP)
//   6  1,3-bisphosphoglycerate (1,3-BPG)
//   7  3-phosphoglycerate (3-PG)
//   8  2-phosphoglycerate (2-PG)
//   9  Phosphoenolpyruvate (PEP)
//   10 Pyruvate
//   11 Lactate
//   12 ATP
//   13 ADP
//   14 NAD+
//   15 NADH
//   16 Pi
//   17 2,3-bisphosphoglycerate (2,3-BPG)
//
// Indices match Rust `MetaboliteIndices` from src/biochemistry/glycolysis.rs.
//
// Mixed-precision strategy: f32 throughout. The species span 0.001–5 mM so
// f32 (24-bit mantissa, ~7 sig figs) is adequate. The CPU side keeps f64
// and the host wrapper does f64↔f32 conversion at upload/download.

const N_SPECIES: u32 = 18u;

// Species indices.
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

// Uniforms — world-level constants and dispatch parameters.
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

// === Enzyme rate functions ==========================================

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
    // k_forward = 10, k_reverse = 10/0.045 = 222.22
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
        5.0, 0.1,        // vmax_f, vmax_r (Keq >> 1, strongly forward)
        0.2, 0.5,        // km_adp (b), km_3pg (p)
        0.002, 0.5,      // ki_a (≈km_bpg), ki_q (≈km_atp)
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
        0.01, 10.0,    // km_nadh (b), km_lac (p)
        0.2, 0.1,      // ki_pyr (a), ki_nad (q)
        pyr, nadh, lac, nad
    );
}

// 2,3-BPG shunt (Rapoport-Luebering).
//
// BPGM: 1,3-BPG → 2,3-BPG (with product inhibition by 2,3-BPG).
fn rate_bpgm(bpg13: f32, bpg23: f32) -> f32 {
    let base = michaelis_menten(0.045, 0.004, bpg13);
    return base / (1.0 + bpg23 / 2.0);
}

// BPGP: 2,3-BPG → 3-PG + Pi (with competitive Pi inhibition + 3-PG
// product inhibition).
fn rate_bpgp(bpg23: f32, pi: f32, pg3: f32) -> f32 {
    let km_apparent = 2.0 * (1.0 + pi / 0.5);
    let base = michaelis_menten(0.0015, km_apparent, bpg23);
    return base / (1.0 + pg3 / 1.0);
}

// === Derivatives ====================================================

// Compute dy/dt for the 18-species glycolytic system (incl. 2,3-BPG shunt).
fn derivatives(y: ptr<function, array<f32, 18>>, dydt: ptr<function, array<f32, 18>>) {
    // Zero output.
    for (var i = 0u; i < N_SPECIES; i = i + 1u) {
        (*dydt)[i] = 0.0;
    }

    // Read state into locals (WGSL doesn't allow indexed reads inside a
    // ptr without explicit temp variables on some backends).
    let glucose = (*y)[GLUCOSE];
    let g6p     = (*y)[G6P];
    let f6p     = (*y)[F6P];
    let fbp     = (*y)[FBP];
    let dhap    = (*y)[DHAP];
    let gap     = (*y)[GAP];
    let bpg     = (*y)[BPG13];
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

    // === Enzymes ===================================================

    // HK: Glucose + ATP → G6P + ADP
    let v_hk = rate_hexokinase(glucose, atp, g6p);
    (*dydt)[GLUCOSE] = (*dydt)[GLUCOSE] - v_hk;
    (*dydt)[ATP]     = (*dydt)[ATP]     - v_hk;
    (*dydt)[G6P]     = (*dydt)[G6P]     + v_hk;
    (*dydt)[ADP]     = (*dydt)[ADP]     + v_hk;

    // GPI: G6P ⇌ F6P
    let v_gpi = rate_gpi(g6p, f6p);
    (*dydt)[G6P] = (*dydt)[G6P] - v_gpi;
    (*dydt)[F6P] = (*dydt)[F6P] + v_gpi;

    // PFK: F6P + ATP → FBP + ADP
    let v_pfk = rate_pfk(f6p, atp);
    (*dydt)[F6P] = (*dydt)[F6P] - v_pfk;
    (*dydt)[ATP] = (*dydt)[ATP] - v_pfk;
    (*dydt)[FBP] = (*dydt)[FBP] + v_pfk;
    (*dydt)[ADP] = (*dydt)[ADP] + v_pfk;

    // Aldolase: FBP ⇌ DHAP + GAP
    let v_ald = rate_aldolase(fbp, dhap, gap);
    (*dydt)[FBP]  = (*dydt)[FBP]  - v_ald;
    (*dydt)[DHAP] = (*dydt)[DHAP] + v_ald;
    (*dydt)[GAP]  = (*dydt)[GAP]  + v_ald;

    // TPI: DHAP ⇌ GAP
    let v_tpi = rate_tpi(dhap, gap);
    (*dydt)[DHAP] = (*dydt)[DHAP] - v_tpi;
    (*dydt)[GAP]  = (*dydt)[GAP]  + v_tpi;

    // GAPDH: GAP + NAD+ + Pi ⇌ 1,3-BPG + NADH
    let v_gapdh = rate_gapdh(gap, nad, pi, bpg, nadh);
    (*dydt)[GAP]   = (*dydt)[GAP]   - v_gapdh;
    (*dydt)[NAD]   = (*dydt)[NAD]   - v_gapdh;
    (*dydt)[PI]    = (*dydt)[PI]    - v_gapdh;
    (*dydt)[BPG13] = (*dydt)[BPG13] + v_gapdh;
    (*dydt)[NADH]  = (*dydt)[NADH]  + v_gapdh;

    // PGK: 1,3-BPG + ADP ⇌ 3-PG + ATP
    let v_pgk = rate_pgk(bpg, adp, pg3, atp);
    (*dydt)[BPG13] = (*dydt)[BPG13] - v_pgk;
    (*dydt)[ADP]   = (*dydt)[ADP]   - v_pgk;
    (*dydt)[PG3]   = (*dydt)[PG3]   + v_pgk;
    (*dydt)[ATP]   = (*dydt)[ATP]   + v_pgk;

    // PGM: 3-PG ⇌ 2-PG
    let v_pgm = rate_pgm(pg3, pg2);
    (*dydt)[PG3] = (*dydt)[PG3] - v_pgm;
    (*dydt)[PG2] = (*dydt)[PG2] + v_pgm;

    // Enolase: 2-PG ⇌ PEP
    let v_eno = rate_enolase(pg2, pep);
    (*dydt)[PG2] = (*dydt)[PG2] - v_eno;
    (*dydt)[PEP] = (*dydt)[PEP] + v_eno;

    // PK: PEP + ADP → Pyruvate + ATP
    let v_pk = rate_pk(pep, adp);
    (*dydt)[PEP] = (*dydt)[PEP] - v_pk;
    (*dydt)[ADP] = (*dydt)[ADP] - v_pk;
    (*dydt)[PYR] = (*dydt)[PYR] + v_pk;
    (*dydt)[ATP] = (*dydt)[ATP] + v_pk;

    // LDH: Pyruvate + NADH ⇌ Lactate + NAD+
    let v_ldh = rate_ldh(pyr, nadh, lac, nad);
    (*dydt)[PYR]  = (*dydt)[PYR]  - v_ldh;
    (*dydt)[NADH] = (*dydt)[NADH] - v_ldh;
    (*dydt)[LAC]  = (*dydt)[LAC]  + v_ldh;
    (*dydt)[NAD]  = (*dydt)[NAD]  + v_ldh;

    // BPGM: 1,3-BPG → 2,3-BPG (Rapoport-Luebering shunt entry).
    let v_bpgm = rate_bpgm(bpg, bpg23);
    (*dydt)[BPG13] = (*dydt)[BPG13] - v_bpgm;
    (*dydt)[BPG23] = (*dydt)[BPG23] + v_bpgm;

    // BPGP: 2,3-BPG → 3-PG + Pi (shunt exit, slow).
    let v_bpgp = rate_bpgp(bpg23, pi, pg3);
    (*dydt)[BPG23] = (*dydt)[BPG23] - v_bpgp;
    (*dydt)[PG3]   = (*dydt)[PG3]   + v_bpgp;
    (*dydt)[PI]    = (*dydt)[PI]    + v_bpgp;

    // === World-level dynamics ======================================

    // External ATP consumption (membrane pumps, etc.).
    (*dydt)[ATP] = (*dydt)[ATP] - u.atp_consumption_mM_per_sec;
    (*dydt)[ADP] = (*dydt)[ADP] + u.atp_consumption_mM_per_sec;

    // GLUT1 facilitated diffusion: equilibrate with external glucose.
    if (u.enable_glucose_transport != 0u) {
        let transport = 0.5 * (u.external_glucose_mM - glucose);
        (*dydt)[GLUCOSE] = (*dydt)[GLUCOSE] + transport;
    }

    // MCT1 lactate export above threshold.
    if (u.enable_lactate_export != 0u && lac > 1.0) {
        let export_rate = 0.1 * (lac - 1.0);
        (*dydt)[LAC] = (*dydt)[LAC] - export_rate;
    }
}

// === RK4 step =======================================================

// Apply max-change clamp + min-concentration floor (mirrors the CPU
// integrator's per-step stability hacks).
fn apply_clamp_and_floor(y_old: f32, dy: f32) -> f32 {
    let dy_clamped = clamp(dy, -u.max_change_mM, u.max_change_mM);
    return max(y_old + dy_clamped, u.min_concentration_mM);
}

@compute @workgroup_size(64)
fn rk4_step(@builtin(global_invocation_id) gid: vec3<u32>) {
    let cell = gid.x;
    if (cell >= u.n_cells) { return; }

    let base = cell * N_SPECIES;

    // Read state into local registers.
    var y: array<f32, 18>;
    for (var i = 0u; i < N_SPECIES; i = i + 1u) {
        y[i] = state[base + i];
    }

    // Inner ms-step loop. Each iteration advances by u.dt_sec.
    for (var step = 0u; step < u.n_steps; step = step + 1u) {
        var k1: array<f32, 18>;
        var k2: array<f32, 18>;
        var k3: array<f32, 18>;
        var k4: array<f32, 18>;
        var y_temp: array<f32, 18>;

        // k1 = f(y)
        derivatives(&y, &k1);

        // k2 = f(y + 0.5*dt*k1)
        for (var i = 0u; i < N_SPECIES; i = i + 1u) {
            y_temp[i] = y[i] + 0.5 * u.dt_sec * k1[i];
        }
        derivatives(&y_temp, &k2);

        // k3 = f(y + 0.5*dt*k2)
        for (var i = 0u; i < N_SPECIES; i = i + 1u) {
            y_temp[i] = y[i] + 0.5 * u.dt_sec * k2[i];
        }
        derivatives(&y_temp, &k3);

        // k4 = f(y + dt*k3)
        for (var i = 0u; i < N_SPECIES; i = i + 1u) {
            y_temp[i] = y[i] + u.dt_sec * k3[i];
        }
        derivatives(&y_temp, &k4);

        // y_new = y + dt/6 * (k1 + 2k2 + 2k3 + k4), clamped + floored.
        let dt6 = u.dt_sec / 6.0;
        for (var i = 0u; i < N_SPECIES; i = i + 1u) {
            let dy = dt6 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
            y[i] = apply_clamp_and_floor(y[i], dy);
        }
    }

    // Write back.
    for (var i = 0u; i < N_SPECIES; i = i + 1u) {
        state[base + i] = y[i];
    }
}
