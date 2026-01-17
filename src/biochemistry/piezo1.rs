//! Piezo1 mechanosensitive channel for RBC mechanotransduction.
//!
//! Piezo1 is a mechanosensitive cation channel expressed in RBCs that opens
//! in response to membrane tension, allowing Ca2+ influx. This Ca2+ elevation
//! triggers downstream signaling including ATP release via Pannexin-1 channels.
//!
//! In the context of RBCs:
//! - Piezo1 mutations cause hereditary xerocytosis (dehydrated stomatocytosis)
//! - Channel opening under shear stress may regulate RBC volume and deformability
//! - Ca2+ influx triggers ATP release which has autocrine/paracrine signaling roles
//!
//! References:
//! - Cahalan SM et al. Science. 2015;348:1261-1264 (RBC Piezo1)
//! - Dyrda A et al. J Biol Chem. 2010;285:28210-28219 (Ca2+ and ATP release)
//! - Coste B et al. Science. 2010;330:55-60 (Piezo1 discovery)
//! - Engelmann B et al. Mol Membr Biol. 2008;25:523-531 (RBC Ca2+ signaling)

use super::MetabolitePool;
use super::pentose_phosphate::RedoxIndices;

/// Piezo1 mechanosensitive channel
///
/// Opens in response to membrane tension, allowing Ca2+ influx.
/// The channel has cooperative activation with Hill-like kinetics.
#[derive(Debug, Clone)]
pub struct Piezo1Channel {
    /// Half-activation tension (pN/nm)
    /// Reference: Cahalan 2015, ~1-2 pN/nm for RBC Piezo1
    pub half_activation_tension_pN_per_nm: f64,

    /// Hill coefficient for cooperative activation
    /// Reference: Lewis 2017, n ~ 2
    pub hill_coefficient: f64,

    /// Maximum Ca2+ permeability (uM/s per uM driving force)
    /// This represents the channel's permeability * open probability coefficient
    pub max_permeability_per_sec: f64,

    /// Extracellular Ca2+ concentration (mM)
    /// Reference: ~1.8 mM in plasma
    pub ca_external_mM: f64,

    /// Basal cytosolic Ca2+ (uM) - resting level
    /// Reference: Engelmann 2008, ~100 nM = 0.1 uM
    pub ca_basal_cytosolic_uM: f64,

    /// Ca2+ extrusion rate constant (1/s)
    /// Represents PMCA pump activity returning Ca2+ to baseline
    pub ca_extrusion_rate_per_sec: f64,

    /// Km for ATP in PMCA (mM)
    /// Reference: ~0.1-0.2 mM for PMCA
    pub pmca_km_atp_mM: f64,

    /// Index for cytosolic Ca2+ in metabolite pool
    idx_ca: usize,

    /// Index for ATP in metabolite pool
    idx_atp: usize,

    /// Index for ADP in metabolite pool
    idx_adp: usize,
}

impl Piezo1Channel {
    /// Create a new Piezo1 channel with default parameters
    pub fn new(indices: &RedoxIndices) -> Self {
        Self {
            // Reference: Cahalan 2015, Lewis 2017
            half_activation_tension_pN_per_nm: 1.5,
            hill_coefficient: 2.0,
            max_permeability_per_sec: 0.5,  // uM/s at saturating open probability
            ca_external_mM: 1.8,
            ca_basal_cytosolic_uM: 0.1,  // 100 nM
            ca_extrusion_rate_per_sec: 0.5,  // Tau ~ 2s for return to baseline
            pmca_km_atp_mM: 0.1,  // Km for ATP in PMCA
            idx_ca: indices.ca2_plus_cytosolic,
            idx_atp: indices.atp,
            idx_adp: indices.adp,
        }
    }

    /// Calculate Piezo1 open probability from membrane tension
    ///
    /// Uses Hill equation:
    /// P_open = (T/T_half)^n / (1 + (T/T_half)^n)
    ///
    /// # Arguments
    /// * `tension_pN_per_nm` - Membrane tension in pN/nm
    ///
    /// # Returns
    /// Open probability (0 to 1)
    pub fn open_probability(&self, tension_pN_per_nm: f64) -> f64 {
        if tension_pN_per_nm <= 0.0 {
            return 0.0;
        }

        let x = tension_pN_per_nm / self.half_activation_tension_pN_per_nm;
        let x_n = x.powf(self.hill_coefficient);

        x_n / (1.0 + x_n)
    }

    /// Calculate Ca2+ influx rate (uM/s)
    ///
    /// Influx depends on:
    /// - Channel open probability (tension-dependent)
    /// - Driving force (extracellular - cytosolic Ca2+)
    ///
    /// # Arguments
    /// * `tension_pN_per_nm` - Membrane tension
    /// * `ca_cytosolic_uM` - Current cytosolic Ca2+ concentration in uM
    ///
    /// # Returns
    /// Ca2+ influx rate in uM/s
    pub fn ca_influx_rate_uM_per_sec(&self, tension_pN_per_nm: f64, ca_cytosolic_uM: f64) -> f64 {
        let p_open = self.open_probability(tension_pN_per_nm);

        // Driving force: [Ca2+]ext (in uM) - [Ca2+]cyt
        let ca_external_uM = self.ca_external_mM * 1000.0;  // Convert mM to uM
        let driving_force = (ca_external_uM - ca_cytosolic_uM).max(0.0);

        // Normalize driving force to avoid excessive rates
        // Using Michaelis-like saturation with Km at ~500 uM
        let driving_term = driving_force / (500.0 + driving_force);

        p_open * self.max_permeability_per_sec * driving_term * 1000.0  // Scale factor
    }

    /// Calculate Ca2+ extrusion rate (uM/s) without ATP dependency
    ///
    /// PMCA pumps return Ca2+ to baseline level with first-order kinetics
    /// above the resting concentration.
    ///
    /// # Arguments
    /// * `ca_cytosolic_uM` - Current cytosolic Ca2+ concentration in uM
    ///
    /// # Returns
    /// Ca2+ extrusion rate in uM/s (positive = Ca2+ leaving cytosol)
    pub fn ca_extrusion_rate_uM_per_sec(&self, ca_cytosolic_uM: f64) -> f64 {
        let ca_elevation = (ca_cytosolic_uM - self.ca_basal_cytosolic_uM).max(0.0);
        self.ca_extrusion_rate_per_sec * ca_elevation
    }

    /// Calculate Ca2+ extrusion rate with ATP dependency (uM/s)
    ///
    /// PMCA is an ATP-dependent pump. This method calculates the extrusion
    /// rate considering ATP availability using Michaelis-Menten kinetics.
    ///
    /// Stoichiometry: 1 Ca2+ extruded = 1 ATP consumed
    ///
    /// # Arguments
    /// * `ca_cytosolic_uM` - Current cytosolic Ca2+ concentration in uM
    /// * `atp_mM` - ATP concentration in mM
    ///
    /// # Returns
    /// Ca2+ extrusion rate in uM/s (positive = Ca2+ leaving cytosol)
    ///
    /// Reference: Carafoli E. Physiol Rev. 1991;71:129-153
    pub fn ca_extrusion_rate_with_atp_uM_per_sec(&self, ca_cytosolic_uM: f64, atp_mM: f64) -> f64 {
        let ca_elevation = (ca_cytosolic_uM - self.ca_basal_cytosolic_uM).max(0.0);
        let atp_term = atp_mM / (self.pmca_km_atp_mM + atp_mM);
        self.ca_extrusion_rate_per_sec * ca_elevation * atp_term
    }

    /// Update cytosolic Ca2+ based on Piezo1 and pumps
    ///
    /// # Arguments
    /// * `metabolites` - Metabolite pool to update
    /// * `tension_pN_per_nm` - Current membrane tension
    /// * `dydt` - Derivatives array to update
    pub fn compute_ca_derivative(
        &self,
        metabolites: &MetabolitePool,
        tension_pN_per_nm: f64,
        dydt: &mut [f64],
    ) {
        let ca_uM = metabolites.get(self.idx_ca);

        // Net Ca2+ flux = influx - extrusion
        let influx = self.ca_influx_rate_uM_per_sec(tension_pN_per_nm, ca_uM);
        let extrusion = self.ca_extrusion_rate_uM_per_sec(ca_uM);

        dydt[self.idx_ca] += influx - extrusion;
    }

    /// Update cytosolic Ca2+ with ATP-dependent PMCA
    ///
    /// This method includes ATP coupling for PMCA:
    /// - Ca2+ extrusion depends on ATP availability
    /// - Each Ca2+ extruded consumes 1 ATP
    ///
    /// # Arguments
    /// * `metabolites` - Metabolite pool to update
    /// * `tension_pN_per_nm` - Current membrane tension
    /// * `dydt` - Derivatives array to update
    pub fn compute_ca_derivative_with_atp(
        &self,
        metabolites: &MetabolitePool,
        tension_pN_per_nm: f64,
        dydt: &mut [f64],
    ) {
        let ca_uM = metabolites.get(self.idx_ca);
        let atp_mM = metabolites.get(self.idx_atp);

        // Ca2+ influx (passive, via Piezo1)
        let influx = self.ca_influx_rate_uM_per_sec(tension_pN_per_nm, ca_uM);

        // Ca2+ extrusion (ATP-dependent PMCA)
        let extrusion = self.ca_extrusion_rate_with_atp_uM_per_sec(ca_uM, atp_mM);

        // Update Ca2+
        dydt[self.idx_ca] += influx - extrusion;

        // PMCA ATP consumption: 1 Ca2+ extruded = 1 ATP consumed
        // Convert uM/s to mM/s (divide by 1000)
        let atp_consumption_mM_per_sec = extrusion / 1000.0;
        dydt[self.idx_atp] -= atp_consumption_mM_per_sec;
        dydt[self.idx_adp] += atp_consumption_mM_per_sec;
    }

    /// Get PMCA ATP consumption rate (mM/s)
    pub fn pmca_atp_consumption_mM_per_sec(&self, ca_cytosolic_uM: f64, atp_mM: f64) -> f64 {
        let extrusion = self.ca_extrusion_rate_with_atp_uM_per_sec(ca_cytosolic_uM, atp_mM);
        extrusion / 1000.0  // Convert uM/s to mM/s
    }
}

/// Calculate ATP release rate triggered by Ca2+ elevation
///
/// In RBCs, elevated cytosolic Ca2+ triggers ATP release through Pannexin-1
/// channels. This has autocrine/paracrine signaling effects.
///
/// Uses Hill-like kinetics with Ca2+ elevation as the signal.
///
/// # Arguments
/// * `ca_cytosolic_uM` - Current cytosolic Ca2+ in uM
/// * `baseline_uM` - Resting Ca2+ level in uM
///
/// # Returns
/// ATP release rate in mM/s
///
/// Reference: Dyrda 2010, Locovei 2006
pub fn atp_release_rate_mM_per_sec(ca_cytosolic_uM: f64, baseline_uM: f64) -> f64 {
    let ca_elevation = (ca_cytosolic_uM - baseline_uM).max(0.0);

    // Half-max activation at ~500 nM elevation (0.5 uM)
    let half_activation_uM = 0.5;

    // Maximum ATP release rate
    // Reference: ~0.05-0.1 mM/s under strong stimulation
    let max_release_mM_per_sec = 0.05;

    // Hill kinetics (n=2 for cooperative activation)
    let x = ca_elevation / half_activation_uM;
    let x_2 = x * x;

    max_release_mM_per_sec * x_2 / (1.0 + x_2)
}

/// Piezo1 mechanotransduction system
///
/// Combines Piezo1 channel with downstream ATP release signaling.
pub struct Piezo1System {
    pub channel: Piezo1Channel,
    /// Enable ATP release triggered by Ca2+
    pub enable_atp_release: bool,
    /// Enable ATP-dependent PMCA (Ca2+ extrusion consumes ATP)
    pub enable_pmca_atp_coupling: bool,
}

impl Piezo1System {
    pub fn new(indices: &RedoxIndices) -> Self {
        Self {
            channel: Piezo1Channel::new(indices),
            enable_atp_release: true,
            enable_pmca_atp_coupling: true,  // Default: PMCA consumes ATP
        }
    }

    /// Compute all derivatives for Piezo1 system
    ///
    /// # Arguments
    /// * `metabolites` - Metabolite pool
    /// * `tension_pN_per_nm` - Current membrane tension
    /// * `dydt` - Derivatives array to update
    pub fn compute_derivatives(
        &self,
        metabolites: &MetabolitePool,
        tension_pN_per_nm: f64,
        dydt: &mut [f64],
    ) {
        // Ca2+ dynamics (with or without ATP coupling)
        if self.enable_pmca_atp_coupling {
            self.channel.compute_ca_derivative_with_atp(metabolites, tension_pN_per_nm, dydt);
        } else {
            self.channel.compute_ca_derivative(metabolites, tension_pN_per_nm, dydt);
        }

        // ATP release if enabled
        if self.enable_atp_release {
            let ca_uM = metabolites.get(self.channel.idx_ca);
            let atp_release = atp_release_rate_mM_per_sec(
                ca_uM,
                self.channel.ca_basal_cytosolic_uM,
            );

            // ATP is released (consumed from cytosol)
            dydt[self.channel.idx_atp] -= atp_release;
        }
    }

    /// Get diagnostic information
    pub fn diagnostics(&self, metabolites: &MetabolitePool, tension_pN_per_nm: f64) -> Piezo1Diagnostics {
        let ca_uM = metabolites.get(self.channel.idx_ca);
        let atp_mM = metabolites.get(self.channel.idx_atp);
        let p_open = self.channel.open_probability(tension_pN_per_nm);
        let ca_influx = self.channel.ca_influx_rate_uM_per_sec(tension_pN_per_nm, ca_uM);

        // Use ATP-dependent extrusion if coupling is enabled
        let ca_extrusion = if self.enable_pmca_atp_coupling {
            self.channel.ca_extrusion_rate_with_atp_uM_per_sec(ca_uM, atp_mM)
        } else {
            self.channel.ca_extrusion_rate_uM_per_sec(ca_uM)
        };

        let atp_release = atp_release_rate_mM_per_sec(ca_uM, self.channel.ca_basal_cytosolic_uM);
        let pmca_atp = self.channel.pmca_atp_consumption_mM_per_sec(ca_uM, atp_mM);

        Piezo1Diagnostics {
            membrane_tension_pN_per_nm: tension_pN_per_nm,
            open_probability: p_open,
            ca_cytosolic_uM: ca_uM,
            ca_influx_uM_per_sec: ca_influx,
            ca_extrusion_uM_per_sec: ca_extrusion,
            ca_net_flux_uM_per_sec: ca_influx - ca_extrusion,
            atp_release_mM_per_sec: atp_release,
            pmca_atp_consumption_mM_per_sec: pmca_atp,
        }
    }
}

/// Diagnostic information for Piezo1 system
#[derive(Debug, Clone)]
pub struct Piezo1Diagnostics {
    pub membrane_tension_pN_per_nm: f64,
    pub open_probability: f64,
    pub ca_cytosolic_uM: f64,
    pub ca_influx_uM_per_sec: f64,
    pub ca_extrusion_uM_per_sec: f64,
    pub ca_net_flux_uM_per_sec: f64,
    pub atp_release_mM_per_sec: f64,
    /// PMCA ATP consumption rate (mM/s)
    pub pmca_atp_consumption_mM_per_sec: f64,
}

impl Piezo1Diagnostics {
    pub fn print_summary(&self) {
        println!("=== Piezo1 System ===");
        println!("Membrane tension:   {:.2} pN/nm", self.membrane_tension_pN_per_nm);
        println!("Open probability:   {:.1}%", self.open_probability * 100.0);
        println!("Cytosolic Ca2+:     {:.2} uM ({:.0} nM)",
            self.ca_cytosolic_uM, self.ca_cytosolic_uM * 1000.0);
        println!("Ca2+ influx:        {:.4} uM/s", self.ca_influx_uM_per_sec);
        println!("Ca2+ extrusion:     {:.4} uM/s", self.ca_extrusion_uM_per_sec);
        println!("Ca2+ net flux:      {:.4} uM/s", self.ca_net_flux_uM_per_sec);
        println!("ATP release:        {:.6} mM/s", self.atp_release_mM_per_sec);
        println!("PMCA ATP usage:     {:.6} mM/s", self.pmca_atp_consumption_mM_per_sec);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::glycolysis::MetaboliteIndices;

    fn create_test_indices() -> RedoxIndices {
        let glycolysis = MetaboliteIndices::default();
        RedoxIndices::new(&glycolysis, 18)
    }

    fn create_test_metabolites(indices: &RedoxIndices) -> MetabolitePool {
        let total = 18 + RedoxIndices::new_metabolite_count();
        let mut pool = MetabolitePool::new(total);

        // Set basal Ca2+ (in uM)
        pool.set(indices.ca2_plus_cytosolic, 0.1);  // 100 nM

        // Set ATP
        pool.set(indices.atp, 2.0);
        pool.set(indices.adp, 0.25);

        pool
    }

    #[test]
    fn test_piezo1_open_probability() {
        let indices = create_test_indices();
        let channel = Piezo1Channel::new(&indices);

        // Zero tension = closed
        assert_eq!(channel.open_probability(0.0), 0.0);

        // At half-activation tension, P_open should be 0.5
        let p_at_half = channel.open_probability(channel.half_activation_tension_pN_per_nm);
        assert!((p_at_half - 0.5).abs() < 0.01, "P_open at T_half should be ~0.5");

        // High tension should approach 1.0
        let p_high = channel.open_probability(10.0);
        assert!(p_high > 0.9, "High tension should give P_open > 0.9");
    }

    #[test]
    fn test_piezo1_ca_influx() {
        let indices = create_test_indices();
        let channel = Piezo1Channel::new(&indices);

        // No tension = no influx
        let influx_zero = channel.ca_influx_rate_uM_per_sec(0.0, 0.1);
        assert_eq!(influx_zero, 0.0);

        // Some tension should give positive influx
        let influx = channel.ca_influx_rate_uM_per_sec(2.0, 0.1);
        assert!(influx > 0.0, "Tension should cause Ca2+ influx");
    }

    #[test]
    fn test_piezo1_ca_extrusion() {
        let indices = create_test_indices();
        let channel = Piezo1Channel::new(&indices);

        // At baseline, no extrusion
        let extrusion_baseline = channel.ca_extrusion_rate_uM_per_sec(
            channel.ca_basal_cytosolic_uM,
        );
        assert!(extrusion_baseline.abs() < 1e-10);

        // Above baseline, extrusion should increase
        let extrusion_elevated = channel.ca_extrusion_rate_uM_per_sec(1.0);  // 1 uM
        assert!(extrusion_elevated > 0.0);
    }

    #[test]
    fn test_atp_release() {
        // At baseline Ca2+, no release
        let release_baseline = atp_release_rate_mM_per_sec(0.1, 0.1);
        assert!(release_baseline.abs() < 1e-10);

        // Elevated Ca2+ should trigger release
        let release_elevated = atp_release_rate_mM_per_sec(1.0, 0.1);  // 1 uM Ca2+
        assert!(release_elevated > 0.0);

        // Higher Ca2+ should give more release
        let release_high = atp_release_rate_mM_per_sec(2.0, 0.1);
        assert!(release_high > release_elevated);
    }

    #[test]
    fn test_piezo1_system_derivatives() {
        let indices = create_test_indices();
        let system = Piezo1System::new(&indices);
        let metabolites = create_test_metabolites(&indices);

        let mut dydt = vec![0.0; 18 + RedoxIndices::new_metabolite_count()];

        // At zero tension, should have minimal effect
        system.compute_derivatives(&metabolites, 0.0, &mut dydt);
        assert!(dydt[indices.ca2_plus_cytosolic].abs() < 1e-6);

        // At high tension, should see Ca2+ increase
        let mut dydt2 = vec![0.0; 18 + RedoxIndices::new_metabolite_count()];
        system.compute_derivatives(&metabolites, 3.0, &mut dydt2);
        assert!(dydt2[indices.ca2_plus_cytosolic] > 0.0, "High tension should increase Ca2+");
    }

    #[test]
    fn test_piezo1_diagnostics() {
        let indices = create_test_indices();
        let system = Piezo1System::new(&indices);
        let metabolites = create_test_metabolites(&indices);

        let diag = system.diagnostics(&metabolites, 2.0);

        assert!(diag.open_probability > 0.0);
        assert!(diag.open_probability < 1.0);
        assert!(diag.ca_cytosolic_uM > 0.0);
    }
}
