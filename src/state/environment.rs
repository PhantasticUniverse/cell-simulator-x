//! Environment state data structures.
//!
//! Represents the external environment surrounding the red blood cell,
//! including plasma conditions and flow properties.

/// Environmental state surrounding the cell
#[derive(Debug, Clone)]
pub struct EnvironmentState {
    /// Plasma conditions
    pub plasma: PlasmaState,
    /// Flow conditions
    pub flow: FlowState,
}

impl Default for EnvironmentState {
    fn default() -> Self {
        Self {
            plasma: PlasmaState::default(),
            flow: FlowState::default(),
        }
    }
}

/// Plasma composition and properties
#[derive(Debug, Clone)]
pub struct PlasmaState {
    /// Oxygen partial pressure (mmHg)
    /// Reference: 95-100 mmHg arterial, 40 mmHg venous
    /// Source: West, Respiratory Physiology, 2012
    pub po2_mmHg: f32,
    /// Carbon dioxide partial pressure (mmHg)
    /// Reference: 40 mmHg arterial, 46 mmHg venous
    /// Source: West, Respiratory Physiology, 2012
    pub pco2_mmHg: f32,
    /// Plasma pH
    /// Reference: 7.35-7.45
    /// Source: Davenport, The ABC of Acid-Base Chemistry, 1974
    pub ph: f32,
    /// Temperature (°C)
    pub temperature_celsius: f32,
    /// Plasma glucose concentration (mM)
    /// Reference: 4.0-6.0 mM fasting
    /// Source: American Diabetes Association standards
    pub glucose_mM: f32,
    /// Plasma osmolarity (mOsm/L)
    /// Reference: 280-295 mOsm/L
    /// Source: Guyton & Hall, Textbook of Medical Physiology
    pub osmolarity_mOsm_per_L: f32,
    /// Plasma viscosity (mPa·s)
    /// Reference: ~1.2 mPa·s at 37°C
    /// Source: Késmárky et al., Clin Hemorheol Microcirc 2008
    pub viscosity_mPa_s: f32,
}

impl Default for PlasmaState {
    fn default() -> Self {
        Self {
            po2_mmHg: 100.0,      // Arterial
            pco2_mmHg: 40.0,      // Arterial
            ph: 7.4,              // Normal arterial
            temperature_celsius: 37.0,
            glucose_mM: 5.0,
            osmolarity_mOsm_per_L: 290.0,
            viscosity_mPa_s: 1.2,
        }
    }
}

/// Flow conditions around the cell
#[derive(Debug, Clone)]
pub struct FlowState {
    /// Shear rate (1/s)
    /// Reference: 100-1000 s⁻¹ in arterioles, up to 5000 s⁻¹ in capillaries
    /// Source: Popel & Johnson, Annu Rev Fluid Mech 2005
    pub shear_rate_per_sec: f32,
    /// Flow velocity (μm/s)
    pub flow_velocity_um_per_sec: f32,
    /// Reynolds number (dimensionless)
    /// Reference: Re << 1 for RBCs (Stokes flow regime)
    pub reynolds_number: f32,
    /// Vessel diameter (μm) - 0 for unbounded flow
    pub vessel_diameter_um: f32,
}

impl Default for FlowState {
    fn default() -> Self {
        Self {
            shear_rate_per_sec: 0.0, // Static initially
            flow_velocity_um_per_sec: 0.0,
            reynolds_number: 0.0,
            vessel_diameter_um: 0.0, // Unbounded
        }
    }
}
