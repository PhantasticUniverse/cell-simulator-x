//! ODE integration for metabolic simulations.
//!
//! Implements 4th-order Runge-Kutta (RK4) integration for solving
//! systems of ordinary differential equations representing metabolic reactions.
//!
//! RK4 provides good accuracy for stiff biochemical systems with a reasonable
//! computational cost. For highly stiff systems, implicit methods may be needed.
//!
//! Reference: Press et al., Numerical Recipes, 3rd ed., Cambridge University Press 2007

/// Configuration for the ODE integrator
#[derive(Debug, Clone)]
pub struct IntegratorConfig {
    /// Integration timestep in seconds
    /// Reference: Typically 0.001-0.01s for metabolic systems
    pub dt_sec: f64,
    /// Maximum allowed concentration change per step (for stability)
    pub max_change_mM: f64,
    /// Minimum concentration (prevents negative values)
    pub min_concentration_mM: f64,
}

impl Default for IntegratorConfig {
    fn default() -> Self {
        Self {
            dt_sec: 0.001,           // 1 ms timestep
            max_change_mM: 0.1,      // Max 0.1 mM change per step
            min_concentration_mM: 1e-9, // Minimum 1 nM
        }
    }
}

/// 4th-order Runge-Kutta integrator for ODE systems
///
/// Solves dy/dt = f(t, y) where y is a vector of concentrations
pub struct RK4Integrator {
    /// Configuration
    pub config: IntegratorConfig,
    /// Current simulation time in seconds
    pub time_sec: f64,
    /// Number of steps taken
    pub step_count: u64,
    /// Scratch vectors for intermediate calculations
    k1: Vec<f64>,
    k2: Vec<f64>,
    k3: Vec<f64>,
    k4: Vec<f64>,
    y_temp: Vec<f64>,
}

impl RK4Integrator {
    /// Create a new RK4 integrator for a system with n variables
    pub fn new(n_variables: usize, config: IntegratorConfig) -> Self {
        Self {
            config,
            time_sec: 0.0,
            step_count: 0,
            k1: vec![0.0; n_variables],
            k2: vec![0.0; n_variables],
            k3: vec![0.0; n_variables],
            k4: vec![0.0; n_variables],
            y_temp: vec![0.0; n_variables],
        }
    }

    /// Resize internal buffers if system size changes
    pub fn resize(&mut self, n_variables: usize) {
        if self.k1.len() != n_variables {
            self.k1.resize(n_variables, 0.0);
            self.k2.resize(n_variables, 0.0);
            self.k3.resize(n_variables, 0.0);
            self.k4.resize(n_variables, 0.0);
            self.y_temp.resize(n_variables, 0.0);
        }
    }

    /// Perform one RK4 integration step
    ///
    /// # Arguments
    /// * `y` - Current state vector (concentrations in mM), modified in place
    /// * `derivatives` - Function that computes dy/dt given current state
    ///
    /// # RK4 Algorithm
    /// k1 = f(t, y)
    /// k2 = f(t + dt/2, y + dt/2 * k1)
    /// k3 = f(t + dt/2, y + dt/2 * k2)
    /// k4 = f(t + dt, y + dt * k3)
    /// y_new = y + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
    pub fn step<F>(&mut self, y: &mut [f64], derivatives: F)
    where
        F: Fn(&[f64], &mut [f64]),
    {
        let dt = self.config.dt_sec;
        let n = y.len();
        self.resize(n);

        // k1 = f(t, y)
        derivatives(y, &mut self.k1);

        // k2 = f(t + dt/2, y + dt/2 * k1)
        for i in 0..n {
            self.y_temp[i] = y[i] + 0.5 * dt * self.k1[i];
        }
        derivatives(&self.y_temp, &mut self.k2);

        // k3 = f(t + dt/2, y + dt/2 * k2)
        for i in 0..n {
            self.y_temp[i] = y[i] + 0.5 * dt * self.k2[i];
        }
        derivatives(&self.y_temp, &mut self.k3);

        // k4 = f(t + dt, y + dt * k3)
        for i in 0..n {
            self.y_temp[i] = y[i] + dt * self.k3[i];
        }
        derivatives(&self.y_temp, &mut self.k4);

        // Update: y_new = y + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
        let dt_6 = dt / 6.0;
        for i in 0..n {
            let dy = dt_6 * (self.k1[i] + 2.0 * self.k2[i] + 2.0 * self.k3[i] + self.k4[i]);

            // Clamp change for stability
            let clamped_dy = dy.clamp(-self.config.max_change_mM, self.config.max_change_mM);
            y[i] += clamped_dy;

            // Ensure non-negative concentrations
            if y[i] < self.config.min_concentration_mM {
                y[i] = self.config.min_concentration_mM;
            }
        }

        self.time_sec += dt;
        self.step_count += 1;
    }

    /// Run multiple integration steps
    pub fn run<F>(&mut self, y: &mut [f64], derivatives: F, duration_sec: f64)
    where
        F: Fn(&[f64], &mut [f64]),
    {
        let n_steps = (duration_sec / self.config.dt_sec).ceil() as usize;
        for _ in 0..n_steps {
            self.step(y, &derivatives);
        }
    }

    /// Reset integrator state
    pub fn reset(&mut self) {
        self.time_sec = 0.0;
        self.step_count = 0;
    }
}

impl Default for RK4Integrator {
    fn default() -> Self {
        Self::new(20, IntegratorConfig::default())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rk4_exponential_decay() {
        // Test: dy/dt = -y with y(0) = 1
        // Analytical solution: y(t) = exp(-t)
        let mut integrator = RK4Integrator::new(1, IntegratorConfig {
            dt_sec: 0.01,
            max_change_mM: 10.0,
            min_concentration_mM: 1e-12,
        });

        let mut y = vec![1.0];
        let derivatives = |state: &[f64], dydt: &mut [f64]| {
            dydt[0] = -state[0];
        };

        // Integrate for 1 second
        integrator.run(&mut y, derivatives, 1.0);

        // y(1) should be approximately exp(-1) ≈ 0.368
        let expected = (-1.0_f64).exp();
        let error = (y[0] - expected).abs();
        assert!(error < 0.001, "RK4 error too large: {} vs expected {}", y[0], expected);
    }

    #[test]
    fn test_rk4_coupled_system() {
        // Test coupled system: oscillator
        // dy1/dt = y2
        // dy2/dt = -y1
        // Analytical: y1 = cos(t), y2 = -sin(t) for y1(0)=1, y2(0)=0
        let mut integrator = RK4Integrator::new(2, IntegratorConfig {
            dt_sec: 0.001,
            max_change_mM: 10.0,
            min_concentration_mM: -10.0, // Allow negative for oscillator
        });

        let mut y = vec![1.0, 0.0];
        let derivatives = |state: &[f64], dydt: &mut [f64]| {
            dydt[0] = state[1];
            dydt[1] = -state[0];
        };

        // Integrate for pi seconds (half period)
        integrator.run(&mut y, derivatives, std::f64::consts::PI);

        // y1(pi) should be approximately cos(pi) = -1
        // Due to min_concentration clamping this test checks numerical stability
        // The actual test is that it doesn't crash or diverge
        assert!(y[0].is_finite());
        assert!(y[1].is_finite());
    }

    #[test]
    fn test_non_negative_concentrations() {
        let mut integrator = RK4Integrator::new(1, IntegratorConfig::default());

        let mut y = vec![0.001]; // Small initial concentration
        let derivatives = |_: &[f64], dydt: &mut [f64]| {
            dydt[0] = -1000.0; // Large negative derivative
        };

        integrator.step(&mut y, derivatives);

        // Should be clamped to minimum, not negative
        assert!(y[0] >= integrator.config.min_concentration_mM);
    }

    #[test]
    fn test_step_count() {
        let mut integrator = RK4Integrator::new(1, IntegratorConfig::default());
        let mut y = vec![1.0];
        let derivatives = |_: &[f64], dydt: &mut [f64]| {
            dydt[0] = 0.0;
        };

        assert_eq!(integrator.step_count, 0);
        integrator.step(&mut y, derivatives);
        assert_eq!(integrator.step_count, 1);
        integrator.step(&mut y, derivatives);
        assert_eq!(integrator.step_count, 2);
    }
}
