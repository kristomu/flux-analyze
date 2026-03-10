#include "pulse_train/flux_decoder.h"
#include "tools.h"

// This is a (relatively simple) dynamic programming-based flux/pulse
// decoder. Let y_i be a guess for the clock value at the ith flux
// transition. Let m_i be the number of zero bits (MFM decoding) given
// x_i and the guessed clock y_i.

// It then seeks to find the clock values y_i that minimize

//		sum i=1...n alpha * (x_i - f(y_i, m_i))^2 +
//		sum i=2...n (1-alpha) * (y_i - y_(i-1))^2 +
//		beta * M(c)

// where f is the function that determines the flux transition value
// we would expect to see if m_i were exactly centered with a clock of
// c_i;

// and M(c) is a count of the number of MFM errors observed during
// the given decoding c_1,...c_n.

// The decoder makes no attempt to manually check for preambles with
// intentional errors, nor of inserting missing bits.

// There's also no memory check. Try this on large sequences at your
// own risk.

// I may need to tweak the MFM error term later so that it's independent
// of the scale of the flux transitions (i.e. treating DD and HD floppies
// equally).

// -----

// This is very rough and also very slow. I should reimplement it properly
// later, or use simpler DP approaches. Use with caution.

class QND_one_mode_DP : public pulse_decoder {
	private:
		double alpha = 0.005;
		double beta = 50.0;

		// Number of flux transitions at the start where we
		// give no penalty for MFM errors. Hack for testing this
		// in combination with ordinal search.
		size_t grace_period = 0;

		// If true, only runs a single chunk (for testing).
		bool do_only_one = true;

		size_t fidelity_err_power = 2;
		size_t smoothing_err_power = 2;

		bool is_mfm_valid(short zero_bits, u_char exclusive_parity) const;

	public:
		using pulse_decoder::get_MFM_train;

		MFM_train_data get_MFM_train(
			const std::vector<int> & fluxes, size_t start_pos,
			size_t end_pos, double & error_out) const;

		std::string name() const {
			return "One-mode DP";
		}

		void set_grace_period(size_t grace_period_in) {
			grace_period = grace_period_in;
		}

		void set_beta(double beta_in) { beta = beta_in; }

		void set_alpha(double alpha_in) { alpha = alpha_in; }

		void set_fidelity_err_power(size_t fidelity_err_power_in) {
			fidelity_err_power = fidelity_err_power_in;
		}

		void set_smoothing_err_power(size_t smoothing_err_power_in) {
			smoothing_err_power = smoothing_err_power_in;
		}


		void set_params(const std::vector<double> & params) {
			alpha = params[0];
			beta = params[1];
			fidelity_err_power = params[2];
			smoothing_err_power = params[3];
		}

		std::vector<param_t> get_parameter_types() const {
			return {PARAM_LOG_REAL, PARAM_LOG_REAL, PARAM_INTEGER, PARAM_INTEGER};
		}

		std::vector<double> get_parameter_min() const {
			return {0, 0, 1, 1};
		}

		std::vector<double> get_parameter_max() const {
			return {1, 100, 10, 10};
		}

		std::vector<double> get_current_params() const {
			return {alpha, beta, (double)fidelity_err_power,
				(double)smoothing_err_power};
		}
};
