#pragma once

#include "flux_decoder.h"

#include <vector>


// Expectation-maximization idea: Let our model fidelity be
//		alpha * sum i = 1..n: (flux_i - step_i * clock_i)^2
// and the smoothing term be
//		(1-alpha) * sum i=2..n: (flux_i - flux_(i-1))^2.

// Then alternate between using acausal EWMA to fit clock_i
// (dividing through by step_i) and fitting step_i (either
// directly as flux_i/clock_i or by rounding).


// Two filter types based on EWMA:
// the usual "filtfilt" two-pass filtering approach,
// and one that averages together forward and backward EWMA.

typedef enum { AFT_FILTFILT, AFT_AVERAGE_FILT } filter_t;

class alternating_EWMA_decoder : public pulse_decoder {
	private:
		double alpha;
		size_t num_iterations;

		filter_t filter_type = AFT_FILTFILT;

		double linear_mix_intercept = 0.5;
		double linear_mix_slope = 0;

		std::vector<double> get_new_clocks_filtfilt(
			const std::vector<int> & fluxes,
			const std::vector<double> & steps) const;

		std::vector<double> get_new_clocks_avg(
			const std::vector<int> & fluxes,
			const std::vector<double> & steps) const;

		std::vector<double> get_new_clocks(
			const std::vector<int> & fluxes,
			const std::vector<double> & steps) const;

		std::vector<double> get_new_steps(
			const std::vector<int> & fluxes,
			const std::vector<double> & clocks,
			int iteration) const;

	public:
		using pulse_decoder::get_MFM_train;

		alternating_EWMA_decoder(double alpha_in, size_t num_iters_in) {
			alpha = alpha_in;
			num_iterations = num_iters_in;
		}

		void set_alpha(double alpha_in) {
			alpha = alpha_in;
		}

		void set_num_iterations(size_t num_iters_in) {
			num_iterations = num_iters_in;
		}

		void set_intercept(double intercept_in) {
			linear_mix_intercept = intercept_in;
		}

		void set_slope(double slope_in) {
			linear_mix_slope = slope_in;
		}

		alternating_EWMA_decoder() : alternating_EWMA_decoder(0.5, 4) {}

		MFM_train_data get_MFM_train(
			const std::vector<int> & fluxes, size_t start_pos,
			size_t end_pos, double & error_out) const;
};
