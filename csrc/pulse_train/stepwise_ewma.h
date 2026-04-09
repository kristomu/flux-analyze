#pragma once

#include "flux_decoder.h"
#include "ewma.h"

#include "ewma_search/interval.h"

#include <vector>

// A stepwise EWMA decoder: it takes a list of pairs of intervals and
// relative flux vector indices. For each index, the interval states
// that some initial clock value between the interval bounds should set
// at that point.

// It's meant to be used in combination with ewma_search() to directly
// detect preambles and fix the clock in case it got skewed by earlier
// drift or outliers. As such, it's not very useful standalone.

// It doesn't support outlier truncation, but it will clamp output
// half-clocks (number of zero bits) to [1...4], just like ordinary
// EWMA does.

class stepwise_EWMA_decoder : public pulse_decoder {
	private:
		double alpha;
		std::vector<std::pair<int, interval> > initial_clocks;

		double tolerance = 1e-10;

		double max_relative_change = -1;

		bool ignore_lower = false;
		bool clamp = true;

		/*bool are_initial_clocks_set() const {
			return !initial_clock.empty(); }*/

		// decodes each piece
		causal_EWMA_clock_decoder piece_decoder;

	public:
		using pulse_decoder::get_MFM_train;

		// TBD. Parameters is going to be... interesting. We can't
		// grid search on this!

		stepwise_EWMA_decoder(double alpha_in) {
			alpha = alpha_in;
		}

		void set_alpha(double alpha_in) {
			alpha = alpha_in;
		}

		void clear_initial_clocks() {
			initial_clocks.clear();
		}

		void set_ignore_lower(bool ignore_lower_in) {
			ignore_lower = ignore_lower_in;
		}

		void set_clamp(bool clamp_in) {
			clamp = clamp_in;
		}

		void set_max_relative_change(double max_relative_change_in) {
			max_relative_change = max_relative_change_in;
		}

		stepwise_EWMA_decoder() : stepwise_EWMA_decoder(0.5) {
			piece_decoder = causal_EWMA_clock_decoder(0.5, 24);
		}

		MFM_train_data get_MFM_train(
			const std::vector<int> & fluxes, size_t start_pos,
			size_t end_pos, double & error_out) const;

		// Initial clock is not considered a parameter here because
		// we can do without it, and fixing an initial clock wouldn't
		// generalize well. (TODO: Is this the right call?)

		void set_params(const std::vector<double> & params) {
			alpha = params[0];
			ignore_lower = (int)params[1];
			clamp = (int)params[2];
			max_relative_change = params[3];

			piece_decoder.set_params(params);
		}

		std::vector<param_t> get_parameter_types() const {
			return {PARAM_LOG_REAL, PARAM_INTEGER, PARAM_INTEGER, PARAM_LOG_REAL};
		}

		std::vector<double> get_parameter_min() const {
			return {0, 0, 0, 1};
		}

		std::vector<double> get_parameter_max() const {
			return {1, 1, 1, 3};
		}

		std::vector<double> get_current_params() const {
			return {alpha, (double)((int)ignore_lower),
							(double)((int)clamp), max_relative_change};
		}

		std::string name() const {
			return "Stepwise EWMA decoder";
		}
};
