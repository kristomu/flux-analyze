#pragma once

#include "flux_decoder.h"

class offset_clock_decoder: public pulse_decoder {
	private:
		double clock_factor, clock_intercept;

	public:
		// Do I want to do this everywhere, or just do a lazy
		// two-step function approach instead? Hmm...
		using pulse_decoder::get_MFM_train;

		void set_clock(double clock_in) { clock_factor = clock_in; }

		void set_intercept(double clock_intercept_in) {
			clock_intercept = clock_intercept_in;
		}

		MFM_train_data get_MFM_train(
			const std::vector<int> & fluxes, size_t start_pos,
			size_t end_pos, double & error_out) const;

		// clock is a parameter of sorts, though an ill-fitting one.
		// Might as well expose it here.
		std::vector<param_t> get_parameter_types() const {
			return {PARAM_REAL, PARAM_REAL};
		}

		std::vector<double> get_parameter_min() const {
			return {-10, 1};
		}

		std::vector<double> get_parameter_max() const {
			return {10, 50}; // For reasonable floppy types?
		}

		std::vector<double> get_current_params() const {
			return {clock_intercept, clock_factor};
		}

		void set_params(const std::vector<double> & params) {
			clock_intercept = params[0];
			clock_factor = params[1];
		}

		std::string name() const {
			return "Offset constant-clock decoder";
		}
};
