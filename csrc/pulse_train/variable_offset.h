#pragma once

#include "flux_decoder.h"

#include <vector>

// Decoders that estimate y(t) = a(t) x(t) + b(t)
// with y(t) being the number of half-clocks (number of zeroes plus one),
// x(t) being the observed flux delay,
// and a(t) and b(t) being slowly varying coefficients.

// It's like EWMA but with a 2D state space.

class GD_variable_offset_decoder : public pulse_decoder {
	private:
		double stepsize;

	public:
		using pulse_decoder::get_MFM_train;

		GD_variable_offset_decoder(double stepsize_in) {
			stepsize = stepsize_in;
		}

		void set_stepsize(double stepsize_in) {
			stepsize = stepsize_in;
		}

		GD_variable_offset_decoder() : GD_variable_offset_decoder(0.5) {}

		MFM_train_data get_MFM_train(
			const std::vector<int> & fluxes, size_t start_pos,
			size_t end_pos, double & error_out) const;

		std::vector<param_t> get_parameter_types() const {
			return {PARAM_LOG_REAL};
		}

		std::vector<double> get_parameter_min() const {
			return {0};
		}

		std::vector<double> get_parameter_max() const {
			return {0.003};
		}

		std::vector<double> get_current_params() const {
			return {stepsize};
		}

		void set_params(const std::vector<double> & params) {
			stepsize = params[0];
		}

		std::string name() const {
			return "Gradient descent-type variable offset decoder";
		}
};