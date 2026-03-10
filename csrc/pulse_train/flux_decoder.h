#pragma once

#include "pulse_train.h"

// Turn a flux record vector into an MFM pulse train by comparing the distance
// between flux reversals. Our job is made more difficult because apparently
// the mapping is nonlinear: one clock length corresponds to 01, 1.5 clock
// lengths to 001, and 2 clock lengths to 0001.

// For this reason, longer delays are undefined. I may handle these later,
// but note that they are out of spec and should never occur on uncorrupted
// floppies with the right clock.

// Now that I'm coming back to the code after two years, it seems
// kinda ugly. Let's try to salvage it. I'll deal with the MFM_train_data
// class later.

// The error_out double is set to a badness-of-fit value: the higher the worse.
// It is used as a quick and dirty clock inference method (see below), but
// is more useful for debugging (e.g. when I implement dewarping, it could
// show how well the signal has been dewarped).

// TODO: Make the error_out a common unit, so that different decoders'
// performance is comparable... but is that even possible? It may be
// better to just remove it... hm. I need to know what it's used for.


// A pulse decoder may have one or more adjustable parameters, such
// as smoothing parameters (for EWMA) or outlier thresholds. The param_t
// type is used in a function that returns what parameters can be set and
// what types they have, so that regions of good performance can be found
// by automated search.

// PARAM_REAL is a real value (double) on a scale from min to max inclusive.
// PARAM_LOG_REAL is the same, but log-spaced (usually this is something
//		multiplicative like a smoothing factor).
// PARAM_INTEGER is an integer from min to max inclusive. Categorical
//		or boolean variables are particular cases.

typedef enum { PARAM_REAL, PARAM_LOG_REAL, PARAM_INTEGER } param_t;

class pulse_decoder {
	public:
		virtual MFM_train_data get_MFM_train(
			const std::vector<int> & fluxes, size_t start_pos,
			size_t end_pos, double & error_out) const = 0;

		MFM_train_data get_MFM_train(
			const std::vector<int> & fluxes,
			double & error_out) const {

			return get_MFM_train(fluxes,
				0, fluxes.size(), error_out);
		}

		// Very wasteful. XXX: Remove later?
		double get_error(const std::vector<int> & fluxes,
			size_t start_pos, size_t end_pos) const {

			double error_out;

			get_MFM_train(fluxes, start_pos, end_pos,
				error_out);

			return error_out;
		}

		virtual ~pulse_decoder() = default;

		// The following functions define the parameter space for
		// the decoder in question, so that we can do grid search
		// without having to know what concrete parameters are
		// available.

		virtual std::vector<param_t> get_parameter_types() const {
			return {};
		}

		virtual std::vector<double> get_parameter_min() const {
			return {};
		}

		virtual std::vector<double> get_parameter_max() const {
			return {};
		}

		// Get the current parameter values.

		virtual std::vector<double> get_current_params() const {
			return {};
		}

		// Set them.

		virtual void set_params(const std::vector<double> & /*params*/) {
			throw std::invalid_argument("flux_decoder: decoder has no parameters");
		}

		virtual std::string name() const = 0;
};