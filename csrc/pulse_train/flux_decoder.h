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
};