#include <algorithm>
#include <limits>
#include <vector>

#include <iostream>

#include "offset_clock.h"
#include "tools.h"

// This is mostly a cut-and-paste from constant_clock_decoder.
// See classical_decoders.cc for more info.

MFM_train_data offset_clock_decoder::get_MFM_train(
		const std::vector<int> & fluxes, size_t start_pos,
		size_t end_pos, double & error_out) const {

	MFM_train_data train;

	train.data.push_back(1);
	train.flux_indices.push_back(start_pos);

	error_out = 0;

	for (size_t i = start_pos; i < end_pos; ++i) {
		// We'll model the nonlinearity like this:
		//		- There's always half a clock's delay before anything happens.
		//		- Then one zero corresponds to half a clock more, two zeroes is
		//			two halves more, and three zeroes is three, and so on.

		double flux_delay = fluxes[i] - clock_intercept;
		int half_clocks = round((flux_delay*2)/clock_factor);

		// Subtract the constant half-clock offset of one, then round the next
		// to get the number of zeroes.
		// HACK. TODO: Fix.
		int zeroes = std::min(20, std::max(0, half_clocks-1));

		// Update Euclidean error.
		double error_term = flux_delay - half_clocks * clock_factor/2.0;
		error_out += error_term * error_term;

		// Ignore RR -- treat it as noise (as the reference VHDL code
		// does). Better might be to just add the delay to the next
		// term. Maybe do that later.
		if (zeroes == 0) { continue; }

		for (int j = 0; j < zeroes; ++j) {
			train.data.push_back(0);
			train.flux_indices.push_back(i);
		}
		train.data.push_back(1);
		train.flux_indices.push_back(i);
	}

	error_out = std::sqrt(error_out / (end_pos - start_pos));

	return train;
}