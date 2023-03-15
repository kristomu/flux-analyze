#include <algorithm>
#include <limits>
#include <vector>

#include <iostream>

#include "pulse_train.h"
#include "tools.h"

// Turn a flux record vector into an MFM pulse train by comparing the distance
// between flux reversals. Our job is made more difficult because apparently
// the mapping is nonlinear: one clock length corresponds to 01, 1.5 clock
// lengths to 001, and 2 clock lengths to 0001.

// For this reason, longer delays are undefined. I may handle these later,
// but note that they are out of spec and should never occur on uncorrupted
// floppies with the right clock.

// The error_out double is set to a badness-of-fit value: the higher the worse.
// It is used as a quick and dirty clock inference method (see below), but
// is more useful for debugging (e.g. when I implement dewarping, it could
// show how well the signal has been dewarped).

// In addition, the nonlinearity might save us sometimes: suppose that we have an
// interval of 2.5 clocks. Then it's likely that this is either one clock followed
// by 1.5 clocks (01001) or the other way around (00101), with some central flux
// reversal having been erased. But I'm not going to do recovery before I've got
// this working on normal images.

// Level one:



MFM_train_data get_MFM_train(double clock,
		const std::vector<int> & fluxes, size_t start_pos,
		size_t end_pos, double & error_out) {

	MFM_train_data train;
	// expected number of bits per reversal, including ones.
	train.data.reserve(fluxes.size() * 8/3);
	train.flux_indices.reserve(fluxes.size() * 8/3);

	train.data.push_back(1);
	train.flux_indices.push_back(start_pos);

	error_out = 0;

	for (size_t i = start_pos; i < end_pos; ++i) {
		// We'll model the nonlinearity like this:
		//		- There's always half a clock's delay before anything happens.
		//		- Then one zero corresponds to half a clock more, two zeroes is
		//			two halves more, and three zeroes is three, and so on.

		int flux_delay = fluxes[i];
		int half_clocks = round((flux_delay*2)/clock);

		// Subtract the constant half-clock offset of one, then round the next
		// to get the number of zeroes.
		int zeroes = std::max(0, half_clocks-1);

		// Update Euclidean error.
		double error_term = flux_delay - half_clocks * clock/2.0;
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

	error_out = std::sqrt(error_out / fluxes.size());

	return train;
}

MFM_train_data get_MFM_train(double clock,
		const std::vector<int> & fluxes, double & error_out) {

	return get_MFM_train(clock, fluxes, 0, fluxes.size(), error_out);
}

double get_MFM_train_error(double clock, const std::vector<int> & fluxes) {

	double error = 0;
	get_MFM_train(clock, fluxes, error);
	return error;
}

// Horribly cut and paste code, but it'll work as a proof of concept.
// MFM decoder with dewarping based on exponential smoothing.
// Better (and less cut and paste) stuff later.
// The exponential smoothing has no right being this powerful, either
// (not that I'm really complaining).

MFM_train_data get_MFM_train_dewarp(double clock,
		const std::vector<int> & fluxes, double & error_out) {

	MFM_train_data train;
	train.data.reserve(fluxes.size() * 8/3);
	train.flux_indices.reserve(fluxes.size() * 8/3);

	train.data.push_back(1);
	train.flux_indices.push_back(0);

	error_out = 0;
	double bias = 0;

	for (size_t i = 0; i < fluxes.size(); ++i) {
		int flux_delay = fluxes[i] - bias;
		int half_clocks = round((flux_delay*2)/clock);

		int zeroes = std::max(0, half_clocks-1);

		double error_term = flux_delay - half_clocks * clock/2.0;
		error_out += error_term * error_term;

		double alpha = 0.03;
                bias = bias * (1-alpha) + error_term * alpha;

		if (zeroes == 0) { continue; }

		for (int j = 0; j < zeroes; ++j) {
			train.data.push_back(0);
			train.flux_indices.push_back(i);
		}
		train.data.push_back(1);
		train.flux_indices.push_back(i);
	}

	error_out = std::sqrt(error_out / fluxes.size());

	return train;
}
