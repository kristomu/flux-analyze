#include <algorithm>
#include <limits>
#include <vector>
#include <cmath>

#include <iostream>

#include "pulse_train.h"

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

MFM_train_t get_MFM_train(double clock,
		const std::vector<int> & fluxes, double & error_out) {

	MFM_train_t sequence_bits;
	// expected number of bits per reversal, including ones.
	sequence_bits.reserve(fluxes.size() * 8/3);
	sequence_bits.push_back(1);

	error_out = 0;

	for (auto flux_delay: fluxes) {
		// We'll model the nonlinearity like this:
		//		- There's always half a clock's delay before anything happens.
		//		- Then one zero corresponds to half a clock more, two zeroes is
		//			two halves more, and three zeroes is three, and so on.

		double half_clocks = (flux_delay*2)/clock;

		// Subtract the constant half-clock offset of one, then round the next
		// to get the number of zeroes.
		int zeroes = std::max(0, (int)std::round(half_clocks - 1));

		// Update Euclidean error.
		// Clamp the zeroes to [1..3] to penalize wrong clock guesses.
		// NOTE: This may produce very misleading errors if the disk is
		// partially corrupted because the mean is not a robust estimator.
		// Deal with it later if required.
		int clamped_zeroes = std::max(1, std::min(3, zeroes));
		double error_term = flux_delay - (clamped_zeroes + 1) * clock/2.0;
		error_out += error_term * error_term;

		for (int j = 0; j < zeroes; ++j) {
			sequence_bits.push_back(0);
		}
		sequence_bits.push_back(1);
	}

	error_out = std::sqrt(error_out / fluxes.size());

	return sequence_bits;
}

double get_MFM_train_error(double clock, const std::vector<int> & fluxes) {

	double error = 0;
	get_MFM_train(clock, fluxes, error);
	return error;
}

// Baseline inference for clocks. This should be very quick and work for
// most non-corrupted floppies without warping. The way this works is that
// it takes the median of three flux delays, then checks the best error it
// can get by assuming that delay roughly corresponds to one, 1.5,
// or 2 clocks.

template<typename T> double median(std::vector<T> vec) {
	// This will only be called for small vectors, so just sort.
	std::sort(vec.begin(), vec.end());

	if (vec.size() % 2 == 0) {
		return 0.5 * (vec[vec.size()/2-1] + vec[vec.size()/2]);
	} else {
		return vec[vec.size()/2];
	}
}


// Ternary search. Fibonacci would be more efficient but also
// a lot more code.
double get_local_clock_optimum(double left, double right,
	const std::vector<int> & flux_sample) {

	if (fabs(right - left) < 0.01) {
		return (right + left) * 0.5;
	}

	double center_left = (2 * left + right) / 3.0;
	double center_right = (2 * right + left) / 3.0;

	if (get_MFM_train_error(center_left, flux_sample) <
		get_MFM_train_error(center_right, flux_sample)) {

		return get_local_clock_optimum(left, center_right,
			flux_sample);
	} else {
		return get_local_clock_optimum(center_left, right,
			flux_sample);
	}
}

double infer_clock(std::vector<int>::const_iterator fluxes_start,
	std::vector<int>::const_iterator fluxes_end) {

	// First find some typical delay by sampling evenly.
	size_t f = fluxes_end - fluxes_start;
	double typical_delay = median(std::vector<int>(
		{fluxes_start[0], fluxes_start[f/3], fluxes_start[2*f/3]}));

	// Assuming a normal MFM encoding, this delay could either be
	// one clock, 1.5 or 2 clocks. We'll test each hypothesis by
	// finding the local optimal clock consistent with it; for that
	// we need some intervals to search within.

	// Create a representation of these intervals. Add some large
	// value so that (guess_intervals[i] ... [i+1]) form intervals
	// that together cover every plausible clock value (except when i=0,
	// which we special-case below).
	std::vector<double> guess_interval = {typical_delay,
		typical_delay * 1.5, typical_delay * 2, typical_delay * 3};

	// Take an even spread of 100 points to make error calculation
	// much quicker. XXX: It might be possible to make this very quick
	// by using a histogram, try that later at some point???
	size_t i, sample_size = 100, data_size = f;
	std::vector<int> flux_sample;
	if (data_size < sample_size) {
		flux_sample = std::vector<int>(fluxes_start, fluxes_end);
	} else {
		for (i = 0; i < data_size; i += data_size/sample_size) {
			flux_sample.push_back(fluxes_start[i]);
		}
	}

	double record_error = std::numeric_limits<double>::infinity();
	double record_clock = 0;

	for (i = 0; i < guess_interval.size() - 1; ++i) {
		double left = guess_interval[i],
			right = guess_interval[i+1];

		// The first interval descends all the way to zero.
		if (i == 0) { left = 0; }

		double guess = get_local_clock_optimum(left, right, flux_sample);
		double error = get_MFM_train_error(guess, flux_sample);
		if (error < record_error) {
			record_error = error;
			record_clock = guess;
		}
	}

	return record_clock;
}

double infer_clock(const std::vector<int> & fluxes) {
	return infer_clock(fluxes.begin(), fluxes.end());
}
