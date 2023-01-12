#include <algorithm>
#include <limits>
#include <vector>
#include <cmath>

#include <iostream>

// Turn a flux record vector into an MFM pulse train by comparing the distance
// between flux reversals. Our job is made more difficult because apparently
// the mapping is nonlinear: one clock length corresponds to 01, 1.5 clock
// lengths to 001, and 2 clock lengths to 0001.

// For this reason, longer delays are undefined. I may handle these later,
// but note that they are out of spec and should never occur on uncorrupted
// floppies with the right clock.

// The error_out double is set to a badness-of-fit value: the higher the worse.
// This might be usable as a very quick clock inference method, but I'll have to
// test more before I know for sure.

// In addition, the nonlinearity might save us sometimes: suppose that we have an
// interval of 2.5 clocks. Then it's likely that this is either one clock followed
// by 1.5 clocks (01001) or the other way around (00101), with some central flux
// reversal having been erased. But I'm not going to do recovery before I've got
// this working on normal images.

// Level one:

std::vector<char> get_MFM_train(double clock,
		const std::vector<int> & fluxes, double & error_out) {

	std::vector<char> sequence_bits;
	sequence_bits.reserve(fluxes.size() * 4 / 3); // expected number of bits per reversal
	sequence_bits.push_back(1);

	error_out = 0;

	for (size_t i = 0; i < fluxes.size(); ++i) {
		// We'll model the nonlinearity like this:
		//		- There's always half a clock's delay before anything happens.
		//		- Then one zero corresponds to half a clock more, two zeroes is
		//			two halves more, and three zeroes is three, and so on.

		double half_clocks = (fluxes[i]*2)/clock;

		// Subtract the constant half-clock offset of one, then round the next
		// to get the number of zeroes.
		int zeroes = std::max(0, (int)std::round(half_clocks - 1));

		// Update Euclidean error.
		// Clamp the zeroes to [1..3] to penalize wrong clock guesses.
		// NOTE: This may produce very misleading errors if the disk is
		// partially corrupted because the mean is not a robust estimator.
		// Deal with it later if required.
		int clamped_zeroes = std::max(1, std::min(3, zeroes));
		double error_term = fluxes[i] - (clamped_zeroes + 1) * clock/2.0;
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
// most non-corrupted floppies. The way this works is that it takes the
// median of three flux delays, then checks the best error it can gets
// by assuming that delay roughly corresponds to one, 1.5, or 2 clocks.

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

double infer_clock(const std::vector<int> & fluxes) {
	// First find some typical delay by sampling evenly.
	size_t f = fluxes.size();
	double typical_delay = median(std::vector<int>({fluxes[0],
		fluxes[f/3], fluxes[2*f/3]}));

	// Create clock estimates based on this delay being exactly one
	// clock, 1.5 or 2 clocks. Add some large value so that
	// guess_intervals[i] ... [i+1] form intervals that together cover
	// every plausible clock value (except when i=0, which we
	// special-case below).
	std::vector<double> guess_interval = {typical_delay,
		typical_delay * 1.5, typical_delay * 2, typical_delay * 3};

	// Take an even spread of 100 points to make error calculation
	// much quicker.
	size_t i, sample_size = 100;
	std::vector<int> flux_sample;
	if (fluxes.size() < sample_size) {
		flux_sample = fluxes;
	} else {
		for (i = 0; i < fluxes.size(); i += fluxes.size()/sample_size) {
			flux_sample.push_back(fluxes[i]);
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