#include <algorithm>
#include <limits>
#include <vector>

#include <iostream>

#include "classical_decoders.h"
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

// TODO: It might be useful to create an "explicit clock vector decoder" that
// takes a list of clock values and creates an MFM_train_data object based on
// it. Then constant-clock is trivial and the "modern" decoders are fairly
// easy too. Only problem is that it's creating and destroying a bunch of
// temporary vectors when not doing that would be faster. Optimize later,
// perhaps.

// Level one:

MFM_train_data constant_clock_decoder::get_MFM_train(
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

	error_out = std::sqrt(error_out / (end_pos - start_pos));

	return train;
}

// Horribly cut and paste code, but it'll work as a proof of concept.
// MFM decoder with dewarping based on exponential smoothing.
// Better (and less cut and paste) stuff later.
// The exponential smoothing has no right being this powerful, either
// (not that I'm really complaining).

// If different friction or warping of the surface were to bias the
// clock, we'd expect the bias to be proportional to the length between
// clock transitions as a longer delay would mean more time to accumulate.
// The dewarp function below instead assumes that the bias is independent
// of length, and bizarrely, that seems to work *better*, at least for
// certain chunks that are otherwise very tough to decode. I don't get it,
// but I can't dispute the results.

MFM_train_data orig_EWMA_causal_clock_decoder::get_MFM_train(
	const std::vector<int> & fluxes, size_t start_pos,
	size_t end_pos, double & RMSE_out) const {

	MFM_train_data train;

	train.data.push_back(1);
	train.flux_indices.push_back(0);

	RMSE_out = 0;
	double bias = 0;

	double clock = initial_clock;

	for (size_t i = start_pos; i < end_pos; ++i) {
		int corrected_flux_delay = fluxes[i] - bias;
		int half_clocks = round((corrected_flux_delay*2)/clock);

		// The error term is the difference between the
		// observed (corrected) flux delay and what we would
		// have seen if the delay had been an exact multiple of
		// the initial clock rate.

		double error_term = corrected_flux_delay - half_clocks * clock/2.0;
		RMSE_out += error_term * error_term;

		bias = bias * (1-alpha) + error_term * alpha;

		int zeroes = std::max(0, half_clocks-1);

		if (zeroes == 0) { continue; }

		for (int j = 0; j < zeroes; ++j) {
			train.data.push_back(0);
			train.flux_indices.push_back(i);
		}
		train.data.push_back(1);
		train.flux_indices.push_back(i);
	}

	RMSE_out = std::sqrt(RMSE_out / (end_pos - start_pos));

	return train;
}

// Trying a reasonably faithful reproduction of the VHDL PLL code.

// Changing ints to doubles makes it worse.

MFM_train_data historical_EWMA_decoder::get_MFM_train(
		const std::vector<int> & fluxes, size_t start_pos,
		size_t end_pos, double & RMSE_out) const {

	MFM_train_data train;
	
	train.data.push_back(1);
	train.flux_indices.push_back(0);

	RMSE_out = std::numeric_limits<double>::infinity();

	train.data.push_back(1);
	train.flux_indices.push_back(start_pos);
 
	int old_clock, new_clock;
	int clock = fluxes[start_pos]; // initial clock estimate ???

	for (size_t i = start_pos; i < end_pos; ++i) {

		// pulse edge, store clock in a variable
		old_clock = clock;

		int length = fluxes[i]; // is this what length is ???

	    // Not entirely sure if this is what's going on
	    // "skip" is just outputting a zero, I think

		bool emit = false;
		int mfm_type = -1;
		new_clock = clock;

		//          0 ...  1 clock/4      -> too short
		//  1 clock/4 ...  5 clock/4      -> short
		//  5 clock/4 ...  7 clock/4      -> medium
		//  7 clock/4 ... 12 clock/4      -> long
		// 12 clock/4 ... infty	          -> too long
		
		if (4 * length < clock) { // too short
			// nothing
		} else if (4 * length < 5 * old_clock) { // 1 us
			new_clock = old_clock - (old_clock - length) / 2;
			mfm_type = 2;
			emit = true;
		} else if (4 * length < 7 * old_clock) { // 1.5 us
			emit = true;
			mfm_type = 3;
			// Some kind of magic; the whole thing breaks down
			// if you change the ints to doubles, for some reason.

			// new_clock = old_clock/2 + 11/32 length
			// but not quite.

			new_clock = old_clock - (old_clock
				- length/2
				- length/4
				+ length/16) / 2;
		} else if (length < 3 * clock) { // 2 us
			emit = true;
			mfm_type = 4;
			new_clock = old_clock - (old_clock - length/2)/2;
		} else { // Too long a duration
			// skip = 5; // ???
		}

		if (emit) {
			for (int j = 0; j < mfm_type-1; ++j) {
				train.data.push_back(0);
				train.flux_indices.push_back(start_pos + i);
			}
			train.data.push_back(1);
			train.flux_indices.push_back(start_pos + i);
		}

		// Update clock
		// The maximum amount that the clock can change is
		// bounded.

		if (new_clock > 9 * old_clock / 8) {
			clock = 9 * old_clock / 8;
		} else if ( new_clock < 7 * old_clock / 8) {
			clock = 7 * old_clock / 8;
		} else {
			clock = new_clock;
		}
	}

	return train;
}