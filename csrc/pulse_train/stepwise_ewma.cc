#include <iostream>
#include <stdexcept>

#include "preambles.h"
#include "stepwise_ewma.h"

#include "ewma_search/ewma_search.h"

// TODO: Compare alpha=0 with ordinal search; where does the
// resulting train differ? They should be the same.

MFM_train_data stepwise_EWMA_decoder::get_MFM_train(
	const std::vector<int> & fluxes, size_t start_pos,
	size_t end_pos, double & error_out) const {

	// Get the A1A1A1 preamble (fix this later)

	/* std::vector<char> short_A1 = {
       	0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1};
		needs to be repeated thrice. In addition, the first
		0 1 should be discarded because we don't know if
		that's a 01, 001, or 0001. Anything will match it.

		So we get the following nonzero sequence:
			3 2 3 2 1
			3 2 3 2 1
			3 2 3 2
		hence the following half-clocks:

	*/

	std::vector<int> half_clock_needle = {
		4, 3, 4, 3, 2,
		4, 3, 4, 3, 2,
		4, 3, 4, 3 };

	error_out = std::numeric_limits<double>::infinity();

	std::vector<int> one_way_truncated_fluxes(
		fluxes.begin() + start_pos,
		fluxes.begin() + end_pos);

	// Look for EWMA matches.

	// Universalizing over alpha would be nice (simulating this by
	// just trying a very large number of candidate alphas makes
	// t78 go from 518 to 525 hits). But the math for doing so is
	// really hairy, getting 1% more matches doesn't sound like
	// it'd be worth the extreme complexity it'd take, and choosing
	// the wrong alpha among every candidate still is sufficient
	// to degrade the t73 performance back to 8 sectors.

	// I think making a variant EWMA dynamic programming approach, that
	// just goes with the value inside the interval that minimizes MFM
	// errors, would provide more bang for the buck, as would doing
	// approximate matching for the preamble (edit distance style).
	
	std::vector<std::pair<int, interval> > matches =
		ewma_search_all(one_way_truncated_fluxes,
			half_clock_needle, alpha, tolerance);

	std::cout << matches.size() << " matches.\n";

	// If no matches were found, hand off to ordinary EWMA.
	if (matches.size() == 0) {
		// TODO: Insert prior of choice here.
		// It shouldn't matter because we didn't find any preambles.
		piece_decoder.set_alpha(0.01);
		piece_decoder.set_initial_clock(22);

		return piece_decoder.get_MFM_train(
			fluxes, start_pos, end_pos, error_out);
	}

	// Pick a "typical" value for the clock so we don't wiggle it
	// about too much. This makes sense because an uncorrupted
	// floppy should have the same clock everywhere.

	double total_clock = 0;

	size_t i;

	for (i = 0; i < matches.size(); ++i) {
		total_clock += (matches[i].second.lower_clock + 
			matches[i].second.upper_clock) / 2.0;
	}

	double mean_clock = total_clock / (double)matches.size();

	// Stitch together the pieces, each using the internal EWMA decoder.

	piece_decoder.set_initial_clock(mean_clock);
	piece_decoder.set_alpha(alpha);

	MFM_train_data composite_train = piece_decoder.get_MFM_train(
		fluxes, start_pos, matches[0].first, error_out);

	for (i = 0; i < matches.size(); ++i) {
		size_t next_start_pos = end_pos;

		// If there's a next match, then our ending position is
		// where it starts.
		if (i+1 < matches.size()) {
			next_start_pos = matches[i+1].first;
		}

		double initial_clock = mean_clock;

		if (initial_clock < matches[i].second.lower_clock) {
			initial_clock = matches[i].second.lower_clock;
		}

		if (initial_clock > matches[i].second.upper_clock) {
			initial_clock = matches[i].second.upper_clock;
		}

		piece_decoder.set_initial_clock(initial_clock);
		piece_decoder.set_alpha(alpha);

		MFM_train_data this_piece = piece_decoder.get_MFM_train(fluxes,
			matches[i].first, next_start_pos, error_out);

		composite_train += this_piece;
	}

	return composite_train;
}
