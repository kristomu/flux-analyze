// Do an ordinal search to more robustly find the start of
// floppy data chunks ("address marks"). The ordinal search
// makes use of the observation that if the magic sequence
// that starts a preamble is, say 1 0 0 1 0 0 0 1, then the
// delay between the first and second flux transition must
// be less than the delay between the second and third (absent
// severe warping).

// The ordinal pattern is just a guideline. For instance, the
// flux pattern above (an increase in delay) might arise from
// 1 0 1 0 0 1. Thus, after having found ordinal candidates, we
// need to try to fit a clock at that position. If that's
// impossible, then it's a false positive; otherwise, we get a
// good idea of what the clock actually *is* at that position.

#include "ordinal_search.h"
#include "rabin_karp.h"
#include "preambles.h"

#include <stdexcept>
#include <algorithm>

#include <iostream>
#include <cmath>

// Gets a run length encoding of zeroes broken up by ones, e.g.
// 00101 becomes "2 1". The start must be a zero and there must be
// a one before the end.
std::vector<char> get_run_length_coding(
	std::vector<char>::const_iterator start,
	std::vector<char>::const_iterator end) {

	char run_length = 0;
	std::vector<char> run_length_sequence;

	for (auto pos = start; pos != end; ++pos) {
		if (*pos == 1) {
			run_length_sequence.push_back(run_length);
			run_length = 0;
		} else {
			++run_length;
		}
	}

	return run_length_sequence;
}

std::vector<char> get_ordinal_search_sequence(
	const std::vector<char> & MFM_train_search_sequence) {

	// The sequence must begin in a one, otherwise we don't
	// know the length of the first flux delay (e.g. the
	// search sequence is "0 1 0 0 1" and we don't know if it's
	// really "1 0 1 0 0 1" or " 1 0 0 0 1 0 0 1").

	if (MFM_train_search_sequence[0] != 1) {
		throw std::logic_error("Can't make ordinal sequence from "
			"pattern starting in 0");
	}

	// 1 means this entry must be greater (longer delay) than
	// the previous one. -1 means shorter. 0 is not allowed yet,
	// as there's no way we can match it with the noise in flux timings.
	std::vector<char> ordinal_sequence;

	// First create the run-length encoding of the search sequence
	// (e.g. 10010001 becomes 2 3). Then get the delta values: these
	// form our search sequence. Skip the initial one.

	std::vector<char> run_length_sequence = get_run_length_coding(
		MFM_train_search_sequence.begin() + 1,
		MFM_train_search_sequence.end());

	std::vector<char> ordinal_search_sequence = get_delta_coding(
		run_length_sequence);

	for (char digit: ordinal_search_sequence) {
		if (digit == 0) {
			throw std::logic_error("Can't make ordinal sequence with "
				"more than one of the same clock length in a row!");
		}
	}

	return ordinal_search_sequence;
}

// Given a match from an ordinal search, check if it's really a match to
// the MFM train search sequence. The logic is based on the observation
// that we can make bins for each clock "band": e.g. if we have a search
// sequence of 1001010001, then the first matching flux delay must be in
// the second band, the second must be in the first, and the third must
// be in the third band. If the ranges of these bands overlap, then it's
// impossible to assign a clock that will put each value in its proper
// band, and we have a false positive.

// For instance, if we're ordinally matching 100101000101 (i.e. flux delay
// gets shorter, then longer, then shorter), and we have a match where the
// flux delays are 30 10 20 11, then that would assign the values 10 and
// 11 to band 1, the value 30 to band 2, and the value 20 to band 3. This
// is a false positive because band 3 can't contain a value
// that's shorter than a value in band 2.

// Note that we can't guarantee that even if the bands are well-separated,
// there's a clock value that fits. The Python version works directly on
// band centers and thus avoids this problem, but I'll see if it's a real
// problem before I change anything.

template<typename T> int max_vec(const T & x) {
	return *std::max_element(x.begin(), x.end());
}

template<typename T> int min_vec(const T & x) {
	return *std::min_element(x.begin(), x.end());
}

double get_clock(const std::vector<char> & MFM_train_search_sequence,
	std::vector<int>::const_iterator match_start,
	std::vector<int>::const_iterator fluxes_end) {

	// A leading run of zeroes could technically be used, but it would be
	// too complex and not worth it.
	if (MFM_train_search_sequence[0] != 1) {
		throw std::logic_error("Can't get clock from fluxes with "
			"pattern starting in 0");
	}

	// First create the run length sequence to get the bands that
	// fluxes should be assigned to.
	std::vector<char> run_length_sequence = get_run_length_coding(
		MFM_train_search_sequence.begin() + 1,
		MFM_train_search_sequence.end());

	int num_bands = max_vec(run_length_sequence) + 1;

	std::vector<std::vector<int> > bands(num_bands);

	// Do a very ugly impression of the zip operator.
	auto match_pos = match_start;
	size_t i;

	for (i = 0; i < run_length_sequence.size() &&
		match_pos != fluxes_end; ++i) {
		// Assign the flux at match_pos to the band it would be
		// associated with if the match is correct.
		bands[run_length_sequence[i]].push_back(*match_pos++);
	}

	// Check for monotonicity (i.e. no interval overlap)
	// Check that no band's minimum value is lower or equal to the
	// previous band's max value.

	std::vector<double> min_in_band(bands.size(), 0),
		max_in_band(bands.size(), 0);

	for (i = 1; i < bands.size(); ++i) {
		// Get the minima and maxima for each band. This is
		// used for clock determination.
		max_in_band[i] = max_vec(bands[i]);
		min_in_band[i] = min_vec(bands[i]);

		if (bands[i].empty() || bands[i-1].empty()) {
			continue;
		}

		if (min_in_band[i] <= max_in_band[i-1]) {
			return -1;
		}
	}

	// Clock determination: the clock must be set so that the
	// minimum and maximum value for are assigned to that band.
	double clock_lower, clock_upper;

	clock_upper = std::numeric_limits<double>::max();
	clock_lower = std::numeric_limits<double>::min();

	// Find the upper and lower bounds for the clock.
	for (size_t band = 1; band < bands.size(); ++band) {
		// The highest clock we can have is so that
		// round(min_in_band / half_clock) = band + 1, i.e.
		// min_in_band / half_clock = band + 0.5 + epsilon,
		// clock = 4 * min_in_band / (2 * band + 1)
		clock_upper = std::min(clock_upper, 4 * min_in_band[band] / (2.0 * band + 1));
		// Similarly, the lowest clock we can have is so that
		// round(max_in_band / half_clock) = band + 1,
		// max_in_band / half_clock = band + 1.5 - epsilon,
		// clock = 4 * max_in_band / (2 * band + 3)
		clock_lower = std::max(clock_lower, 4 * max_in_band[band] / (2.0 * band + 3));
	}

	// Because neither clock_lower nor clock_upper will work, we
	// must require that there's some value in between, so weak
	// inequality wouldn't work.
	if (clock_lower < clock_upper) {
		// Empirical weighting, seems to work well.
		return 0.6 * clock_lower + 0.4 * clock_upper;
	} else {
		// There are distinct bands but it's not possible to fit
		// a single clock to it. This needs Python-style manually
		// defined bands. For now, signal the impossibility.
		return -2;
	}
}

// TODO? Somehow generalize so that it doesn't just take preambles.
// E.g. when we start looking for damaged IDAMs, we'll want to get
// results that are given by the IBM_preamble structure.

// Also TODO: using the offset to offset the preamble so that it matches
// the offset we gave the ordinal preamble... is kind of ugly.

std::vector<match_with_clock> filter_matches(
	const std::vector<int> & flux_transitions,
	const std::vector<search_result> & possible_matches,
	const IBM_preamble & preamble_info) {

	std::vector<match_with_clock> out;
	match_with_clock true_match;

	for (const search_result & result: possible_matches) {
		double clock =	get_clock(
			preamble_info.get_offset_preamble_by_ID(result.ID),
			flux_transitions.begin() + result.idx, flux_transitions.end());

		if (clock < 0) {
			continue; // false positive or impossible to fit clock
		}

		true_match.match_location = result.idx;
		true_match.estimated_clock = clock;

		out.push_back(true_match);
	}

	return out;
}