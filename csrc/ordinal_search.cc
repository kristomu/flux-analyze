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
#include "preambles.h"

#include <stdexcept>
#include <algorithm>

#include <iostream>
#include <cmath>

// The ordinal search process consists of two phases: In the first phase,
// we do delta coding (e.g. 100101 becomes "second is shorter than first").
// In the second phase, we try to fit a clock value to the result, thus
// eliminating certain false positives and giving the MFM train decoder
// the information it needs to try to decode the rest of the sector.

// When fitting a clock to a particular ordinal match, there are two types
// of constraint. The former happens at the edges and we'll call them
// "at least" constraints, while the latter happens everywhere else.

// For instance, consider that the MFM train to be matched is
// 0010100010100. Then since we don't know how long the first
// MFM run is - it could be 1001 or 10001 - we only know that the
// first pulse delay must, when decoded, produce two or more zeroes.
// Similarly, for the last MFM run, we match 100, which could be
// either 1001 or 10001, so the situation is the same. For every
// other MFM run (101, 10001, and 101), we know exactly how many
// zeroes a correct clock should produce.
class bound {
	public:
		int run;
		bool at_least; // if true, at-least constraint; if false, exact.

		bound(int run_in, bool al_in) {
			run = run_in;
			at_least = al_in;
		}

		bool operator<(const bound & rhs) { return run < rhs.run; }

		operator int() const { return run; }
};

// Gets a run length encoding of zeroes broken up by ones, e.g.
// 00101 becomes "2 1".
std::vector<bound> get_run_length_coding(
	std::vector<char>::const_iterator start,
	std::vector<char>::const_iterator end) {

	int run_length = 0;
	std::vector<bound> run_length_sequence;

	// the first run of zeroes is an at-least constraint
	// because we don't know how long it actually is. (Making
	// assumptions about the longest possible run would
	// complicate things unnecessarily.)
	bool at_least = true;

	for (auto pos = start; pos != end; ++pos) {
		if (*pos == 1) {
			run_length_sequence.push_back(
				bound(run_length, at_least));
			run_length = 0;
			at_least = false;
		} else {
			++run_length;
		}
	}

	// The last run length is also an at-least constraint,
	// if there was any after the last one.
	if (run_length > 0) {
		run_length_sequence.push_back(
			bound(run_length, true));
	}

	return run_length_sequence;
}

// Returns the ordinal search sequence as well as an offset
// indicating how far into the actual thing we're searching for
// that an ordinal match would be.
ordinal_pattern get_ordinal_search_sequence(
	const std::vector<char> & MFM_train_search_sequence) {

	// The format of the needle is: 1 means this entry must be greater
	// (longer delay) than the previous one. -1 means shorter. 0 is not
	// allowed yet, as there's no way we can match it with the noise in
	// flux timings -- it would need a wildcard search.
	ordinal_pattern out;

	// First create the run-length encoding of the search sequence
	// (e.g. 10010001 becomes 2 3). Then get the delta values: these
	// form our search sequence.

	std::vector<bound> run_length_sequence = get_run_length_coding(
		MFM_train_search_sequence.begin(),
		MFM_train_search_sequence.end());

	// Skip past every at_least value, because we don't know what
	// its value will be, we can't subtract from it either.

	for (out.offset = 0; out.offset < run_length_sequence.size() &&
		run_length_sequence[out.offset].at_least; ++out.offset) {}

	if (out.offset >= run_length_sequence.size()-1) {
		throw std::invalid_argument("Ordinal search requires at "
			"least two ones in the MFM train search pattern.");
	}

	// Since ordinal search is based on differences between
	// adjacent flux delay values, it must necessarily start
	// matching at the second proper flux.
	++out.offset;

	for (size_t i = out.offset; i < run_length_sequence.size(); ++i) {
		// Skip any one-sided (at-least) constraints at the end as
		// we don't know their exact value.
		if (run_length_sequence[i].at_least) { continue; }
		out.needle.push_back(sign(run_length_sequence[i].run -
			run_length_sequence[i-1].run));
	}

	for (char digit: out.needle) {
		if (digit == 0) {
			throw std::logic_error("Can't make ordinal sequence with "
				"more than one of the same clock length in a row!");
		}
	}

	return out;
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

	auto match_pos = match_start;

	// First create the run length sequence to get the bands that
	// fluxes should be assigned to.
	std::vector<bound> run_length_sequence = get_run_length_coding(
		MFM_train_search_sequence.begin(),
		MFM_train_search_sequence.end());

	int num_bands = max_vec(run_length_sequence) + 1;

	std::vector<std::vector<int> > bands(num_bands);
	std::vector<std::vector<int> > at_least_bands(num_bands);

	// Do a very ugly impression of the zip operator.
	size_t i;

	for (i = 0; i < run_length_sequence.size() &&
		match_pos != fluxes_end; ++i) {

		if (run_length_sequence[i].at_least) {
			// Assign the flux at match_pos to the lowest band it
			// can be associated with if the match is correct.
			at_least_bands[run_length_sequence[i].run].push_back(
				*match_pos++);
		} else {
			// Assign the flux at match_pos to the band it would be
			// associated with if the match is correct.
			bands[run_length_sequence[i].run].push_back(
				*match_pos++);
		}
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
		clock_upper = std::min(clock_upper,
			4 * min_in_band[band] / (2.0 * band + 1));
		// Similarly, the lowest clock we can have is so that
		// round(max_in_band / half_clock) = band + 1,
		// max_in_band / half_clock = band + 1.5 - epsilon,
		// clock = 4 * max_in_band / (2 * band + 3)
		clock_lower = std::max(clock_lower,
			4 * max_in_band[band] / (2.0 * band + 3));
	}

	// For every one-sided band value, the clock rate must be sufficiently
	// low that the flux delay associated with this band value is at least
	// at that band (it may be higher).

	for (size_t band = 1; band < at_least_bands.size(); ++band) {
		for (auto pulse_delay: at_least_bands[band]) {
			clock_upper = std::min(clock_upper,
				4 * pulse_delay / (2.0 * band + 1));
		}
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

		// This may also be triggered if it's impossible
		// to pass the "at least this band" constraints; in that
		// case we're dealing with an actual false positive. I may
		// need to distinguish these later.
		return -2;
	}
}

// TODO? Somehow generalize so that it doesn't just take preambles.
// E.g. when we start looking for damaged IDAMs, we'll want to get
// results that are given by the IBM_preamble structure.

// Also TODO: using the offset to offset the preamble so that it matches
// the offset we gave the ordinal preamble... is kind of ugly.

std::vector<match_with_clock> get_flux_matches(
	const std::vector<int> & flux_transitions,
	const std::vector<search_result> & possible_matches,
	const IBM_preamble & preamble_info) {

	std::vector<match_with_clock> out;
	match_with_clock true_match;

	for (const search_result & result: possible_matches) {
		// Due to various effects, an ordinal match at x corresponds
		// to a flux match at x - offset. Get this offset.
		// See get_ordinal_search_sequence for more info.
		size_t offset = preamble_info.get_ordinal_offset_by_ID(
			result.ID);

		// If the match was right at the start, then the drive started
		// reading right at the edge of either a preamble or something
		// that looks a lot like it. Because we don't have the previous
		// byte, we can't tell if the match is real, so skip.
		if (offset > result.idx) {
			continue;
		}

		// Ugly HACK to deal with false positives.
		// I do not recommend doing this... we instead need to
		// explicitly introduce every preamble and then make
		// the ordinal search only search up till the first repeated
		// bit (where the strategy will no longer work). TODO as part
		// of the design cleanup.
		// Ultimately we also need to support overlapping sectors to be
		// robust to false positives.
		std::vector<char> amended_preamble =
			preamble_info.get_preamble_by_ID(result.ID);
		// The upper nibble following the preamble will always be high,
		// but ordinal search can't handle that. Manually push four
		// high MFM bits to deal with it.
		for (int i = 0; i < 4; ++i) {
			amended_preamble.push_back(0);
			amended_preamble.push_back(1);
		}

		double clock =	get_clock(amended_preamble,
			flux_transitions.begin() + result.idx - offset,
			flux_transitions.end());

		if (clock < 0) {
			continue; // false positive or impossible to fit clock
		}

		true_match.match_location = result.idx - offset;
		true_match.offset = offset;
		true_match.estimated_clock = clock;

		out.push_back(true_match);
	}

	return out;
}