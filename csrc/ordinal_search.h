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

#pragma once

#include "tools.h"
#include "rabin_karp.h"
#include "ordinal_pattern.h"
#include "preambles.h"
#include <vector>

struct match_with_clock {
	size_t match_location;
	double estimated_clock;
};

// Returns the sign values of successive differences, e.g. 
// 3 2 2 3 becomes -1 0 1.

template<typename T> std::vector<char> get_delta_coding(
	const T & in_vector) {
	std::vector<char> out;

	for (size_t i = 1; i < in_vector.size(); ++i) {
		out.push_back(sign(in_vector[i] - in_vector[i-1]));
	}

	return out;
}

ordinal_pattern get_ordinal_search_sequence(
	const std::vector<char> & MFM_train_search_sequence);

double get_clock(const std::vector<char> & MFM_train_search_sequence,
	std::vector<int>::const_iterator match_start,
	std::vector<int>::const_iterator fluxes_end);

// This filters out a bunch of possible false positives.

std::vector<match_with_clock> filter_matches(
	const std::vector<int> & flux_transitions,
	const std::vector<search_result> & possible_matches,
	const IBM_preamble & preamble_info);