#pragma once

// Given a needle (integer half-clocks) and a smoothing parameter value alpha,
// determine if there exists a position in the haystack (flux transitions) and
// an initial clock interval so that starting an EWMA decoder at that position
// with any clock inside that interval, with the given alpha, will return the
// needle sequence of half-clocks (offset MFM train values).

// This could then be used to find preambles (A1A1A1, C2C2C2) without needing
// to know the clock prior to the preamble; which in turn should be more
// robust to corruption, outliers, etc. that could throw off an ordinary
// floppy PLL.

#include <vector>

#include "interval.h"
#include "ewma.h"

const interval IMPOSSIBLE_INTERVAL = interval(1, -1);

// Can't be found due to numerical imprecision
const interval INDETERMINATE_INTERVAL = interval(2, -2);

const std::pair<int, interval> NO_EWMA_MATCH = {-1, IMPOSSIBLE_INTERVAL};

// Return a pair of (index into haystack, interval) giving the first match
// for the needle half-clock sequence in haystack. If there is none, return
// -1 and some interval that should be disregarded (currently
// IMPOSSIBLE_INTERVAL).

std::pair<int, interval> ewma_search(
	const std::vector<int> & haystack,
	const std::vector<int> & half_clock_needle,
	size_t haystack_pos,
	double alpha, double tolerance);

std::pair<int, interval> ewma_search(
	const std::vector<int> & haystack,
	const std::vector<int> & half_clock_needle,
	double alpha, double tolerance);

std::vector<std::pair<int, interval> > ewma_search_all(
	const std::vector<int> & haystack,
	const std::vector<int> & half_clock_needle,
	double alpha, double tolerance);
