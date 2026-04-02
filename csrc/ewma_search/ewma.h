#pragma once

// This file contains functions implementing a floppy EWMA/PLL-type decoder,
// in particular determining if starting at a haystack of flux transitions with
// a given initial clock and smoothing parameter alpha can produce a sequence
// of half-clock outputs that matches a given needle.

// To recap, the EWMA decoder works like this:

//		- Starting with an initial clock, repeatedly:
//			- Round the current flux transition delay times 2/clock
//			  to get the number of half-clocks for that flux transition
//			- Let the ideal clock be the clock that would put the current
//			  half-clock output in dead center of the band for that
//			  half-clock (i.e. minimially surprising location).
//			- Update the clock estimate by (1-alpha) * clock + alpha * ideal clock.

// The tests check the output half-clocks, comparing them to the desired
// output.

#include <vector>

#include "interval.h"

// Return a boolean indicating if, given initial clock and alpha,
// and haystack starting at haystack_pos, EWMA decoding returns the
// needle.

bool is_valid_ewma(double initial_clock,
	const std::vector<int> & haystack,
	const std::vector<int> & half_clock_needle,
	size_t haystack_pos, double alpha);

std::pair<bool, bool> is_valid_ewma(
	std::pair<int, interval> search_result,
	const std::vector<int> & haystack,
	const std::vector<int> & half_clock_needle,
	double alpha);

// Return a more granular output: either that we match the needle,
// or that the given initial clock is too low or too high, or that
// the current clock proves that no clock can produce the needle.

// The idea is that we take the given input clock and run it through the
// whole EWMA procedure while keeping track of whether any given iterate
// is below the lower bound for the desired number of half-clocks,
// or above the upper bound.

// The only difference from plain EWMA decoding is that we always assume
// that prior iterations hit the right number of half-clocks; this works
// because it's true when we're at the right clock, and it reduces the
// degree to which earlier results can upset later ones, thus making
// the function more well-behaved.

typedef enum {TOO_LOW, VALID, TOO_HIGH, IMPOSSIBLE} margin_dir;

margin_dir ewma_margin_direction(double initial_clock,
	const std::vector<int> & haystack,
	const std::vector<int> & half_clock_needle,
	size_t haystack_pos, double alpha);

