#include <limits>
#include <math.h>

#include "ewma.h"

bool is_valid_ewma(
	double initial_clock,
	const std::vector<int> & haystack,
	const std::vector<int> & half_clock_needle,
	size_t haystack_pos,
	double alpha) {

	double clock = initial_clock;

	bool possible_match = true;

	for (size_t j = 0; j < half_clock_needle.size() && possible_match; ++j) {

		int half_clocks = round(haystack[haystack_pos+j] * 2.0 / clock);

		if (half_clocks != half_clock_needle[j]) {
			possible_match = false;
			continue;
		}

		double ideal_clock = (2.0 * haystack[haystack_pos+j]) / half_clocks;

		clock = (1-alpha) * clock + alpha * ideal_clock;
	}

	return possible_match;
}

// Determine whether the input clock fails because it is too low
// or too high. The idea is that we take the given input clock and
// run it through the whole (modified) EWMA procedure while keeping track
// of whether any given iterate is below the lower bound for the desired
// number of half-clocks, or above the upper bound.

std::pair<double, double> ewma_margin (
	double initial_clock,
	const std::vector<int> & haystack,
	const std::vector<int> & half_clock_needle,
	size_t haystack_pos,
	double alpha) {

	// Use greater precision internally than the provided parameters
	// so we're less likely to get numerical precision problems.
	long double clock = initial_clock;

	bool possible_match = true;

	long double lower_divergence = std::numeric_limits<long double>::infinity(),
		upper_divergence = std::numeric_limits<long double>::infinity();

	for (size_t j = 0; j < half_clock_needle.size() && possible_match; ++j) {

		// For round(haystack[haystack_pos+j] * 2.0 / clock) to equal
		// half_clock_needle[j], it must be strictly between the
		// lower and upper bounds given below; they correspond to the
		// term inside the round() expression being half_clock_needle[j] + 0.5
		// and half_clock_needle[j] - 0.5 respectively.

		long double lower_bound = (2.0 * haystack[haystack_pos+j])/(half_clock_needle[j] + 0.5);
		long double upper_bound = (2.0 * haystack[haystack_pos+j])/(half_clock_needle[j] - 0.5);

		lower_divergence = std::min(lower_divergence, clock-lower_bound);
		upper_divergence = std::min(upper_divergence, upper_bound - clock);
		
		long double ideal_clock = (2.0 * haystack[haystack_pos+j]) / half_clock_needle[j];

		clock = (1-alpha) * clock + alpha * ideal_clock;
	}

	return {lower_divergence, upper_divergence};
}

margin_dir ewma_margin_direction(double initial_clock,
	const std::vector<int> & haystack,
	const std::vector<int> & half_clock_needle,
	size_t haystack_pos,
	double alpha) {

	std::pair<double, double> margins = ewma_margin(
		initial_clock, haystack, half_clock_needle,
		haystack_pos, alpha);

	if (margins.first > 0 && margins.second > 0) {
		return VALID; // The given clock is a valid solution.
	}

	// If either margin is exactly zero, resolve by checking if
	// the given initial clock is actually valid. (This should fix
	// some numerical precision problems in extreme cases, but may
	// also be overkill.)

	if (margins.first >= 0 && margins.second >= 0
		&& is_valid_ewma(initial_clock, haystack, half_clock_needle,
			haystack_pos, alpha)) {

		return VALID;
	}

	if (margins.first <= 0 && margins.second <= 0 ) {
		return IMPOSSIBLE; // The clock proves that no solution exists.
	}

	if (margins.first <= 0) {
		return TOO_LOW;
	}

	return TOO_HIGH;
}

std::pair<bool, bool> is_valid_ewma(
	std::pair<int, interval> search_result,
	const std::vector<int> & haystack,
	const std::vector<int> & half_clock_needle,
	double alpha) {

	return {
		is_valid_ewma(search_result.second.lower_clock,
			haystack, half_clock_needle, search_result.first, alpha),

		is_valid_ewma(search_result.second.upper_clock,
			haystack, half_clock_needle, search_result.first, alpha)
	};
}
