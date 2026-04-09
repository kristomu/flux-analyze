#include <stdexcept>
#include <vector>

#include <math.h>

#include "ewma_search.h"

// Suppose that we have confirmed that there's an EWMA match
// starting at haystack[haystack_pos] with the given alpha. This
// function finds the maximal closed interval of valid choices
// for the initial clock, up to the given absolute tolerance,
// subject to numerical precision limits. If there's no actual
// match, it will return IMPOSSIBLE_INTERVAL.

interval get_boundary(
	const std::vector<int> & haystack,
	const std::vector<int> & half_clock_needle,
	size_t haystack_pos, double alpha, double tolerance) {

	// TODO: If needle is empty, return the full infinite interval.
	// (Once I've got the tests properly implemented.)

	// The interval of valid clocks must be inside
	// the first interval (almost by definition). So
	// start with it and then use binary search to find
	// the largest interval inside.

	interval current_interval;

	current_interval.shrink_bounds(
		(2.0 * haystack[haystack_pos])/(half_clock_needle[0] + 0.5),
		(2.0 * haystack[haystack_pos])/(half_clock_needle[0] - 0.5));

	// First search: find *some* point that's valid.

	margin_dir margin = IMPOSSIBLE;
	double low = current_interval.lower_clock,
		high = current_interval.upper_clock, mid,
		initial_valid_point;

	// This triggers if divisions stop working and we thus
	// can't get within the desired absolute tolerance. This
	// can happen e.g. with very large absolute values.
	bool numerical_error = false;

	while (margin != VALID) {

		mid = (low + high) / 2;

		margin = ewma_margin_direction(mid, haystack,
			half_clock_needle, haystack_pos, alpha);

		if (margin == TOO_LOW) {
			if (low == mid) {
				numerical_error = true;
			}
			low = mid;
		} else {
			if (high == mid) {
				numerical_error = true;
			}
			high = mid;
		}

		if (margin == IMPOSSIBLE || high-low < 1e-10) {
			return IMPOSSIBLE_INTERVAL; // or exception?
		}

		if (numerical_error) {
			return INDETERMINATE_INTERVAL;
		}
	}

	initial_valid_point = mid;

	// Find the lower bound.
	double invalid_point = current_interval.lower_clock,
		valid_point = initial_valid_point;

	numerical_error = false;

	while (!numerical_error && fabs(valid_point - invalid_point) > tolerance) {

		mid = (invalid_point + valid_point) / 2;

		if (is_valid_ewma(mid, haystack, half_clock_needle,
			haystack_pos, alpha)) {

			if (mid == valid_point) {
				numerical_error = true;
			}

			valid_point = mid;
		} else {
			if (mid == invalid_point) {
				numerical_error = true;
			}

			invalid_point = mid;
		}

	}

	current_interval.lower_clock = valid_point;

	// Find the upper bound.
	invalid_point = current_interval.upper_clock;
	valid_point = initial_valid_point;

	numerical_error = false;

	while (!numerical_error && fabs(valid_point - invalid_point) > tolerance) {

		mid = (invalid_point + valid_point) / 2;

		if (is_valid_ewma(mid, haystack, half_clock_needle,
			haystack_pos, alpha)) {

			if (mid == valid_point) {
				numerical_error = true;
			}

			valid_point = mid;
		} else {
			if (mid == invalid_point) {
				numerical_error = true;
			}

			invalid_point = mid;
		}

	}

	current_interval.upper_clock = valid_point;

	return current_interval;
}

std::pair<int, interval> ewma_search(
	const std::vector<int> & haystack,
	const std::vector<int> & half_clock_needle,
	size_t haystack_pos,
	double alpha, double tolerance) {

	// If the needle is too long, we obviously can't have a match.
	if (haystack_pos + half_clock_needle.size() > haystack.size()) {
		return {-1, IMPOSSIBLE_INTERVAL};
	}

	// If alpha values are too large or not even finite, then bail.
	if (!finite(alpha) || fabs(alpha) > 1) {
		throw std::invalid_argument("Alpha out of bounds.");
	}

	// Check for needle values of zero that would otherwise cause
	// division by zero problems. They should never appear in MFM
	// streams. (It is possible to get around this later if needed
	// for e.g. copy protection.)
	for (size_t i = 0; i < half_clock_needle.size(); ++i) {
		if (half_clock_needle[i] == 0) {
			throw std::invalid_argument("Needle half-clock values"
			" of zero are not supported.");
		}

		// Ditto negative values. This is technically possible to
		// deal with by recording the sign value and flipping
		// bounds, but I neve use negative needle values, so...
		if (half_clock_needle[i] < 0) {
			throw std::invalid_argument("Negative needle values"
			" are not supported.");
		}

	}

	for (size_t i = haystack_pos; i < haystack.size() - half_clock_needle.size(); ++i) {

		interval current_interval;

		bool possible_match = true;

		for (size_t j = 0; j < half_clock_needle.size() && possible_match; ++j) {

			// Let L_i be the lower bound of the clock at i,
			// and U_i be its upper bound. Let y_i be the number
			// of half-clocks at position i.

			// y_i = round(2 * flux(i)/clock(i)).

			// So needle(j) - 0.5 < 2 * flux(i)/clock(i) < needle(j) + 0.5

			// (2 flux(i))/(needle(j) + 0.5) < clock(i) < (2 flux(i)/(needle(j) - 0.5)

			// So L_i = std::max(L_i, (2 flux(i))/(needle(j) + 0.5))
			// U_i = std::min(U_i, (2 flux(i)/(needle(j) - 0.5))

			// These are the bounds also described in ewma.cc.

			// Strictly speaking, we would need to add an epsilon to shrink the
			// interval to the interior, i.e. to make the endpoints also valid.
			// This means that the interval we derive may not be entirely
			// correct. Due to numerical error, we can't just reverse the EWMA
			// calculations either; get_interval will be used to find the proper
			// interval on match instead.

			current_interval.shrink_bounds(
				(2.0 * haystack[i+j])/(half_clock_needle[j] + 0.5),
				(2.0 * haystack[i+j])/(half_clock_needle[j] - 0.5));

			possible_match &= current_interval.valid();

			// Now, consistent with the description in ewma.cc of how the
			// EWMA decoder works, we need to update the whole interval
			// based on the influence of the current "ideal clock" 

			// clock^(i) = (2 * flux(i)) / needle(j),

			// which would make the flux dead-on target.

			// This can be applied directly to each interval endpoint because
			// the transformation is montone in clock(i).

			double ideal_clock = (2.0 * haystack[i+j]) / half_clock_needle[j];

			current_interval.lower_clock = (1-alpha) *
				current_interval.lower_clock + alpha *
				ideal_clock;

			current_interval.upper_clock = (1-alpha) *
				current_interval.upper_clock + alpha *
				ideal_clock;

			possible_match &= current_interval.valid();
		}

		if (!current_interval.is_set || !possible_match) {
			continue;
		}

		// We found a match. Use binary search to find its interval,
		// and be sure to skip false positives.

		current_interval = get_boundary(haystack,
			half_clock_needle, i, alpha, tolerance);

		if (current_interval.valid()) {
			return {i, current_interval};
		}
	}

	return {-1, IMPOSSIBLE_INTERVAL};
}

std::pair<int, interval> ewma_search(
	const std::vector<int> & haystack,
	const std::vector<int> & half_clock_needle,
	double alpha, double tolerance) {

	return ewma_search(haystack, half_clock_needle, 0,
		alpha, tolerance);
}

std::vector<std::pair<int, interval> > ewma_search_all(
	const std::vector<int> & haystack,
	const std::vector<int> & half_clock_needle,
	double alpha, double tolerance) {

	size_t next = 0;

	std::vector<std::pair<int, interval> > matches;

	while (next < haystack.size()) {
		std::pair<int, interval> next_match =
			ewma_search(haystack, half_clock_needle,
				next, alpha, tolerance);

		if (next_match.first != -1) {
			next = next_match.first + 1;
			matches.push_back(next_match);
		} else {
			// We've searched through everything, so
			// next position is at the end of the
			// haystack.
			next = haystack.size();
		}
	}

	return matches;
}
