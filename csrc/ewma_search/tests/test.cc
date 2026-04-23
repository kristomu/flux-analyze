
#include <vector>
#include <iostream>
#include <iterator>
#include <stdexcept>

#include <math.h>
#include <stdlib.h>

#include <gtest/gtest.h>

#include "ewma_search/interval.h"
#include "ewma_search/ewma.h"
#include "ewma_search/ewma_search.h"

// For testing - to determine if the sophisticated search
// algorithm matches iff there is a match. (Up to some margin
// of confidence; test can never decide the matter with certainty,
// of course.)
std::pair<int, interval> brute_force_ewma_search(
	const std::vector<int> & haystack,
	const std::vector<int> & half_clock_needle,
	double alpha) {

	if (half_clock_needle.size() > haystack.size() || half_clock_needle.size() == 0) {
		return NO_EWMA_MATCH;
	}

	for (size_t i = 0; i < half_clock_needle.size(); ++i) {
		if (half_clock_needle[i] == 0) {
			throw std::invalid_argument("Needle half-clock values"
			" of zero are not supported.");
		}
		if (half_clock_needle[i] < 0) {
			throw std::invalid_argument("Negative needle values"
			" are not supported.");
		}
	}

	for (size_t i = 0; i < haystack.size() - half_clock_needle.size(); ++i) {

		bool found_something = false;
		interval matching_clock;
		matching_clock.lower_clock = 1000; // some very large positive value
		matching_clock.upper_clock = -1000; // ditto but negative

		for (int initial_clock_idx = 1; initial_clock_idx < 10000; ++initial_clock_idx) {

			double initial_clock = initial_clock_idx/100.0;
			double clock = initial_clock;
			bool possible_match = true;

			for (size_t j = 0; j < half_clock_needle.size() && possible_match; ++j) {

				double half_clocks = round(haystack[j+i] * 2.0 / clock);

				if (half_clocks != half_clock_needle[j]) {
					possible_match = false;
					continue;
				}

				double ideal_clock = (2.0 * haystack[i+j]) / half_clocks;

				clock = (1-alpha) * clock + alpha * ideal_clock;
			}

			if (possible_match) {
				found_something = true;
				matching_clock.grow_bounds(initial_clock, initial_clock);
			}
		}
		if (found_something) {
			return {i, matching_clock};
		}
	}

	return NO_EWMA_MATCH;
}

std::pair<int, interval> half_brute_force_ewma_search(
	const std::vector<int> & haystack,
	const std::vector<int> & half_clock_needle,
	double alpha) {

	if (half_clock_needle.size() > haystack.size() || half_clock_needle.size() == 0) {
		return NO_EWMA_MATCH;
	}

	for (size_t i = 0; i < half_clock_needle.size(); ++i) {
		if (half_clock_needle[i] == 0) {
			throw std::invalid_argument("Needle half-clock values"
			" of zero are not supported.");
		}
		if (half_clock_needle[i] < 0) {
			throw std::invalid_argument("Negative needle values"
			" are not supported.");
		}
	}

	for (size_t i = 0; i < haystack.size() - half_clock_needle.size(); ++i) {

		margin_dir test = IMPOSSIBLE;

		interval matching_clock;
		matching_clock.lower_clock = 1e-8;
		matching_clock.upper_clock = 1e-8;

		bool found_something = false;

		while (ewma_margin_direction(matching_clock.lower_clock,
			haystack, half_clock_needle, i, alpha) == TOO_LOW) {
			matching_clock.lower_clock *= 2;
		}

		matching_clock.lower_clock /= 2.0;
		matching_clock.upper_clock = matching_clock.lower_clock;

		do {
			matching_clock.upper_clock *= 2;

			test = ewma_margin_direction(matching_clock.upper_clock,
			haystack, half_clock_needle, i, alpha);
		} while (test != TOO_HIGH && test != IMPOSSIBLE);

		if (test == IMPOSSIBLE) {
			continue;
		}

		double b = matching_clock.lower_clock;
		double a = (matching_clock.upper_clock - matching_clock.lower_clock);

		for (int initial_clock_idx = 1; initial_clock_idx < 10000; ++initial_clock_idx) {

			double initial_clock = a * initial_clock_idx/10000.0 + b;
			double clock = initial_clock;
			bool possible_match = true;

			for (size_t j = 0; j < half_clock_needle.size() && possible_match; ++j) {

				double half_clocks = round(haystack[j+i] * 2.0 / clock);

				if (half_clocks != half_clock_needle[j]) {
					possible_match = false;
					continue;
				}

				double ideal_clock = (2.0 * haystack[i+j]) / half_clocks;

				clock = (1-alpha) * clock + alpha * ideal_clock;
			}

			if (possible_match) {
				found_something = true;
				matching_clock.grow_bounds(initial_clock, initial_clock);
			}
		}
		if (found_something) {
			return {i, matching_clock};
		}
	}

	return NO_EWMA_MATCH;
}

void run_test(const std::vector<int> & haystack,
	const std::vector<int> & needle, double alpha) {

	std::cout << "\n\nNow testing: Haystack:\t";
	std::copy(haystack.begin(), haystack.end(), std::ostream_iterator<int>(std::cout, " "));
	std::cout << "\nNow testing: Needle:\t";
	std::copy(needle.begin(), needle.end(), std::ostream_iterator<int>(std::cout, " "));
	std::cout << "\n";
	std::cout << "alpha = " << alpha << "\n";

	//haystack = needle;

	double tolerance = 1e-10;

	std::pair<int, interval> search_result = /*ewma_search(
		haystack, needle, alpha, tolerance);*/half_brute_force_ewma_search(
		haystack, needle, alpha);

	std::pair<int, interval> brute_result = brute_force_ewma_search(
		haystack, needle, alpha);

	std::cout << "Search: " << search_result.first << " ("
		<< search_result.second.lower_clock << " - "
		<< search_result.second.upper_clock
		<< ")\t Brute: " << brute_result.first << " ("
		<< brute_result.second.lower_clock << " - "
		<< brute_result.second.upper_clock << ")\n";

	// No validity to check here, move along.
	if (search_result.first == -1 && brute_result.first == -1) {
		return;
	}

	std::pair<bool, bool> search_validity, brute_validity;

	if (search_result.first == -1) {
		search_validity = {false, false}; // Well, N/A
	} else {
		search_validity = is_valid_ewma(search_result,
			haystack, needle, alpha);
	}

	if (brute_result.first == -1) {
		brute_validity = {false, false};
	} else {
		brute_validity = is_valid_ewma(brute_result,
			haystack, needle, alpha);
	}

	if (search_result.first != -1 && brute_result.first == -1) {
		std::cout << "Search found something brute couldn't.\n";
		return;
	}

	std::cout << "Search interval validity ";
	if (search_validity.first) { std::cout << "OK\t"; 
	} else { std::cout << "Nope\t";}
	if (search_validity.second) { std::cout << "OK\n"; 
	} else { std::cout << "Nope\n";}

	std::cout << "Brute interval validity ";
	if (brute_validity.first) { std::cout << "OK\t"; 
	} else { std::cout << "Nope\t";}
	if (brute_validity.second) { std::cout << "OK\n"; 
	} else { std::cout << "Nope\n";}

	double brute_lower_margin = search_result.second.lower_clock - brute_result.second.lower_clock,
		brute_upper_margin = brute_result.second.upper_clock - search_result.second.upper_clock;

	if (brute_lower_margin > tolerance || brute_upper_margin > tolerance) {
		std::cout << "Warning: Brute-force returns larger interval than EWMA search!\n";
		std::cout << "Brute is " << std::max(brute_lower_margin, brute_upper_margin) << " wider.\n";
	}

	if ((brute_result.first != -1 && search_result.first == -1) || 
		(brute_result.first >= 0 && brute_result.first < search_result.first)) {
		throw std::runtime_error("Brute found result that search didn't.\n");
	}
}

int test_fuzz_input(const char *Data, long long Size) {
	// The data we need is
	//		alpha				8 bytes			double
	//		haystack length		2 bytes			u_short
	//		needle length		2 bytes			u_short
	//		haystack			k * 2 bytes		vector<int>
	//		needle				j * 1 byte		vector<char>

	// Too small?
	if (Size < 10 + 2) {
		return 0;
	}

	double alpha = *(((double *)Data));
	unsigned short haystack_len = *((unsigned short *)(Data + 8));
	unsigned short needle_len = *((unsigned short *)(Data + 10));

	// If it exceeds the size, then clamp to most of the string,
	// but leave some for the needle.
	while (haystack_len * 4 + needle_len >= Size - 12) {
		if (haystack_len == 0 || needle_len == 0) {
			return 0;
		}
		--haystack_len;
		--needle_len;
	}

	if (!finite(alpha)) {
		return 0; // out of spec alpha
				// Hack to fix T2, but the problem itself remains.
	}

	// out of spec alpha, scale it down
	if (fabs(alpha) > 1) {
		while (fabs(alpha) > 1) { alpha /= 2.0;}
	}

	std::vector<int> haystack;

	for (int byte_off = 0; byte_off < haystack_len * 4; byte_off += 4) {
		haystack.push_back(*((int *)(Data + 12 + byte_off)));
	}

	// Everything else goes to the needle.

	std::vector<int> needle;
	double tolerance = 1e-10;

	for (int byte = 12 + haystack_len * 4; byte < 12 + haystack_len * 4 + needle_len; ++byte) {
		char next_byte = *((char *)(Data + byte));
		//if (next_byte == 0) { next_byte = 1; } // not legal, see <T1>
		needle.push_back(next_byte);
	}

	/*std::cout << "\n\nNow testing: Haystack:\t";
	std::copy(haystack.begin(), haystack.end(), std::ostream_iterator<int>(std::cout, " "));
	std::cout << "\nNow testing: Needle:\t";
	std::copy(needle.begin(), needle.end(), std::ostream_iterator<int>(std::cout, " "));
	std::cout << "\n";
	std::cout << "Alpha: " << alpha << "\n";*/

	//std::cout << "Sophisticated search" << std::endl;

	std::pair<int, interval> search_result, brute_result;

	try {
		search_result = ewma_search(
			haystack, needle, alpha, tolerance);
	} catch (const std::invalid_argument & e) {
		// anticipated error
		search_result = NO_EWMA_MATCH;
	}

	//std::cout << "Brute force" << std::endl;

	try {
		brute_result = half_brute_force_ewma_search(
			haystack, needle, alpha);
		if (brute_result.first == -1) {
			brute_result = brute_force_ewma_search(
			haystack, needle, alpha);
		}
	} catch (const std::invalid_argument & e) {
		// anticipated error
		brute_result = NO_EWMA_MATCH;
	}

	if (brute_result.first == -1) {
		return 0;
	}

	if (brute_result.first != -1 && search_result.first == -1) {
		std::cout << "Brute force found something (I)\n";
		std::cout << brute_result.first << " vs " << search_result.first << "\n";
		run_test(haystack, needle, alpha);
		abort(); // Error: Brute found something that search didn't.
	}

	if (brute_result.first < search_result.first) {
		std::cout << "Brute force found earlier match\n";
		std::cout << "Brute: " << brute_result.first << " Search: " << search_result.first << "\n";
		run_test(haystack, needle, alpha);
		abort(); // Error: Brute found an earlier match than search.
	}

	double brute_lower_margin = search_result.second.lower_clock - brute_result.second.lower_clock,
		brute_upper_margin = brute_result.second.upper_clock - search_result.second.upper_clock;

	if (brute_result.first == search_result.first && 
		(brute_lower_margin > tolerance || brute_upper_margin > tolerance)) {
		std::cout << "Tolerance exceeded (lower " << brute_lower_margin << ", upper " << brute_upper_margin << ")\n";
		run_test(haystack, needle, alpha);
		abort(); // One of the bounds on brute is much better than search's.
	}

	return 0;
}

/*int main() {
	run_selected_tests();
}*/
