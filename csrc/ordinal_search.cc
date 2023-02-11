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
#include <stdexcept>

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
	// form our search sequence.

	char run_length = 0;
	std::vector<char> run_length_sequence;

	for (size_t i = 1; i < MFM_train_search_sequence.size(); ++i) {
		if (MFM_train_search_sequence[i] == 1) {
			run_length_sequence.push_back(run_length);
			run_length = 0;
		} else {
			++run_length;
		}
	}

	std::vector<char> ordinal_search_sequence = delta_code(
		run_length_sequence);

	for (char digit: ordinal_search_sequence) {
		if (digit == 0) {
			throw std::logic_error("Can't make ordinal sequence with "
				"more than one of the same clock length in a row!");
		}
	}

	return ordinal_search_sequence;
}

// TODO here: clock determination. See Python code for how to do this.