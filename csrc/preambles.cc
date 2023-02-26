// To generate ordinal search sequences.
#include "preambles.h"
#include "ordinal_search.h"
#include <stdexcept>

IBM_preamble::IBM_preamble() {
	// The A1 sequence is the short A1 repeated three times,
	// and ditto for the C2 sequence.

	A1_sequence.reserve(3 * short_A1.size());
	C2_sequence.reserve(3 * short_C2.size());
	for (int i = 0; i < 3; ++i) {
		std::copy(short_A1.begin(), short_A1.end(),
			std::back_inserter(A1_sequence));
		std::copy(short_C2.begin(), short_C2.end(),
			std::back_inserter(C2_sequence));
	}

	// The first bit of both sequences is a zero, and we
	// can't create an ordinal sequence from something starting
	// in a zero. So start from index one instead.
	offset_A1_sequence = std::vector<char>(
		A1_sequence.begin() + 1, A1_sequence.end());
	offset_C2_sequence = std::vector<char>(
		C2_sequence.begin() + 1, C2_sequence.end());
	ordinal_A1_sequence = get_ordinal_search_sequence(offset_A1_sequence);
	ordinal_C2_sequence = get_ordinal_search_sequence(offset_C2_sequence);
}

std::vector<char> IBM_preamble::get_preamble_by_ID(int ID) const {
	switch(ID) {
		case PREAMBLE_ID_A1: return A1_sequence;
		case PREAMBLE_ID_C2: return C2_sequence;
		default:
			throw std::logic_error("IBM_preamble: unknown preamble ID");
	}
}

std::vector<char> IBM_preamble::get_offset_preamble_by_ID(int ID) const {
	switch(ID) {
		case PREAMBLE_ID_A1: return offset_A1_sequence;
		case PREAMBLE_ID_C2: return offset_C2_sequence;
		default:
			throw std::logic_error("IBM_preamble: unknown preamble ID");
	}
}