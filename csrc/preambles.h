// Higher order formatting for IBM floppies, from
// https://www-user.tu-chemnitz.de/~heha/basteln/PC/usbfloppy/floppy.chm/

// The start of a chunk is marked by a synchronization header. There are
// four chunk types:

// Abbrev. Full name                 Byte seq.  What it is
// IAM	   Index Address Mark        C2C2C2 FC, start of a track
// IDAM	   ID Address Mark           A1A1A1 FE, start of a sector
// DAM	   Data Address Mark         A1A1A1 FB, start of sector data
//		   Deleted Data Address Mark A1A1A1 F8,	start of deleted sector data

// To find the chunks, we need to find the header. The following class
// defines them.

// Because we want to be robust, we're going to only search for the short
// MFM sequences; we don't make any assumptions about the gap data that
// precedes them, even though real world floppies usually have padding that
// helps the phase-locked loop lock on.

// The sequence vectors are given as MFM train bit sequences, i.e.
// *before* MFM decoding. Hence we don't have to deal with detecting
// out-of-band errors.

// Maybe I don't even need the full 3x. (Consider that later.)

#pragma once
#include <vector>

// To generate ordinal search sequences.
#include "ordinal_search.h"

class IBM_preamble {
	public:
		std::vector<char> short_A1 = {
			0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1};
		std::vector<char> short_C2 = {
			0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0};
		std::vector<char> A1_sequence, C2_sequence;
		std::vector<char> ordinal_A1_sequence, ordinal_C2_sequence;

		IBM_preamble();
};

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
	ordinal_A1_sequence = get_ordinal_search_sequence(std::vector<char>(
		A1_sequence.begin() + 1, A1_sequence.end()));
	ordinal_C2_sequence = get_ordinal_search_sequence(std::vector<char>(
		C2_sequence.begin() + 1, C2_sequence.end()));
}