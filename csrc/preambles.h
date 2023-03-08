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

// The whole "offset" thing is also rather ugly. Either fix or explain better
// why it's there. Best would of course be to make ordinal search handle it
// transparently, offsetting back results.

#pragma once
#include <vector>
#include "ordinal_pattern.h"

// Search result IDs.
const int PREAMBLE_ID_A1 = 0, PREAMBLE_ID_C2 = 1;

class IBM_preamble {
	public:
		std::vector<char> short_A1 = {
			0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1};
		std::vector<char> short_C2 = {
			0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0};
		std::vector<char> A1_sequence, C2_sequence;
		std::vector<char> offset_A1_sequence, offset_C2_sequence;
		ordinal_pattern ordinal_A1_sequence, ordinal_C2_sequence;

		IBM_preamble();

		std::vector<char> get_preamble_by_ID(int ID) const;
		std::vector<char> get_offset_preamble_by_ID(int ID) const;
};