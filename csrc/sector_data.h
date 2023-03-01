#pragma once

#include <vector>
#include "pulse_train.h"

// Level two.

// Higher order formatting for IBM floppies, from
// https://www-user.tu-chemnitz.de/~heha/basteln/PC/usbfloppy/floppy.chm/

class sector_data {
	public:
		std::vector<unsigned char> decoded_data;
		std::vector<unsigned char> errors;
		std::vector<size_t> MFM_train_indices;
};

// Since we want to reset the MFM state at every preamble so that
// we're guaranteed that errors don't propagate indefinitely, we
// decode one chunk at a time (from one preamble to the next).
// Hence start_index and end_index, which give the boundaries of
// these chunks.
sector_data decode_MFM_train(const MFM_train_data & MFM_train,
	size_t start_index, size_t end_index);
