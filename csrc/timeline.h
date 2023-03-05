#pragma once
#include "pulse_train.h"
#include "sector_data.h"
#include <stdexcept>
#include <list>

// Boundary conditions might be a problem... how do I get around that
// short of storing the various decoding functions' internal state at
// the very step where each timeslice ends? And then necessarily at
// *every* step to handle split and merge operations? Eh... however
// we handle it, it's gonna be ugly.

// A concrete example of this problem: suppose chunk 1 is corrupted
// so that chunk 2 starts its A1A1A1 at what would be the 7th bit.
// Then to be robust, we must start chunk 2 at a new byte boundary,
// but that means that the sector data's final bit must be removed.

// Perhaps it would be better to only have the flux data chunks and
// then to decode on demand for each chunk? But translating back
// (cutting off what hasn't been decoded yet) would be a serious pain.

// And look at that, I jump right into such a problem, because the
// preambles aren't flux-aligned... the QND solution would probably
// be to have a preamble offset for this inside the timeslice, which
// we'll need anyway to distinguish timeslices with preambles from
// those without.

class timeslice {
	public:
		// TODO, disambiguate from flux_record or refactor it.
		std::vector<int> flux_data;
		MFM_train_data mfm_train;
		sector_data sec_data;

		// These offsets give the relative position of the beginning
		// of each chunk, relative to the first chunk in the timeline.
		size_t flux_data_begin = 0;
		size_t mfm_train_begin = 0;
		size_t sector_data_begin = 0;

		// Various information relating to the contents of this slice
		// goes here. A slice should ideally only contain at most one
		// thing of interest.

		// For now, if something has a preamble, then the sector data
		// starts at that point. TODO: Fix this. Probably involving an
		// iterator into the bit field of sector_data or somesuch.

		bool has_preamble = false;
		size_t preamble_offset; // in relative mfm_train offset terms.

		size_t flux_data_end() const {
			// Keeping the offsets arranged properly is tricky, so
			// let's add a bunch of tests to make sure it breaks early.
			size_t candidate = flux_data_begin + flux_data.size();

			if (candidate < *mfm_train.flux_indices.rbegin() +
				flux_data_begin) {
				throw std::logic_error("Timeslice ordering inconsistency:"
					" flux data end");
			}
			return candidate;
		}

		size_t mfm_train_end() const {
			size_t candidate = mfm_train_begin + mfm_train.data.size();

			if (candidate < *sec_data.MFM_train_indices.rbegin() +
				mfm_train_begin) {
				throw std::logic_error("Timeslice ordering inconsistency:"
					" MFM train end");
			}

			return candidate;
		}

		size_t sector_data_end() const {
			return sector_data_begin + sec_data.decoded_data.size();
		}
};

class timeline {
	public:
		std::list<timeslice> timeslices;

		// Do all the accounting with offsets. NOTE: The previous
		// timeslice must have been fully decoded (at all three
		// levels).
		void insert(timeslice & next);

		// TODO: Some way to split a timeslice (for the refinement phase when
		// excluding everything that belongs to a certain address mark chunk).
		// It needs to update the linear (offset) view.
};