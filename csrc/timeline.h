#pragma once
#include "pulse_train.h"
#include "sector_data.h"
#include <stdexcept>
#include <list>
#include <set>

#include <stdlib.h>

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

enum slice_status {TS_UNKNOWN,
	TS_PREAMBLE_FOUND,
	TS_DECODED_UNKNOWN, // unknown preambles, say
	TS_DECODED_BAD,
	TS_DECODED_OK,
	TS_TRUNCATED};

// When splitting a timeslice, we have to determine which of the parts
// gets to retain its status.
// The strategies are PRESERVE_FIRST, PRESERVE_LAST, and PRESERVE_NEITHER.

enum split_strategy {PRESERVE_FIRST,
	PRESERVE_LAST, PRESERVE_NEITHER};

class timeslice {
	public:
		// TODO, disambiguate from flux_record or refactor it.
		std::vector<int> flux_data;
		MFM_train_data mfm_train;
		sector_data sec_data;

		// This will be used to connect auxiliary data to a timeslice
		// and enable iteration through still-unknown/non-decoded
		// timeslices.
		uint64_t ID = 0;

		// These offsets give the relative position of the beginning
		// of each chunk, relative to the first chunk in the timeline.
		size_t flux_data_begin = 0;
		size_t mfm_train_begin = 0;
		size_t sector_data_begin = 0;

		// TODO: Use this to connect to auxiliary data (usually address
		// marks) so we don't pollute this structure too much.
		uint64_t uuid = 0;

		// Various information relating to the contents of this slice
		// goes here. A slice should ideally only contain at most one
		// thing of interest.

		slice_status status = TS_UNKNOWN;

		// We need to store the clock estimate, if any. (This may be
		// replaced by explicit bands later.) -1 indicates no clock
		// estimate is known. (A neat solution would be to make this
		// mutable and then call, lazy style, a clock estimator if
		// the clock is unknown.)
		double clock_value = -1; // Unknown clock.

		// For now, if something has a preamble, then the sector data
		// starts at that point to avoid problems that would otherwise
		// ensue from a lack of byte alignment. TODO? It's kind of a hack.
		// Could I fix this? Perhaps with an iterator into a bit field.
		// ... it's an awful lot of complexity just to make a point, though.

		size_t preamble_offset; // in relative mfm_train offset terms.

		size_t flux_data_end() const {
			// Keeping the offsets arranged properly is tricky, so
			// let's add a bunch of tests to make sure it breaks early.
			size_t candidate = flux_data_begin + flux_data.size();

			if (!mfm_train.flux_indices.empty() &&
				candidate < *mfm_train.flux_indices.rbegin() +
				flux_data_begin) {
				throw std::logic_error("Timeslice ordering inconsistency:"
					" flux data end doesn't match MFM train");
			}
			return candidate;
		}

		size_t mfm_train_end() const {
			size_t candidate = mfm_train_begin + mfm_train.data.size();

			if (!sec_data.MFM_train_indices.empty() &&
				candidate < *sec_data.MFM_train_indices.rbegin() +
				mfm_train_begin) {
				throw std::logic_error("Timeslice ordering inconsistency:"
					" MFM train end doesn't match sector data");
			}

			return candidate;
		}

		size_t sector_data_end() const {
			return sector_data_begin + sec_data.decoded_data.size();
		}
};

class timeline {
	private:
		// Inserts a timeslice to the end of the line. It also handles
		// ll the accounting with offsets. NOTE: The previous timeslice
		// must have decoded flux data, TODO: fix that.
		void insert(timeslice & next, int ID);

		// This finds a random ID that's not in use.
		// TODO: Properly randomize, this is a sketch.
		int get_unused_ID() const;

	public:
		std::list<timeslice> timeslices;
		std::set<uint64_t> used_IDs;

		// This function sets the ID of the inserted timeslice to a
		// random value. ()
		void insert(timeslice & next);

		// This splits the given timeslice at the given (zero-indexed) sector
		// data point.
		void split(std::list<timeslice>::iterator & to_split,
			size_t first_byte_after, split_strategy strategy);

		// Update indices in case some data structures were changed.
		// I kind of feel like this is ugly because we're not supposed to
		// reach into the internals of the timeslices, because that can lead
		// to the views going out of sync. On the other hand, if we don't
		// reach into the internals, then we'll need tons of getters and
		// setters. Eugh.
		// The positions are really just memoized data.
		void update(std::list<timeslice>::iterator & to_update);
};