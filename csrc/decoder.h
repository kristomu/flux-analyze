#pragma once

#include <map>
#include <vector>

#include "address_marks.h"
#include "pulse_train/pulse_train.h"
#include "sector_data.h"
#include "preambles.h"

#include "timeline.h"

class decoded_tracks {
	public:
		int last_track = 0;
		int last_decoded_sector = 0;

		std::map<IDAM, DAM> sector_data;

		// Add sectors from another track to this one.
		decoded_tracks& operator+= (const decoded_tracks & rhs) {
			last_track = std::max(last_track, rhs.last_track);
			last_decoded_sector = std::max(last_decoded_sector,
				rhs.last_decoded_sector);

			sector_data.insert(rhs.sector_data.begin(),
				rhs.sector_data.end());

			return *this;
		}
};

// For keeping track of decoding performance.

class decoder_stats {
	public:
		size_t failures;
		size_t recovered_sectors;
		size_t unique_metadata_chunks;
		size_t total_timeslices;
};

class decoder {
	private:
		// Returns whether a given DAM has a sufficiently
		// close IDAM that it can be matched with.
		bool has_close_IDAM(bool has_last_IDAM,
			const address_mark & candidate_DAM,
			const address_mark & candidate_IDAM) const;

		address_mark deserialize(
			std::vector<unsigned char> & raw_bytes,
			size_t byte_stream_start, bool has_last_IDAM,
			const address_mark last_IDAM,
			bool verbose) const;

	public:
		// Adds OK sectors to the given decoded_tracks
		// structure.
		decoder_stats decode(timeline & line_to_decode,
			decoded_tracks & decoded, bool verbose);

		// For debugging.
		void dump_all_to_file(
			const timeline & line_to_dump,
			std::string data_filename,
			std::string error_filename) const;

		// Dumps the timeslices' sector data to files
		// starting in prefix. NOTE: This might be
		// complete garbage for unknowns because we generally
		// don't know what the clock is if the timeslice is unknown.
		void dump_sector_files(const timeline & line_to_dump,
			std::string prefix) const;

		// Dumps the image according to the decoded tracks to the
		// files "file_prefix.img" for the image and "file_prefix.mask"
		// for the bitmask file that shows which sectors were
		// actually decoded (0xFF for a byte that, at the same
		// position in prefix.img, was decoded; 0x00 for a byte
		// that wasn't.)
		// Maybe this should be a global function instead? TODO,
		// find out.
		void dump_image(const decoded_tracks & d_tracks,
			std::string file_prefix, int tracks,
			int heads, int sectors_per_track,
			int default_sector_size) const;
};