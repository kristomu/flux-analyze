#pragma once

#include <map>
#include <vector>

#include "address_marks.h"
#include "pulse_train/pulse_train.h"
#include "sector_data.h"
#include "preambles.h"

#include "timeline.h"

#include <iostream>

// For keeping track of decoding performance.

std::vector<bool> boolean_or(std::vector<bool> a,
	std::vector<bool> b);

size_t count(const std::vector<bool> & a);

class decoder_stats {
	public:
		std::vector<bool> is_sector_recovered;

		size_t failures;
		size_t num_recovered_sectors;
		size_t num_sectors;
		// These will have to wait until I find a proper
		// union operation for them. TODO.

		//size_t unique_metadata_chunks;
		//size_t total_timeslices;

		// The parameter below is a giant HACK. Fix later.
		double flux_fraction_good = -1;

		decoder_stats& operator+= (const decoder_stats & rhs) {
			failures += rhs.failures;

			is_sector_recovered = boolean_or(is_sector_recovered,
				rhs.is_sector_recovered);
			num_recovered_sectors = count(is_sector_recovered);
			num_sectors = is_sector_recovered.size();

			// This is not allowed because we can't add two fractions
			// without knowing their relative contribution to the whole.
			if (flux_fraction_good != -1) {
				throw std::runtime_error("Can't add decoder_stats with "
					"known fraction of fluxes encoded as good!");
			} else {
				flux_fraction_good = rhs.flux_fraction_good;
			}

			return *this;
		}
};

class decoded_tracks {
	public:
		int last_track = 0;
		int last_decoded_sector = 0;

		std::map<IDAM, DAM> sector_data;
		std::map<int, decoder_stats> stats_per_track;

		// Add sectors from another track to this one.
		decoded_tracks& operator+= (const decoded_tracks & rhs) {
			last_track = std::max(last_track, rhs.last_track);
			last_decoded_sector = std::max(last_decoded_sector,
				rhs.last_decoded_sector);

			sector_data.insert(rhs.sector_data.begin(),
				rhs.sector_data.end());

			for (const std::pair<int, decoder_stats>
				track_stat : rhs.stats_per_track) {
				stats_per_track[track_stat.first] += track_stat.second;
			}

			return *this;
		}
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
		// Some of these parameters are pretty hacky. Fix
		// later (TODO).
		decoder_stats decode(timeline & line_to_decode,
			decoded_tracks & out_decoded,
			bool verbose, bool show_stats);

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