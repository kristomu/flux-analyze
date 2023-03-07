#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include <map>
#include <set>

#include "address_marks.h"
#include "pulse_train.h"
#include "sector_data.h"
#include "rabin_karp.h"
#include "preambles.h"

#include "timeline.h"

class decoded_tracks {
	public:
		int last_track = 0;
		int last_decoded_sector = 0;

		std::map<IDAM, DAM> sector_data;
};

class decoder {
	private:
		std::vector<address_mark> address_marks;

		// TODO later: something that keeps this cache from going
		// stale if the user inputs one MFM train when doing calc_preamble,
		// and then another with decode.
		bool preamble_locations_calculated;
		std::vector<size_t> preamble_locations;

		// QND fix later: gets an index into the MFM train
		// consistent with the beginning of the idx-th preamble,
		// or the end of the array if idx is too large.
		size_t get_pos_by_idx(
			const std::vector<char> & MFM_train, size_t idx) const {

			if (!preamble_locations_calculated) {
				throw std::runtime_error("Preambles not located yet");
			}

			if (idx < preamble_locations.size()) {
				return preamble_locations[idx];
			} else {
				return MFM_train.size();
			}
		}

	public:
		// Adds OK sectors to the given decoded_tracks
		// structure.
		void decode(timeline & line_to_decode,
			decoded_tracks & decoded);

		// For debugging.
		void dump_to_file(const timeline & line_to_dump,
			std::string data_filename,
			std::string error_filename) const;

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

// This is still in need of some refactoring. sector_data should also have
// a total error count so that we can determine how many errors there are
// in a given chunk decoding. TODO?

void decoder::decode(timeline & line_to_decode, decoded_tracks & decoded) {
	int CRC_failures = 0;

	bool has_last_admark = false;
	address_mark last_admark;

	std::vector<std::pair<address_mark, address_mark> > full_sectors_found;
	size_t last_start = 0, i = 0;
	size_t unknowns = 0;

	for (timeslice & ts: line_to_decode.timeslices) {
		std::cout << std::endl;
		sector_data sd = ts.sec_data;
		// Quick and dirty (very ugly) way of separating out the
		// vector belonging to this chunk. TODO: Fix later: either
		// actually store multiple vectors in the MFM_train (probably
		// best) or change all the signatures.
		auto this_chunk = sd.decoded_data;
		size_t start = ts.sector_data_begin;

		address_mark admark;
		admark.set_address_mark_type(this_chunk);
		admark.byte_stream_index = start;

		size_t datalen;

		// Try to deserialize. (Maybe this
		// should be in address_mark instead? TODO?)
		try {
			ts.status = TS_DECODED_UNKNOWN;

			switch(admark.mark_type) {
				case A_IAM: break;
				case A_IDAM:
					admark.idam.set(this_chunk);
					break;

				// Handle data AMs.
				case A_DAM:
				case A_DDAM:
					// If we have a preceding IDAM and the gap between the
					// IDAM starts and the DAM starts is too small to fit
					// another DAM, then we assume the previous IDAM goes
					// with this DAM.
					if (has_last_admark &&
						last_admark.mark_type == A_IDAM &&
						start - last_admark.byte_stream_index <
							DAM::minimum_length()) {

						datalen = last_admark.idam.datalen;
						std::cout << "Matching IDAM found" << std::endl;
					} else {
						std::cout << "Guessing" << std::endl;
						datalen = 512; // seems to be standard for floppies.
					}
					admark.dam.set(this_chunk, datalen);
					if (!admark.dam.CRC_OK) {
						++ CRC_failures;
					}

					if (has_last_admark && last_admark.mark_type == A_IDAM) {
						full_sectors_found.push_back(
							std::pair<address_mark, address_mark>
							(last_admark, admark));
					}

					// We've been working on the DAM; if it's actually a
					// DDAM, swap them around.
					if (admark.mark_type == A_DDAM) {
						std::swap(admark.dam, admark.ddam);
					}

					break;
				default: break;
			}
			try {
				if (admark.is_OK()) {
					ts.status = TS_DECODED_OK;
				} else {
					ts.status = TS_DECODED_BAD;
				}
			// Unknown address mark, can't tell.
			} catch (std::runtime_error & e) {}
		} catch (std::out_of_range & e) {
			ts.status = TS_TRUNCATED;
			// The address mark is truncated, so just skip it.
			continue;
		}

		if (admark.mark_type != A_UNKNOWN) {
			std::cout << "Debug: Byte length: " << admark.byte_length() << std::endl;
			size_t admark_MFM_start = ts.sec_data.MFM_train_indices[0],
				admark_MFM_end = ts.sec_data.MFM_train_indices[admark.byte_length()];
			std::cout << "Debug: MFM train indices: " << admark_MFM_start <<
				" to " << admark_MFM_end << std::endl;
			size_t admark_flux_start = ts.mfm_train.flux_indices[admark_MFM_start],
				admark_flux_end = ts.mfm_train.flux_indices[admark_MFM_end];
			std::cout << "Debug: flux stream indices: " << admark_flux_start <<
				" to " << admark_flux_end << " out of " << ts.flux_data.size() << std::endl;
		} else {
			++unknowns;
		}

		std::cout << start;
		if (i++ > 0) {
			std::cout << " (+" << start - last_start << ") ";
		}
		std::cout << "\t";
		admark.print_info();
		std::cout << "\n";

		last_start = start;
		last_admark = admark;
		has_last_admark = true;
	}

	std::cout << "Unknown AMs detected: " << unknowns << std::endl;

	// This set contains the IDAMs corresponding to sectors with
	// valid data.
	std::set<IDAM> unique_OK_sectors, all_sectors;

	// Go through all the found sectors to count them.
	for (auto & sector_constituents: full_sectors_found) {

		if (sector_constituents.second.mark_type == A_DDAM) {
			std::cout << "Warning: found DDAM. This is not "
				"currently supported." << std::endl;
			continue;
		}

		IDAM idam = sector_constituents.first.idam;
		DAM dam = sector_constituents.second.dam;

		if (unique_OK_sectors.find(idam) != unique_OK_sectors.end()) {
			continue; // already seen it
		}

		// Maybe also check for MFM errors??? It's a weak check but should
		// discard stuff that just happens to have a colliding bad CRC...
		if (idam.CRC_OK && dam.CRC_OK) {
			unique_OK_sectors.insert(idam);

			// Add to the decoded structure.
			// TODO??? add IDAMs with blank DAMs if the IDAM was
			// found but the data wasn't???

			decoded.sector_data[idam] = dam;
			decoded.last_track = std::max(decoded.last_track,
				idam.track);
			decoded.last_decoded_sector = std::max(
				decoded.last_decoded_sector, idam.sector);
		}

		// Don't insert IDAMs with bad CRC; their sector metadata could be
		// scrambled and refer to something that doesn't exist.
		if (idam.CRC_OK) {
			all_sectors.insert(idam);
		}
	}

	// Guess at the number of sectors, assuming the sectors start at 1.
	size_t num_sectors = all_sectors.size();

	for (const IDAM & idam: unique_OK_sectors) {
		num_sectors = std::max(num_sectors, (size_t) idam.sector);

		std::cout << "OK sector found: ";
		idam.print_info();
		std::cout << "\n";
	}

	std::cout << "CRC failures: " << CRC_failures << std::endl;
	std::cout << "Sectors recovered: " << unique_OK_sectors.size() << " out of " << num_sectors << std::endl;
	std::cout << "Unique sector metadata chunks: " << all_sectors.size() << std::endl;
	std::cout << "Total timeslices: " << line_to_decode.timeslices.size() << std::endl;
}

void decoder::dump_to_file(const timeline & line_to_dump,
	std::string data_filename, std::string error_filename) const {

	std::ofstream data_out(data_filename, std::ios::out | std::ios::binary),
		error_out(error_filename, std::ios::out | std::ios::binary);

	// This is not technically true, due to the way timeslices' sector data are
	// aligned to preambles. However, the "right" way to do it (just dump all the
	// bits one after the other) would make it much harder to see data in the dump
	// files. So fortuitiously, the side effect that sector datas are byte-aligned
	// with preambles serves us well!
	for (const timeslice & ts: line_to_dump.timeslices) {
		data_out.write((char *)ts.sec_data.decoded_data.data(),
			ts.sec_data.decoded_data.size());
		error_out.write((char *)ts.sec_data.errors.data(),
			ts.sec_data.errors.size());
	}
	data_out.close();
	error_out.close();
}

// remove this ugly hack pronto. Refactor. decoder should be more
// persistent and supply stats like these on demand.
void decoder::dump_image(const decoded_tracks & d_tracks,
	std::string file_prefix, int tracks, int heads,
	int sectors_per_track, int default_sector_size) const {

	// These values are hard-coded for IBM floppies for now, e.g.
	// sectors starting at 1. Fix this later if required.

	std::vector<char> image;
	std::vector<char> mask;
	IDAM lookup;

	int real_sectors = (int)std::max(sectors_per_track,
		d_tracks.last_decoded_sector);

	size_t OK_sectors = 0, all_sectors = 0;

	for (lookup.track = 0; lookup.track < tracks; ++lookup.track) {
		for (lookup.head = 0; lookup.head < heads; ++lookup.head) {
			for(lookup.sector = 1; lookup.sector <= real_sectors;
					++lookup.sector) {
				auto pos = d_tracks.sector_data.find(lookup);

				++all_sectors;

				// If nothing was found, then add a blank region
				// to the image and to the mask.
				if (pos == d_tracks.sector_data.end()) {
					std::cout << "Couldn't find " << lookup.track << ", " << lookup.head << ", " << lookup.sector << std::endl;
					std::vector<char> blank(default_sector_size, 0);
					std::copy(blank.begin(), blank.end(),
						std::back_inserter(image));
					std::copy(blank.begin(), blank.end(),
						std::back_inserter(mask));
					continue;
				}

				++OK_sectors;

				std::vector<char> OK_mask(pos->second.data.size(), 0xFF);
				std::copy(pos->second.data.begin(), pos->second.data.end(),
					std::back_inserter(image));
				std::copy(OK_mask.begin(), OK_mask.end(),
					std::back_inserter(mask));

			}
		}
	}

	std::ofstream image_file(file_prefix + ".img"),
		mask_file(file_prefix + ".mask");

	std::copy(image.begin(), image.end(),
		std::ostream_iterator<char>(image_file));
	std::copy(mask.begin(), mask.end(),
		std::ostream_iterator<char>(mask_file));
}
