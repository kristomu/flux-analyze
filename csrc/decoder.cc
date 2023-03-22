#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include <map>
#include <set>

#include "address_marks.h"
#include "pulse_train.h"
#include "sector_data.h"
#include "preambles.h"

#include "timeline.h"
#include "tools.h"

class decoded_tracks {
	public:
		int last_track = 0;
		int last_decoded_sector = 0;

		std::map<IDAM, DAM> sector_data;
};

class decoder {
	private:
		address_mark deserialize(
			std::vector<unsigned char> & raw_bytes,
			size_t byte_stream_start, bool has_last_IDAM,
			const address_mark last_IDAM) const;

	public:
		// Adds OK sectors to the given decoded_tracks
		// structure.
		void decode(timeline & line_to_decode,
			decoded_tracks & decoded);

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

address_mark decoder::deserialize(
	std::vector<unsigned char> & raw_bytes, size_t byte_stream_start,
	bool has_last_IDAM, const address_mark last_IDAM) const {

	address_mark admark;
	admark.set_address_mark_type(raw_bytes);
	admark.byte_stream_index = byte_stream_start;

	size_t datalen; // Data length for DAMs.

	switch(admark.mark_type) {
		case A_IAM: break;
		case A_IDAM:
			admark.idam.set(raw_bytes);
			break;

		// Handle data AMs.
		case A_DAM:
		case A_DDAM:
			// If we have a preceding IDAM and the gap between the
			// IDAM starts and the DAM starts is too small to fit
			// another DAM, then we assume the previous IDAM goes
			// with this DAM.
			if (has_last_IDAM &&
				byte_stream_start - last_IDAM.byte_stream_index <
					DAM::minimum_length()) {

				datalen = last_IDAM.idam.datalen;
				std::cout << "Matching IDAM found" << std::endl;
			} else {
				// TODO: If this fails, set the timeslice to
				// TS_UNKNOWN instead of TS_TRUNCATED or
				// TS_DECODED_BAD, as there could be some other
				// data length that would work.
				std::cout << "Guessing" << std::endl;
				datalen = 512; // seems to be standard for floppies.
			}
			admark.dam.set(raw_bytes, datalen);

			// We've been working on the DAM; if it's actually a
			// DDAM, swap them around.
			if (admark.mark_type == A_DDAM) {
				std::swap(admark.dam, admark.ddam);
			}

			break;
		default: break;
	}

	return admark;
}


// This is still in need of some refactoring. sector_data should also have
// a total error count so that we can determine how many errors there are
// in a given chunk decoding. TODO?

void decoder::decode(timeline & line_to_decode, decoded_tracks & decoded) {
	int failures = 0;

	// Used for linking DAMs to IDAMs to determine what
	// sector a given DAM belongs to.
	bool has_last_IDAM = false;
	address_mark last_IDAM;

	std::vector<std::pair<address_mark, address_mark> > full_sectors_found;
	std::map<slice_id_t, address_mark> address_marks;

	size_t last_start = 0, start = 0, i = 0;
	size_t unknowns = 0;

	for (auto ts_pos = line_to_decode.timeslices.begin();
		ts_pos != line_to_decode.timeslices.end(); ++ts_pos) {

		std::cout << std::endl;
		address_mark admark;
		start = ts_pos->sector_data_begin;

		try {
			admark = deserialize(ts_pos->sec_data.decoded_data,
				ts_pos->sector_data_begin, has_last_IDAM,
				last_IDAM);
		} catch (std::out_of_range & e) {
			// Truncated
			ts_pos->status = TS_TRUNCATED;
			continue;
		}

		if (admark.mark_type == A_IDAM) {
			last_IDAM = admark;
			has_last_IDAM = true;
		}

		// Only if a data mark.
		if ((admark.mark_type == A_DAM || admark.mark_type == A_DDAM)
			&& has_last_IDAM) {
			full_sectors_found.push_back(
				std::pair<address_mark, address_mark>
				(last_IDAM, admark));
		}

		switch (admark.is_OK()) {
			case YES: ts_pos->status = TS_DECODED_OK; break;
			case NO:
				ts_pos->status = TS_DECODED_BAD;
				++failures;
				break;
			case MAYBE: ts_pos->status = TS_DECODED_UNKNOWN; break;
		}

		// Associate the address mark with the current
		// timeslice's ID.
		address_marks[ts_pos->ID] = admark;

		if (admark.mark_type != A_UNKNOWN) {
			std::cout << "Debug: Byte length: " << admark.byte_length() << std::endl;
			size_t admark_MFM_start = ts_pos->sec_data.MFM_train_indices[0];

			if (admark.byte_length() == ts_pos->sec_data.MFM_train_indices.size()) {
				std::cout << "Debug: MFM train indices: " << admark_MFM_start
					<< " and out." << std::endl;
			} else {
				size_t admark_MFM_end = ts_pos->sec_data.MFM_train_indices[admark.byte_length()];
				std::cout << "Debug: MFM train indices: " << admark_MFM_start
					<< " to " << admark_MFM_end << std::endl;
				size_t admark_flux_start = ts_pos->mfm_train.flux_indices[admark_MFM_start],
					admark_flux_end = ts_pos->mfm_train.flux_indices[admark_MFM_end];
				std::cout << "Debug: flux stream indices: " << admark_flux_start
					<< " to " << admark_flux_end << " out of "
					<< ts_pos->flux_data.size() << std::endl;

				// Partition off the part that we still don't know what is, but
				// only if the current chunk was decoded properly, because otherwise
				// insertion or deletions may lead to splitting off too much.

				if (ts_pos->status == TS_DECODED_OK) {
					line_to_decode.split(ts_pos, admark.byte_length(),
						PRESERVE_FIRST);
				}
			}
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
	}

	// Show stats for unknown timeslices.

	for (auto ts_pos = line_to_decode.timeslices.begin();
		ts_pos != line_to_decode.timeslices.end(); ++ts_pos) {

		if (ts_pos->status != TS_UNKNOWN) { continue; }

		std::cout << "Unknown timeslice at " << ts_pos->sector_data_begin << " to " <<
			ts_pos->sector_data_end() << " (" << ts_pos->sector_data_end()-
			ts_pos->sector_data_begin << " bytes.)" << std::endl;
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
		if (idam.is_OK() != YES) {
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

	std::cout << "Address mark failures: " << failures << std::endl;
	std::cout << "Sectors recovered: " << unique_OK_sectors.size() << " out of " << num_sectors << std::endl;
	std::cout << "Unique sector metadata chunks: " << all_sectors.size() << std::endl;
	std::cout << "Total timeslices: " << line_to_decode.timeslices.size() << std::endl;
}

void decoder::dump_all_to_file(const timeline & line_to_dump,
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

void decoder::dump_sector_files(const timeline & line_to_dump,
	std::string prefix) const {

	size_t i = 0;

	for (const timeslice & ts: line_to_dump.timeslices) {
		std::string prefix_ext = prefix + "_" + itos(i) + "_";
		switch(ts.status) {
			case TS_UNKNOWN: prefix_ext += "unknown"; break;
			case TS_TRUNCATED: prefix_ext += "truncated"; break;
			case TS_DECODED_OK: prefix_ext += "decoded_ok"; break;
			case TS_DECODED_BAD: prefix_ext += "decoded_bad"; break;
			case TS_DECODED_UNKNOWN: prefix_ext += "decoded_unknown"; break;
			case TS_PREAMBLE_FOUND: prefix_ext += "preamble_found"; break;
			default: prefix_ext += "status_error"; break;
		}

		std::ofstream data_out(prefix_ext + ".dat",
			std::ios::out | std::ios::binary),
			error_out(prefix_ext + ".mask",
			std::ios::out | std::ios::binary),
			flux_out(prefix_ext + ".flux",
			std::ios::out | std::ios::binary),
			mfm_out(prefix_ext + ".mfm", std::ios::out);

		data_out.write((char *)ts.sec_data.decoded_data.data(),
			ts.sec_data.decoded_data.size());
		error_out.write((char *)ts.sec_data.errors.data(),
			ts.sec_data.errors.size());

		// Implicit buffered copy: this is faster than dumping one
		// int at a time to an ostream.
		std::vector<char> flux_data_bytes;
		std::copy(ts.flux_data.begin(), ts.flux_data.end(),
			std::back_inserter(flux_data_bytes));
		flux_out.write(flux_data_bytes.data(),
			flux_data_bytes.size());

		for (char bit: ts.mfm_train.data) {
			if (bit == 0) {
				mfm_out << "0";
			} else {
				mfm_out << "1";
			}
		}

		++i;
	}
}

// TODO: clean this hack up. Refactor. decoder should be more
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
			// Check if we can't get *any* sectors for this track.
			// If so, we'll just print an error message once instead
			// of spamming the console.
			bool keep_quiet = true;
			for(lookup.sector = 1; lookup.sector <= real_sectors
					&& keep_quiet; ++lookup.sector) {
				auto pos = d_tracks.sector_data.find(lookup);
				keep_quiet &= (pos == d_tracks.sector_data.end());
			}

			if (keep_quiet) {
				std::cout << "Couldn't find any sector. t " <<
					lookup.track << ", h " << lookup.head << std::endl;
			}

			for(lookup.sector = 1; lookup.sector <= real_sectors;
				++lookup.sector) {

				auto pos = d_tracks.sector_data.find(lookup);
				++all_sectors;

				// If nothing was found, then add a blank region
				// to the image and to the mask.
				if (pos == d_tracks.sector_data.end()) {
					if (!keep_quiet) {
						std::cout << "Couldn't find " << lookup.track << ", " << lookup.head << ", " << lookup.sector << std::endl;
					}
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

	image_file.write(image.data(), image.size());
	mask_file.write(mask.data(), mask.size());
}
