#include <iostream>
#include <fstream>
#include <vector>
#include <set>

#include "address_marks.h"
#include "pulse_train.h"
#include "sector_data.h"
#include "rabin_karp.h"
#include "preambles.h"

#include "timeline.h"

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
		void decode(const timeline & line_to_decode);

		// For debugging.
		void dump_to_file(const timeline & line_to_dump,
			std::string data_filename,
			std::string error_filename) const;
};

// This is still in need of some refactoring. sector_data should also have
// a total error count so that we can determine how many errors there are
// in a given chunk decoding. TODO?

void decoder::decode(const timeline & line_to_decode) {
	int CRC_failures = 0;

	bool has_last_admark = false;
	address_mark last_admark;

	std::vector<std::pair<address_mark, address_mark> > full_sectors_found;
	size_t last_start = 0, i = 0;

	for (const timeslice & ts: line_to_decode.timeslices) {
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
				// with this DAM. The "128" here is the minimum length
				// of the data in a DAM.
				if (has_last_admark &&
					last_admark.mark_type == A_IDAM && 
					start - last_admark.byte_stream_index < 128) {

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

				break;
			default: break;
		}

		std::cout << start;
		if (i++ > 0) {
			std::cout << " (+" << start - last_start << ")";
		}
		std::cout << "\t";
		admark.print_info();
		std::cout << "\n";

		last_start = start;
		last_admark = admark;
		has_last_admark = true;
	}

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


/*
def decode_floppy_struct(dec_bytes, last_idam_datalen):
	signature = dec_bytes[:4]
	floppy_info = {}

	# http://dmweb.free.fr/files/Atari-Copy-Protection-V1.4.pdf p. 12
	data_length = [128, 256, 512, 1024]

	# Gracefully handle partially recovered headers.
	if signature == b'\xa1\xa1\xa1\xfe':
		floppy_info["header"] = "IDAM"
		if len(dec_bytes) > 4:	floppy_info["track"] = dec_bytes[4]
		if len(dec_bytes) > 5: floppy_info["head"] = dec_bytes[5]
		if len(dec_bytes) > 6: floppy_info["sector"] = dec_bytes[6]
		if len(dec_bytes) > 7: floppy_info["datalen"] = data_length[dec_bytes[7] & 3]
		if len(dec_bytes) > 9: floppy_info["CRC"] = struct.unpack('>H', dec_bytes[8:10])[0]
		if len(dec_bytes) > 9: floppy_info["CRC_OK"] = crc16(dec_bytes[:8], 0, 8) == floppy_info["CRC"]
		if len(dec_bytes) < 10: floppy_info["truncated"] = True

		return floppy_info

	if signature == b'\xa1\xa1\xa1\xfb':
		floppy_info["header"] = "DAM"

	if signature == b'\xa1\xa1\xa1\xf8':
		floppy_info["header"] = "DAM (deleted)"

	if (signature == b'\xa1\xa1\xa1\xfb') or (signature == b'\xa1\xa1\xa1\xf8'):
		# If this is the first sector, i.e. we don't know the last IDAM
		# data length, recurse on every possible data length and return OK
		# on the first one with CRC OK.
		if last_idam_datalen == 0 or last_idam_datalen is None:
			for i in reversed(data_length):
				speculative_data_chunk = decode_floppy_struct(dec_bytes, i)
				if "CRC_OK" in speculative_data_chunk and \
					speculative_data_chunk["CRC_OK"]:
					return speculative_data_chunk

		# Every attempt gave a bad result. Abort.
		if last_idam_datalen == 0 or last_idam_datalen is None:
			raise IndexError("Could not determine DAM data length!")

		floppy_info["data"] = dec_bytes[4:4+last_idam_datalen]
		if len(dec_bytes) > 4+last_idam_datalen+2: floppy_info["CRC"] = struct.unpack('>H', dec_bytes[4+last_idam_datalen:4+last_idam_datalen+2])[0]
		if len(dec_bytes) > 4+last_idam_datalen+2: floppy_info["CRC_OK"] = crc16(dec_bytes[:4+last_idam_datalen], 0, 4+last_idam_datalen) == floppy_info["CRC"]
		if len(dec_bytes) < 4+last_idam_datalen+3: floppy_info["truncated"] = True

		return floppy_info

	if signature == b'\xc2\xc2\xc2\xfc':
		floppy_info["header"] = "IAM"

		return floppy_info

	print([hex(x) for x in signature])

	raise KeyError
	*/