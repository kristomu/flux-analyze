#include <iostream>
#include <fstream>
#include <vector>
#include <set>

#include "address_marks.h"
#include "pulse_train.h"
#include "rabin_karp.h"
#include "preambles.h"

// Level two.

// Higher order formatting for IBM floppies, from
// https://www-user.tu-chemnitz.de/~heha/basteln/PC/usbfloppy/floppy.chm/

// I kind of feel like this class is doing double duty... some of its contents
// ought to be de facto singletons, but I can't quite see how to do it
// properly.
class sector_data {
	public:
		std::vector<unsigned char> decoded_data;
		std::vector<unsigned char> errors;
		std::vector<size_t> MFM_train_indices;

		// Since we want to reset the MFM state at every preamble so that
		// we're guaranteed that errors don't propagate indefinitely, we
		// decode one chunk at a time (from one preamble to the next).
		// Hence start_index and end_index, which give the boundaries of
		// these chunks.
		void decode_and_add_MFM(const MFM_train_data & MFM_train,
			size_t start_index, size_t end_index);
};

void sector_data::decode_and_add_MFM(const MFM_train_data & MFM_train,
			size_t start_index, size_t end_index) {

	// If N denotes "no reversal" and R denotes "reversal", then
	// the MFM state machine is as follows:
	//		1		is represented by 		NR
	//		0		is represented by		RN		if the last bit was zero
	//		0		is represented by		NN		if the last bit was one

	// everything else is an error, though (IIRC) the preambles depend
	// on RN and NN always being decoded to a 0. RR is evidence of either
	// the wrong clock or a spurious magnetic flux change. For lack of
	// anything better, I'll change it to a 0 if the last bit was zero,
	// and a 1 otherwise. (I need a better idea of what kind of corruption
	// may happen.)

	// I'm using 1 and 0 literals for one and zero bits; note that this
	// is exactly opposite of the C tradition (0 is true, 1 is false).

	unsigned char current_char = 0, current_char_errors = 0;
	size_t bits_output = 0;

	// We don't know the MFM train data prior to the first bit,
	// so we don't know if the "last" decoded bit would've been a
	// zero or a one. Therefore, be optimistic and never report an
	// error for the first bit.
	bool beginning = true;
	char last_bit = 0;
	auto pos = MFM_train.data.begin() + start_index,
		MFM_train_end = MFM_train.data.begin() + end_index;
	size_t printed_char_began_here = start_index;

	char bits[2];

	if (end_index - start_index > MFM_train.data.size()) {
		throw std::runtime_error("decode_and_add_MFM: invalid selection"
			" of MFM_train bits to decode!");
	}

	while (pos != MFM_train_end && pos != MFM_train.data.end()) {
		// Clock two bits.
		bits[0] = *pos++;
		if (pos == MFM_train_end) { continue; }
		bits[1] = *pos++;

		// Handle RR by aborting. This sequence should never be emitted
		// from the MFM train decoder, and thus it showing up here
		// would be an error.
		int clock_pair = bits[0] * 2 + bits[1];

		switch(clock_pair) {
			case 0: // NN
				if (!beginning && last_bit != 1) {
					++current_char_errors;
				}
				last_bit = 0;
				break;
			case 1: // NR
				++current_char;
				last_bit = 1;
				break;
			case 2: // RN
				if (!beginning && last_bit != 0) {
					++current_char_errors;
				}
				last_bit = 0;
				break;
			case 3: // RR
				// This shouldn't happen.
				throw std::runtime_error("RR stabilization failed!");
				break;
		}
		bits_output++;
		beginning = false;
		if (bits_output == 8) {
			// Add to the index list the first MFM train bit
			// that contributed to this char, and set it for
			// the next char to add next time around.
			MFM_train_indices.push_back(printed_char_began_here);
			printed_char_began_here = (pos - MFM_train.data.begin()) - start_index;

			// Add data and errors, and reset.
			decoded_data.push_back(current_char);
			errors.push_back(current_char_errors);
			current_char = 0;
			current_char_errors = 0;
			bits_output = 0;
		} else {
			current_char <<= 1;
			current_char_errors <<= 1;
		}
	}
}

class decoder {
	private:
		IBM_preamble preambles;
		rabin_karp preamble_search;

		std::vector<address_mark> address_marks;
		sector_data sd;

		// TODO later: something that keeps this cache from going
		// stale if the user inputs one MFM train when doing calc_preamble,
		// and then another with decode.
		bool preamble_locations_calculated;
		std::vector<size_t> preamble_locations;

		std::vector<size_t> get_preamble_locations(
			std::vector<char> & MFM_train) const;

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
		decoder();

		void calc_preamble_locations(const MFM_train_data & MFM_train);
		void decode(const MFM_train_data & MFM_train);

		// For debugging.
		void dump_to_file(const MFM_train_data & MFM_train,
			std::string data_filename,
			std::string error_filename) const;
};

// Relying on something initialized by the constructor of one object
// in the constructor of another is kind of icky, but it seems to work.
decoder::decoder() : preambles(), preamble_search(preambles.A1_sequence) {
	preamble_search.add(preambles.C2_sequence);

	preamble_locations_calculated = false;
}


void decoder::calc_preamble_locations(
	const MFM_train_data & MFM_train) {

	preamble_locations_calculated = true;
	preamble_locations = preamble_search.find_matches(
		MFM_train.data);
}

// This is getting very ugly and ripe for some serious refactoring.
// But how?
// Repeatedly calling decode_MFM is kind of ugly too, though on the other
// hand, I can spare the decoding cost.

// Decoding all the MFM stuff in one go is appealing because the decoder's
// state (i.e. whether the last bit was 0 or 1) would then carry through.
// On the other hand, decoding only what we need makes the mechanism more
// flexible so that we could repair the chunks that fail to decode while
// leaving all the other chunks alone. I don't think it matters much either
// way; I just have to make a decision.

// Looking at how Python did it may be instructive: it basically batch
// decodes and then has one vector for errors and one for data. This would
// also make it easy to count errors where there shouldn't be any (e.g. inside
// data blocks).

void decoder::decode(const MFM_train_data & MFM_train) {
	// First get the preamble locations. (Memoize)
	if (!preamble_locations_calculated) {
		calc_preamble_locations(MFM_train);
	}

	sd = sector_data(); // clear the sector data.

	std::cout << "Sequences found: " << preamble_locations.size() << std::endl;

	int CRC_failures = 0;

	// Then decode four bytes to determine the mark
	// we're dealing with.
	// Then decode enough bytes to get the header,
	// then enough to get the data if any.
	// We should probably carry through index markers for this, but
	// it's pretty hard to do right...

	// Decode the whole thing, collecting indices into the bytestream
	// on the way.

	std::vector<size_t> address_mark_starts;
	size_t i;

	for (i = 0; i < preamble_locations.size(); ++i) {
		// TODO: More coherent comments.
		// Add the start position of the sector_data data, because
		// that's where the chunk belonging to this address mark
		// starts.
		address_mark_starts.push_back(sd.decoded_data.size());


		size_t cur_MFM_idx = get_pos_by_idx(MFM_train.data, i), 
			next_MFM_idx = get_pos_by_idx(MFM_train.data, i+1);

		// Add this chunk to the decoding.
		sd.decode_and_add_MFM(MFM_train, cur_MFM_idx, next_MFM_idx);
	}

	// Not entirely true but it's easier to avoid off-by-ones this way.
	address_mark_starts.push_back(sd.decoded_data.size());

	bool has_last_admark = false;
	address_mark last_admark;

	std::vector<std::pair<address_mark, address_mark> > full_sectors_found;

	for (i = 0; i < address_mark_starts.size()-1; ++i) {
		size_t start = address_mark_starts[i],
			next = address_mark_starts[i+1];
		// Quick and dirty (very ugly) way of separating out the
		// vector belonging to this chunk. TODO: Fix later: either
		// actually store multiple vectors in the MFM_train (probably
		// best) or change all the signatures.
		std::vector<unsigned char> this_chunk(&sd.decoded_data[start],
			&sd.decoded_data[next]);
		std::vector<unsigned char> preamble(&sd.decoded_data[start],
			&sd.decoded_data[start]+4);

		address_mark admark;
		admark.set_address_mark(preamble);
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
		if (i > 0) {
			std::cout << " (+" << start - address_mark_starts[i-1] << ")";
		}
		std::cout << "\t";
		admark.print_info();
		std::cout << "\n";

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
	std::cout << "Total preambles: " << preamble_locations.size() << std::endl;
}

void decoder::dump_to_file(const MFM_train_data & MFM_train,
	std::string data_filename, std::string error_filename) const {

	if (!preamble_locations_calculated) {
		throw std::runtime_error("Preambles not located yet");
	}

	std::ofstream fout(data_filename, std::ios::out | std::ios::binary);
	fout.write((char *)sd.decoded_data.data(), sd.decoded_data.size());
	fout.close();
	fout = std::ofstream(error_filename, std::ios::out | std::ios::binary);
	fout.write((char *)sd.errors.data(), sd.errors.size());
	fout.close();
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