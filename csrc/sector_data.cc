#include <iostream>
#include <fstream>
#include <vector>

#include "pulse_train.h"
#include "rabin_karp.h"
#include "preambles.h"

// Level two.

// Higher order formatting for IBM floppies, from
// https://www-user.tu-chemnitz.de/~heha/basteln/PC/usbfloppy/floppy.chm/

// The start of a chunk is marked by a synchronization header. There are
// four types:

// Abbrev. Full name                 Byte seq.  What it is
// IAM	   Index Address Mark        C2C2C2 FC, start of a track
// IDAM	   ID Address Mark           A1A1A1 FE, start of a sector
// DAM	   Data Address Mark         A1A1A1 FB, start of sector data
//		   Deleted Data Address Mark A1A1A1 F8,	start of deleted sector data

// Each address mark is followed by a different header structure,
// though the IAM has no header at all. Note that the byte sequence
// contains a deliberate error to signal out-of-band that it's an address
// mark, not ordinary data.

enum address_mark_t { A_IAM, A_IDAM, A_DAM, A_DDAM, A_UNKNOWN };

// I kind of feel like this class is doing double duty... some of its contents
// ought to be de facto singletons, but I can't quite see how to do it
// properly.
class sector_data {
	public:
		std::vector<unsigned char> decoded_data;
		std::vector<unsigned char> errors;
		std::vector<size_t> MFM_train_indices;

		// limit_length is used to terminate decoding early: if it's set to true,
		// the decoding will stop after max_output_length bytes have been
		// reconstructed. This is used to reset the state between each
		// preamble.
		void decode_and_add_MFM(const MFM_train_data & MFM_train,
			size_t start_index, size_t length);
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

// http://nerdlypleasures.blogspot.com/2015/11/ibm-pc-floppy-disks-deeper-look-at-disk.html
const int MAX_DATALEN_IDX = 6;
const std::array<int, MAX_DATALEN_IDX+1> IDAM_datalen = {
	128, 256, 512, 1024, 2048, 4096, 8192};

class IDAM {		// ID Address Mark
	public:
		// ID Address Mark info. Note these are raw data
		// so e.g. datalen would be an index into a vector. TODO:
		// Abstract that away, either here or in address_mark.

		int track, head, sector, datalen;
		unsigned short CRC;

		// The address mark doesn't include this, but since we're only
		// reading, we just check the validity of the CRC at read time.
		bool CRC_OK = false;
		bool truncated;

		// By convention, the raw bytes start at the IDAM preamble (A1A1A1FE).
		void set(std::vector<unsigned char> & raw_bytes);
};


// Reads off an unsigned value in most significant byte first format.
unsigned int msb_to_int(const std::vector<unsigned char> & data,
	std::vector<unsigned char>::const_iterator pos, int num_bytes) {

	unsigned int out = 0;

	for (int i = 0; i < num_bytes; ++i) {
		if (pos == data.end()) { return out; } // early stopping
		out <<= 8;
		out += (int)*pos++;
	}

	return out;
}

unsigned int msb_to_int(const std::vector<unsigned char> & data,
	size_t idx, int num_bytes) {
	return msb_to_int(data, data.begin()+idx, num_bytes);
}

// https://stackoverflow.com/a/44947877/6183577
// But as this is CCITT-CRC16, the initial CRC value is -1 (0xFFFF), not 0.
unsigned short crc16(std::vector<unsigned char>::const_iterator start,
	std::vector<unsigned char>::const_iterator end) {

	short i, crc = -1;

	for (auto pos = start; pos != end; ++pos) {
		crc ^= (int)*pos << 8;
		i = 8;
		do {
			if (crc & 0x8000) {
				crc = crc << 1 ^ 0x1021;
			} else {
				crc = crc << 1;
			}
		} while (--i);
	}

	return crc;
}

void IDAM::set(std::vector<unsigned char> & raw_bytes) {

	track = raw_bytes[4];
	head = raw_bytes[5];
	sector = raw_bytes[6];

	// Quick and dirty way to deal with out-of-bounds/corrupted values.
	// Also do a lookup to IDAM_datalen to turn this into an actual data chunk
	// length.
	datalen = IDAM_datalen[std::min((int)raw_bytes[7], MAX_DATALEN_IDX)];

	CRC = msb_to_int(raw_bytes, 8, 2);
	CRC_OK = crc16(raw_bytes.begin(), raw_bytes.begin()+8) == CRC;

	truncated = raw_bytes.size() < 10;
}

class DAM {		// Data Address Mark (also used for Deleted Data)
	public:
		std::vector<char> data;
		unsigned short CRC;

		// Auxiliary info
		bool CRC_OK;
		bool deleted; // Is it a DAM or a DDAM?

		void set(std::vector<unsigned char> & raw_bytes, int datalen);
};

void DAM::set(std::vector<unsigned char> & raw_bytes, int datalen) {

	deleted = (raw_bytes[3] == 0xF8);
	data = std::vector<char>(raw_bytes.begin() + 4,
		raw_bytes.begin() + 4 + datalen);

	CRC = msb_to_int(raw_bytes, raw_bytes.begin() + 4 + datalen, 2);
	CRC_OK = crc16(raw_bytes.begin(), raw_bytes.begin() + 4 + datalen) == CRC;

}

class IAM {		// Index Address Mark
	public:
		// Reference offset would go here; the IAM doesn't actually contain
		// any data or metadata.
};

// It's easier to put all the address mark types into one class
// than deal with the problems of reconciling strong typing with a
// list of different types of address marks. This is doubly true because
// a Data Address Mark's metadata is listed in the IDAM that precedes it;
// so we need to know the relative location of DAMs and IDAMs.

class address_mark {
	public:
		address_mark_t mark_type;

		IDAM idam;
		IAM iam;
		DAM dam, ddam;

		// TODO? Some kind of reference of where this was located all
		// the way back to the pulse train, so that we can cross off
		// successfully decoded chunks and work by a process of
		// elimination.

		void set_address_mark(
			const std::vector<unsigned char> & header_bytes);

};

void address_mark::set_address_mark(const std::vector<unsigned char> & header_bytes) {
	mark_type = A_UNKNOWN;
	if(header_bytes == std::vector<u_char>({0xC2, 0xC2, 0xC2, 0xFC})) {
		mark_type = A_IAM;
	}
	if(header_bytes == std::vector<u_char>({0xA1, 0xA1, 0xA1, 0xFE})) {
		mark_type = A_IDAM;
	}
	if(header_bytes == std::vector<u_char>({0xA1, 0xA1, 0xA1, 0xFB})) {
		mark_type = A_DAM;
	}
	if(header_bytes == std::vector<u_char>({0xA1, 0xA1, 0xA1, 0xF8})) {
		mark_type = A_DDAM;
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
		switch(admark.mark_type) {
			case A_IAM:
				std::cout << "Index address mark at " <<
					start;
				if (i > 0) {
					std::cout << " (+" << start - address_mark_starts[i-1] << ")";
				}
				std::cout << "\n";
				break;
			case A_IDAM:
				std::cout << "ID address mark at "  <<
					start;
				if (i > 0) {
					std::cout << " (+" << start - address_mark_starts[i-1] << ")";
				}
				std::cout << "\n";
				
				admark.idam.set(this_chunk);
				std::cout << "\ttrack: " << admark.idam.track << ", head: "
					<< admark.idam.head << ", sector: " << admark.idam.sector
					<< ", data len: " << admark.idam.datalen << " CRC: "
					<< admark.idam.CRC << " ";
				if (admark.idam.CRC_OK) { std::cout << "CRC OK"; } else {
					++ CRC_failures;
					std::cout << "CRC bad";
				}
				std::cout << "\n";
				break;
			case A_DAM:
				std::cout << "Data address mark at "  <<
					start;
				if (i > 0) {
					std::cout << " (+" << start - address_mark_starts[i-1] << ")";
				}
				std::cout << "\n";
				// The stuff below is a giant hack based on that floppy
				// sectors are usually 512 bytes. I really need to
				// mass decode either every address mark region or
				// the whole record instead of decoding piecemeal.
				// Decoding region wise has the benefit that an MFM
				// error will only propagate so far, even in the worst
				// case. I also need to somehow connect DAM and DDAMs
				// with their preceding IDAMs.
				admark.dam.set(this_chunk, 512);
				std::cout << "\tCRC: " << admark.dam.CRC << " ";
				if (admark.dam.CRC_OK) { std::cout << "CRC OK\n"; } else {
					std::cout << "CRC bad\n";
					++ CRC_failures;
				}
				break;
			case A_DDAM:
				std::cout << "Deleted data address mark at "  <<
					start;
				if (i > 0) {
					std::cout << " (+" << start - address_mark_starts[i-1] << ")";
				}
				std::cout << "\n";
				// See above.
				admark.ddam.set(this_chunk, 512);
				std::cout << "\tCRC: " << admark.ddam.CRC << " ";
				if (admark.ddam.CRC_OK) { std::cout << "CRC OK\n"; } else {
					std::cout << "CRC bad\n";
					++ CRC_failures;
				}
				break;
			default:
				std::cout << "???? Something else at " <<
					start << std::endl;
				break;
		}
	}

	std::cout << "CRC failures: " << CRC_failures << std::endl;
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