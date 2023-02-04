#include <iostream>
#include <vector>

#include "rabin_karp.h"

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

// To find these magic preambles, which we need to find the address marks,
// let's first define some auxiliary classes.

// Because we want to be robust, we're going to only search for the short
// MFM sequences; we don't make any assumptions about the gap data that
// precedes them, even though real world floppies usually have padding that
// helps the phase-locked loop lock on.

// The following are given as bit sequences, i.e. *before* MFM decoding.
// Hence we don't have to deal with detecting out-of-band errors.

class IBM_preamble {
	public:
		std::vector<char> short_A1 = {
			0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1};
		std::vector<char> short_C2 = {
			0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0};
		std::vector<char> A1_sequence, C2_sequence;

		IBM_preamble();
};

IBM_preamble::IBM_preamble() {
	// The A1 sequence is the short A1 repeated three times,
	// and ditto for the C2 sequence.

	A1_sequence.reserve(3 * short_A1.size());
	C2_sequence.reserve(3 * short_C2.size());
	for (int i = 0; i < 3; ++i) {
		std::copy(short_A1.begin(), short_A1.end(),
			std::back_inserter(A1_sequence));
		std::copy(short_C2.begin(), short_C2.end(),
			std::back_inserter(C2_sequence));
	}
}

// I kind of feel like this class is doing double duty... some of its contents
// ought to be de facto singletons, but I can't quite see how to do it
// properly.
class sector_data {
	public:
		std::vector<unsigned char> decoded_data;
		std::vector<unsigned char> errors;

		// limit_length is used to terminate decoding early: if it's set to true,
		// the decoding will stop after max_output_length bytes have been
		// reconstructed. This is used for preamble checks.
		void decode_MFM(const std::vector<char>::const_iterator & MFM_train_start,
			const std::vector<char>::const_iterator & MFM_train_end,
			bool limit_length, size_t max_output_length);

		void decode_MFM(const std::vector<char>::const_iterator & MFM_train_start,
			const std::vector<char>::const_iterator & MFM_train_end) {

			decode_MFM(MFM_train_start, MFM_train_end, false, 0);
		}
};

void sector_data::decode_MFM(const std::vector<char>::const_iterator & MFM_train_start,
	const std::vector<char>::const_iterator & MFM_train_end,
	bool limit_length, size_t max_output_length) {

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

	decoded_data.clear();
	errors.clear();
	unsigned char current_char = 0, current_char_errors = 0;
	size_t bits_output = 0;
	auto pos = MFM_train_start;

	// We don't know the MFM train data prior to the first bit,
	// so we don't know if the "last" decoded bit would've been a
	// zero or a one. Therefore, be optimistic and never report an
	// error for the first bit.
	bool beginning = true;
	char last_bit = 0;

	char bits[2];

	while (pos != MFM_train_end) {
		if (limit_length && decoded_data.size() >= max_output_length) {
			return;
		}
		// Clock two bits.
		bits[0] = *pos++;
		if (pos == MFM_train_end) { continue; }
		bits[1] = *pos++;

		// Handle RR. There's a problem here that RR could always
		// be NR, so it's kind of an unsatisfactory answer to the
		// problem. We'd need a nondeterministic automaton
		// or something...
		int clock_pair = bits[0] * 2 + bits[1];
		if (clock_pair == 3) {
			if (last_bit == 0) {
				clock_pair = 2; // Make it RN
			} else {
				clock_pair = 1; // Make it NR.
			}
		}

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

class IDAM {		// ID Address Mark
	public:
		int track, head, sector, datalen;
		unsigned short CRC;

		bool truncated;

		// TODO? Some kind of reference of where this was located all the way back to
		// the pulse train, so that we can work by a process of elimination.
};

class DAM {		// Data Address Mark
	public:
		std::vector<char> data;
		unsigned short CRC;

		// Auxiliary info
		bool CRC_OK;
};

class IAM {		// Index Address Mark
	public:
		// Reference offset would go here; the IAM doesn't actually contain
		// any data or metadata.
};

class decoder {
	private:
		IBM_preamble preambles;
		rabin_karp preamble_search;

		std::vector<DAM> DAMs;
		std::vector<IDAM> IDAMs;
		std::vector<IAM> IAMs;

	public:
		decoder();

		std::vector<size_t> get_preamble_locations(
			std::vector<char> & MFM_train) const;

		void decode(std::vector<char> & MFM_train);
};

// Relying on something initialized by the constructor of one object
// in the constructor of another is kind of icky, but it seems to work.
decoder::decoder() : preambles(), preamble_search(preambles.A1_sequence) {
	preamble_search.add(preambles.C2_sequence);
}


std::vector<size_t> decoder::get_preamble_locations(
	std::vector<char> & MFM_train) const {

	return preamble_search.find_matches(MFM_train);
}

void decoder::decode(std::vector<char> & MFM_train) {
	// First get the preamble locations.
	std::vector<size_t> preamble_locations =
		get_preamble_locations(MFM_train);

	sector_data sd;

	for (size_t i = 0; i < preamble_locations.size(); ++i) {
		std::vector<char>::const_iterator pos = MFM_train.begin() +
			preamble_locations[i];
		std::vector<char>::const_iterator next_pos;
		if (i < preamble_locations.size()-1) {
			next_pos = MFM_train.begin() + preamble_locations[i+1];
		} else {
			next_pos = MFM_train.end();
		}

		// Decode four bytes to determine what mark we're dealing
		// with. (TODO? Iterators? but then what about error/OOB signaling?)
		sd.decode_MFM(pos, next_pos, true, 4);

		// TODO: Move into IBM_preamble, make it return an enum so we can use
		// switch/case or something nicer.
		bool identified = false;
		if(sd.decoded_data == std::vector<u_char>({0xC2, 0xC2, 0xC2, 0xFC})) {
			std::cout << "Index address mark at " <<
				preamble_locations[i] << std::endl;
			identified = true;
		}
		if(sd.decoded_data == std::vector<u_char>({0xA1, 0xA1, 0xA1, 0xFE})) {
			std::cout << "ID address mark at " <<
				preamble_locations[i] << std::endl;
			identified = true;
		}
		if(sd.decoded_data == std::vector<u_char>({0xA1, 0xA1, 0xA1, 0xFB})) {
			std::cout << "Data address mark at " <<
				preamble_locations[i] << std::endl;
			identified = true;
		}
		if(sd.decoded_data == std::vector<u_char>({0xA1, 0xA1, 0xA1, 0xF8})) {
			std::cout << "Deleted data address mark at " <<
				preamble_locations[i] << std::endl;
			identified = true;
		}
		if (!identified) {
			std::cout << "???? Something else at " <<
				preamble_locations[i] << std::endl;
		}
	}
	// Then decode four bytes to determine the mark
	// we're dealing with.
	// Then decode enough bytes to get the header,
	// then enough to get the data if any.
	// We should probably carry through index markers for this, but
	// it's pretty hard to do right...
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