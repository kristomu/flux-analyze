#include <iostream>
#include <vector>

// Level two.

class sector_data {
	public:
		std::vector<char> decoded_data;
		std::vector<char> errors;

		void decode_MFM(const std::vector<char>::const_iterator & MFM_train_start,
			const std::vector<char>::const_iterator & MFM_train_end);
};

void sector_data::decode_MFM(const std::vector<char>::const_iterator & MFM_train_start,
	const std::vector<char>::const_iterator & MFM_train_end) {

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