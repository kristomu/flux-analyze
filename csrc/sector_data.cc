// Level two.

#include "sector_data.h"

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