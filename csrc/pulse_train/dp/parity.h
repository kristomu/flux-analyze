#pragma once

/*	Detecting what MFM transitions are allowed is pretty easy from the pairs of sequence
	bits alone. Assuming the first pair has already been verified as being part of an
	allowed sequence, we have

		XX 01	is valid (01 encodes for a 1 bit and can be preceded by anything)
		01 00	is valid, but no other prefix is valid
		00 10	is valid, and
		10 10	is valid,

	but no other combination is valid.

	However, in order to make in-band give a penalty on MFM transitions that aren't
	allowed, we need to keep track of when something disallowed happens even as each
	decoded pulse sequence contributes a different number of bits. Consider, for
	instance, an incoming 0001 sequence:

			prior bits    | new bits
			... 01 00 01  | 00 01
			even parity	  | even parity

	This is allowed, but

			prior bits    | new bits
			.... 01 00 1	0 00 1?
			odd parity    | odd parity

	here our bitstream stops at an odd number of bits, so we can't tell from the
	prior bits alone if everything is allowed. There is a trick, though; because
	every (legal) pulse length is at least one long, it produces at least one
	zero bit before the trailing edge. So we must have

			.... 01 00 1	0 00 10

	I'll put more information in a separate parity document, but the upshot
	is that, as long as we have no much-too-long or much-too-short flux
	transitions, then an MFM error occurs *only* at odd parity with a pulse
	length of three.
*/


// Checks if the sequence 1, a zeroes, 1, b zeroes, breaks any rules, when
// the sequence 1, a zeroes has already been tested and can be assumed to
// not have broken any rules. This is done by padding with a one on the left
// and a zero on the right if required. Parity is from the perspective of the
// bits so far, including the trailing edge but not any of the zero bits
// counted by zero_bits nor the present pulse's trailing edge.

// See parity.md for more info.
bool is_mfm_valid(short zero_bits, u_char exclusive_parity) {
	// (Trivia: it actually doesn't matter if we use exclusive or
	// inclusive parity here because any pulse with three zero bits
	// emits four in total, and 4 mod 2 == 0, so parity doesn't
	// change. That said, doing so would be supremely confusing.)

	return !(exclusive_parity == 1 && zero_bits == 3);
}

// Given an inclusive parity, get the corresponding exclusive
// parity.
u_char get_exclusive_parity(short zero_bits,
	u_char inclusive_parity) {

	// The +1 is for the trailing edge's one bit.
	int signed_parity = inclusive_parity - (zero_bits + 1);

	while (signed_parity < 0) {
		signed_parity += 2;
	}

	return signed_parity % 2;
}

// And vice versa.

u_char get_inclusive_parity(short zero_bits,
	u_char exclusive_parity) {

	return (exclusive_parity + zero_bits + 1) % 2;
}
