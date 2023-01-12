// Standalone FluxEngine IBM decoder with some tricks
// from the Python version... hopefully faster and
// more accurate. Perhaps multiple flux file support ???

#include <algorithm>
#include <iostream>
#include <sqlite3.h>

#include <fstream>
#include <string>
#include <stdexcept>
#include <vector>
#include <cmath>

#include <zlib.h>

#include "rabin_karp.h"
#include "flux_record.h"

// https://stackoverflow.com/a/14612943
int sign(int x) {
	return (x > 0) - (x < 0);
}

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

// Because we want to be robust, we're going to only search for the short
// MFM sequences; we don't make any assumptions about the gap data that
// precedes them.

// The following are given as bit sequences, i.e. *before* MFM decoding.
// Hence we don't have to deal with detecting out-of-band errors.

// // //

// Turn a sequence (e.g. 1 0 1 0 0 1) into a run length list (1 2). Run
// length lists are used for the ordinal search approach.

std::vector<int> sequence_to_run_lengths(const std::vector<char> & sequence) {
	// The sequence must start and end in an 1, because otherwise
	// we don't know if, say "0 1" isn't really an excerpt from
	// "0 0 1"; ditto, we don't know if "1 0" isn't really an excerpt
	// from "1 0 0 ". So first establish this.

	//if (sequence.size() == 0) { return {}; }

	size_t begin, end;

	for (begin = 0; begin < sequence.size() &&
		sequence[begin] != 1; ++begin);

	for (end = sequence.size()-1; end+1 > 0 &&
		sequence[end] != 1; --end);

	std::vector<int> run_lengths;

	int cur_run_length = 0;

	for (size_t i = begin+1; i <= end; ++i) {
		if (sequence[i] == 0) {
			++cur_run_length;
		} else {
			run_lengths.push_back(cur_run_length);
			cur_run_length = 0;
		}
	}

	return run_lengths;
}

class MFM_magic_bytes {
	public:
		std::vector<char> short_A1 = {
			0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1};
		std::vector<char> short_C2 = {
			0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0};
		std::vector<char> A1_sequence, C2_sequence;
		std::vector<int> A1_run_lengths, C2_run_lengths;

		MFM_magic_bytes();
};

MFM_magic_bytes::MFM_magic_bytes() {
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

	A1_run_lengths = sequence_to_run_lengths(A1_sequence);
	C2_run_lengths = sequence_to_run_lengths(C2_sequence);
}

void ordinal_search(const flux_record & flux) {

	// First generate a delta array, which consists of the sign value
	// of adjacent values of the flux delay. If a delay is increasing,
	// then the delta is +1; if it's decreasing, it's -1; if it's the
	// same, it's 0.
	
	std::vector<int> flux_delta(flux.fluxes.size()-1, 0);

	for (size_t i = 0; i < flux_delta.size(); ++i) {
		flux_delta[i] = sign(flux.fluxes[i+1] - flux.fluxes[i]);
		std::cout << flux_delta[i] << " ";
	}

	std::cout << std::endl;
}

// PROCESSING

// Turn a flux delta vector into an MFM pulse train by comparing the distance
// between flux reversals. Our job is made more difficult because apparently
// the mapping is nonlinear: one clock length corresponds to 01, 1.5 clock
// lengths to 001, and 2 clock lengths to 0001.

// For this reason, longer delays are undefined. I may handle these later,
// but note that they are out of spec and should never occur on uncorrupted
// floppies with the right clock.

// TODO: When searching for patterns (i.e. not correcting errors), it may be
// quicker to make a run-length encoding of this (e.g. 10010001 = 2 3, and then
// search for run-length encoded stuff instead, since both the haystack and
// needle will be shorter. Do that perhaps?)

// The error_out double is set to a badness-of-fit value: the higher the worse.
// This might be usable as a very quick clock inference method, but I'll have to
// test more before I know for sure.

// In addition, the nonlinearity might save us sometimes: suppose that we have an
// interval of 2.5 clocks. Then it's likely that this is either one clock followed
// by 1.5 clocks (01001) or the other way around (00101), with some central flux
// reversal having been erased. But I'm not going to do recovery before I've got
// this working on normal images.

// Level one:

std::vector<char> get_MFM_train(double clock,
		const std::vector<int> fluxes, double & error_out) {

	std::vector<char> sequence_bits;
	sequence_bits.reserve(fluxes.size() * 4 / 3); // expected number of bits per reversal
	sequence_bits.push_back(1);

	error_out = 0;

	for (size_t i = 0; i < fluxes.size(); ++i) {
		// We'll model the nonlinearity like this:
		//		- There's always half a clock's delay before anything happens.
		//		- Then one zero corresponds to half a clock more, two zeroes is
		//			two halves more, and three zeroes is three, and so on.

		double half_clocks = (fluxes[i]*2)/clock;

		// Subtract the constant half-clock offset of one, then round the next
		// to get the number of zeroes.
		int zeroes = std::max(0, (int)std::round(half_clocks - 1));

		// Update Euclidean error.
		// Clamp the zeroes to [1..3] to penalize wrong clock guesses.
		// NOTE: This may produce very misleading errors if the disk is
		// partially corrupted because the mean is not a robust estimator.
		// Deal with it later if required.
		int clamped_zeroes = std::max(1, std::min(3, zeroes));
		double error_term = fluxes[i] - (clamped_zeroes + 1) * clock/2.0;
		error_out += error_term * error_term;

		for (int j = 0; j < zeroes; ++j) {
			sequence_bits.push_back(0);
		}
		sequence_bits.push_back(1);
	}

	error_out = std::sqrt(error_out / fluxes.size());

	return sequence_bits;
}

// Level two:

class MFM_data {
	public:
		std::vector<char> decoded_data;
		std::vector<char> errors;

		void decode_MFM(const std::vector<char>::const_iterator & MFM_train_start,
			const std::vector<char>::const_iterator & MFM_train_end);
};

void MFM_data::decode_MFM(const std::vector<char>::const_iterator & MFM_train_start,
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
	// and a 1 otherwise.

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

/* And then the plan is something like:
 * 	- Keep a list of MFM intro markers
 *	- Generate ordinal delta patterns
 *	- Use these as needles and the flux delta as haystack
 *	- For each match, do a narrowing down to determine the clock 
 *	- Try to decode the match. If OK, all's well. Otherwise... */

/* Then later add fancy stuff like warping (PLL), possible lining up
 * multiple incomplete reads - e.g. if #1 is 2 2 8 and #2 is 2 5 4, then
 * #1 = ..X..X.......X
 * #2 = ..X.....X....X
 * which could suggest
 *    = ..X..X..X....X
 * as a true composition. */


int main() {

	test_rabin_karp();

	std::vector<flux_record> flux_records = get_flux_record(
		"../tracks/MS_Plus_disk3_warped_track.flux", true);

	// Some preliminary testing goes here.

	MFM_magic_bytes magic;

	for (const flux_record & f: flux_records) {
		double error;

		std::vector<char> sequence = get_MFM_train(23.5,
				f.fluxes, error);

		std::cout << f.track << ", " << f.side << " error = " << error << "\n";

		// Do a simple search for A1 and C1 patterns.
		rabin_karp rk(magic.A1_sequence);
		rk.add(magic.C2_sequence);

		std::vector<size_t> a1_positions = rk.find_matches(sequence); //search(sequence, magic.A1_sequence);
		std::cout << "Sequences found: " << a1_positions.size() << std::endl;
		std::cout << std::endl;

		MFM_data md;
		md.decode_MFM(sequence.begin() + a1_positions[1], sequence.begin()+a1_positions[2]);

		std::ofstream fout("data.dat", std::ios::out | std::ios::binary);
		fout.write(md.decoded_data.data(), md.decoded_data.size());
		fout.close();
		fout = std::ofstream("errors.dat", std::ios::out | std::ios::binary);
		fout.write(md.errors.data(), md.errors.size());
		fout.close();

	}
	return 0;
}