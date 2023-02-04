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

#include "pulse_train.h"
#include "sector_data.cc"

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

		double clock = infer_clock(f.fluxes);

		std::cout << "Guessed clock: " << clock << std::endl;

		std::vector<char> sequence = get_MFM_train(clock,
				f.fluxes, error);

		std::cout << f.track << ", " << f.side << " error = " << error << "\n";

		// Do a simple search for A1 and C1 patterns.
		rabin_karp rk(magic.A1_sequence);
		rk.add(magic.C2_sequence);

		std::vector<size_t> a1_positions = rk.find_matches(sequence); //search(sequence, magic.A1_sequence);
		std::cout << "Sequences found: " << a1_positions.size() << std::endl;
		std::cout << std::endl;

		sector_data sd;
		sd.decode_MFM(sequence.begin() + a1_positions[1], sequence.begin()+a1_positions[2]);

		std::ofstream fout("data.dat", std::ios::out | std::ios::binary);
		fout.write(sd.decoded_data.data(), sd.decoded_data.size());
		fout.close();
		fout = std::ofstream("errors.dat", std::ios::out | std::ios::binary);
		fout.write(sd.errors.data(), sd.errors.size());
		fout.close();

	}
	return 0;
}