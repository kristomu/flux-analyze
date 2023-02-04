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

#include "flux_record.h"

#include "pulse_train.h"
#include "sector_data.cc"

// https://stackoverflow.com/a/14612943
int sign(int x) {
	return (x > 0) - (x < 0);
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

	for (const flux_record & f: flux_records) {
		double error;

		double clock = infer_clock(f.fluxes);

		std::cout << "Guessed clock: " << clock << std::endl;

		std::vector<char> sequence = get_MFM_train(clock,
				f.fluxes, error);

		std::cout << f.track << ", " << f.side << " error = " << error << "\n";

		decoder IBM_decoder;

		std::vector<size_t> a1_positions = IBM_decoder.get_preamble_locations(sequence);
		std::cout << "Sequences found: " << a1_positions.size() << std::endl;
		std::cout << std::endl;

		sector_data sd;
		sd.decode_MFM(sequence.begin() + a1_positions[1], sequence.begin()+a1_positions[2]);

		std::ofstream fout("data.dat", std::ios::out | std::ios::binary);
		fout.write((char *)sd.decoded_data.data(), sd.decoded_data.size());
		fout.close();
		fout = std::ofstream("errors.dat", std::ios::out | std::ios::binary);
		fout.write((char *)sd.errors.data(), sd.errors.size());
		fout.close();

		IBM_decoder.decode(sequence);

	}
	return 0;
}