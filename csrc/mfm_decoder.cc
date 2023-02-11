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
#include "tools.h"

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

	// Quick and dirty ordinal search experiments
	IBM_preamble pram;
	rabin_karp ordinal_preamble_search(pram.ordinal_A1_sequence);
	ordinal_preamble_search.add(pram.ordinal_C2_sequence);

	for (const flux_record & f: flux_records) {

		// Quick and dirty ordinal search experiments,
		// currently returns just about exactly twice the number
		// of locations as the "direct" way of guessing the
		// clock then searching.
		std::vector<char> ordinal_flux = delta_code(f.fluxes);
		std::vector<size_t> ordinal_locations =
			ordinal_preamble_search.find_matches(ordinal_flux);

		// TODO: Clock clustering to weed out false positives.

		double error;

		double clock = infer_clock(f.fluxes);

		std::cout << "Guessed clock: " << clock << std::endl;

		std::vector<char> sequence = get_MFM_train(clock,
				f.fluxes, error);

		std::cout << f.track << ", " << f.side << " error = " << error << "\n";

		decoder IBM_decoder;

		IBM_decoder.decode(sequence);
		IBM_decoder.dump_to_file(sequence, 1, 2, "data.dat", "errors.dat");
	}

	return 0;
}