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

	decoder IBM_decoder;

	for (const flux_record & f: flux_records) {

		std::cout << "Checking track " << f.track << ", side " << f.side << std::endl;

		// Get locations/offsets into the flux where a magic preamble (A1A1A1
		// or C2C2C2) might be found. Currently only A1A1A1 until I get
		// Rabin-Karp to report match type.

		std::vector<char> ordinal_flux = get_delta_coding(f.fluxes);
		std::vector<size_t> ordinal_locations =
			ordinal_preamble_search.find_matches(ordinal_flux);
		std::vector<match_with_clock> matches =
			filter_matches(f.fluxes, ordinal_locations,
				std::vector<char>(pram.A1_sequence.begin()+1,
				pram.A1_sequence.end()));

		for (size_t j = 0; j < matches.size(); ++j) {
			// We want to decode everything from this preamble to the
			// next one.
			match_with_clock m = matches[j];

			// We subtract one because the ordinal search starts from 1,
			// and the magic preamble starts with a zero. Thus there has
			// to be a flux transition before the ordinal match. I would
			// prefer to make this more elegant, but I'm not quite sure how.
			// TODO?
			size_t start_idx = matches[j].match_location - 1,
				end_idx = f.fluxes.size();

			if (j < matches.size()-1) {
				end_idx = matches[j+1].match_location;
			}

			std::cout << "Found " << m.match_location << " with clock " <<
				m.estimated_clock << std::endl;

			double error;

			std::vector<char> sequence = get_MFM_train(m.estimated_clock,
					f.fluxes, start_idx, end_idx, error);

			IBM_decoder.decode(sequence);
			IBM_decoder.dump_to_file(sequence, 1, 2, "data.dat", "errors.dat");
		}
	}

	return 0;
}