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
#include "ordinal_search.h"
#include "sector_data.h"
#include "tools.h"

#include "timeline.h"
#include "decoder.cc"

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
	IBM_preamble preamble_info;
	rabin_karp ordinal_preamble_search(
		preamble_info.ordinal_A1_sequence, PREAMBLE_ID_A1);
	ordinal_preamble_search.add(
		preamble_info.ordinal_C2_sequence, PREAMBLE_ID_C2);

	rabin_karp preamble_search(preamble_info.A1_sequence,
		PREAMBLE_ID_A1);
	preamble_search.add(preamble_info.C2_sequence,
		PREAMBLE_ID_C2);

	for (const flux_record & f: flux_records) {

		std::vector<int> fluxes = f.fluxes;
		decoder IBM_decoder;

		std::cout << "Checking track " << f.track << ", side " << f.side << std::endl;

		// Get locations/offsets into the flux where a magic preamble (A1A1A1
		// or C2C2C2) might be found.

		std::vector<char> ordinal_flux = get_delta_coding(fluxes);
		std::vector<search_result> ordinal_locations =
			ordinal_preamble_search.find_matches(ordinal_flux);
		std::vector<match_with_clock> matches =
			filter_matches(fluxes, ordinal_locations,
				preamble_info);

		timeline floppy_line;
		// A linear sequence made up of each decoded chunk concatenated
		// in order.
		timeslice next;

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
				end_idx = fluxes.size();

			if (j < matches.size()-1) {
				end_idx = matches[j+1].match_location;
			}

			std::cout << "Found " << m.match_location << " with clock " <<
				m.estimated_clock << std::endl;

			double error;

			next.flux_data = std::vector<int>(f.fluxes.begin() + start_idx,
				f.fluxes.begin() + end_idx);
			next.mfm_train = get_MFM_train(m.estimated_clock,
					next.flux_data, error);

			// HACK HACK: Find the offset from the start of the MFM train
			// to the start of the preamble. (TODO: lots of optimization
			// can be done here. Not yet though.)
			std::vector<search_result> preamble_locations =
				preamble_search.find_matches(next.mfm_train.data);
			if (preamble_locations.empty()) {
				throw std::logic_error("Found preamble but then couldn't!"
					" What's going on?");
			}
			next.has_preamble = true;
			next.preamble_offset = preamble_locations[0].idx;

			next.sec_data = decode_MFM_train(next.mfm_train,
				next.preamble_offset, next.mfm_train.data.size());

			floppy_line.insert(next);
		}

		IBM_decoder.decode(floppy_line);
		IBM_decoder.dump_to_file(floppy_line, "data.dat", "errors.dat");
	}

	return 0;
}