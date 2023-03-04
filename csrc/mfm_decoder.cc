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


// Determine the approximate period of the flux signal (i.e. recordings
// from consecutive passes around the floppy disk). This works by taking
// snippets from the flux record, then doing a Rabin-Karp search, and
// determining the distance between adjacent hits, but in a way that
// doesn't make outliers or unusually common sequences skew the result.

// We currently need to cast everything to char because rabin_karp only
// does char. TODO: Fix that later with templating.

// Reasonable results for a HD floppy should be something around
// 100k-133k, because a clock is 1.5 us and the standard IBM spin rate
// is 300 RPM, i.e. 0.2s per rotation. The theoretical period is 67k if
// all three-clocks, 88k if all two-clocks, but those are very improbable.

// It works OK on good floppies, but returns an absolutely bizarre result
// for RETRY-77-4-t69.0.flux.

double find_approximate_period(const std::vector<int> & in_fluxes) {
	int snippet_length = 10;
	int snippets_to_check = 360;

	// Any snippet that has fewer than this number of hits
	// won't be taken into account for the period calculation.
	size_t hit_threshold = 5;

	// Anything with this short an estimated period is also
	// discarded for being obviously wrong.
	size_t period_threshold = 128;

	// Ugly casting.
	std::vector<char> fluxes(in_fluxes.begin(), in_fluxes.end());

	rabin_karp snippet_search(snippet_length);
	size_t snippet_num = 0, i = 0;

	for (i = 0; i < fluxes.size();
			i += fluxes.size()/snippets_to_check) {

		std::vector<char> flux_snippet(fluxes.begin() + i,
			fluxes.begin() + i + snippet_length);
		// The snippets may match some earlier snippets, so we need
		// to escape the error that would cause.
		try {
			snippet_search.add(flux_snippet, snippet_num);
		} catch (std::invalid_argument & e) {
			continue;
		}
		++snippet_num;
	}

	// Get the hits and dump their locations into a vector.

	std::vector<std::vector<size_t> > snippet_hits(snippet_num,
		std::vector<size_t>());

	std::vector<search_result> results = snippet_search.find_matches(fluxes);
	for (const search_result & result: results) {
		snippet_hits[result.ID].push_back(result.idx);
	}

	// Extract the median period estimates for each snippet.
	std::vector<double> median_periods;

	for (const std::vector<size_t> & hits_this_snippet : snippet_hits) {
		if (hits_this_snippet.size() < hit_threshold) {
			continue;
		}

		std::vector<int> periods;
		for (i = 1; i < hits_this_snippet.size(); ++i) {
			// The distance between two occurrences of the same snippet
			// produces an estimate of the period.
			periods.push_back(hits_this_snippet[i] -
				hits_this_snippet[i-1]);
		}
		double median_period_est = median(periods);
		if (median_period_est > period_threshold) {
			median_periods.push_back(median_period_est);
		}
	}

	if (median_periods.empty()) {
		throw std::runtime_error("Couldn't guess period.");
	}

	return median(median_periods);
}

int main(int argc, char ** argv) {

	test_rabin_karp();

	std::string flux_filename;

	if (argc < 2) {
		flux_filename = "../tracks/MS_Plus_disk3_warped_track.flux";

		std::cout << "Flux file to analyze was not specified.\n";
		std::cout << "Choosing " << flux_filename << " as default.\n";
	} else {
		flux_filename = argv[1];
		std::cout << "Analyzing flux file " << flux_filename << "\n";
	}

	std::vector<flux_record> flux_records = 
		get_flux_record(flux_filename, true);

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

	decoded_tracks decoded;

	for (const flux_record & f: flux_records) {

		std::vector<int> fluxes = f.fluxes;
		decoder IBM_decoder;

		std::cout << "Checking track " << f.track << ", side " << f.side << std::endl;
		try {
			double flux_period = find_approximate_period(fluxes);
			std::cout << "Estimated period: " << flux_period << std::endl;
		} catch (std::runtime_error & e) {
			std::cout << "Estimated period: [UNKNOWN]" << std::endl;
		}

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
				// This happens when the clock estimate is wrong.
				throw std::logic_error("Found preamble but then couldn't!"
					" What's going on?");
			}
			next.has_preamble = true;
			next.preamble_offset = preamble_locations[0].idx;

			next.sec_data = decode_MFM_train(next.mfm_train,
				next.preamble_offset, next.mfm_train.data.size());

			floppy_line.insert(next);
		}

		IBM_decoder.decode(floppy_line, decoded);
		IBM_decoder.dump_to_file(floppy_line, "data.dat", "errors.dat");
	}

	// Very quick and dirty hard-coded values, fix later. TODO
	decoder IBM_decoder_final;
	IBM_decoder_final.dump_image(decoded, "output", 80, 18, 512);

	return 0;
}