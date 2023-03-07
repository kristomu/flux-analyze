// Standalone FluxEngine IBM decoder with some tricks
// from the Python version... hopefully faster and
// more accurate. Perhaps multiple flux file support ???

#include <algorithm>
#include <iostream>

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

// I think what I need to do is use R-K to search for pieces from a timeslice
// to find the next timeslice that's like it, because the periodicity isn't
// exact (due to missing or inserted flux transitions), and thus it's going
// to be harder than this... What I really want is to do an approximate
// search of a timeslice to find out where it's repeated.

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
	size_t last_track = 0;

	if (flux_records.empty()) {
		std::cout << "Error: flux database contains no data!" << std::endl;
		return -1;
	}

	for (const flux_record & f: flux_records) {

		std::vector<int> fluxes = f.fluxes;
		decoder IBM_decoder;

		std::cout << "Checking track " << f.track << ", side " << f.side << std::endl;
		last_track = std::max(last_track, (size_t)f.track);

		// Get locations/offsets into the flux where a magic preamble (A1A1A1
		// or C2C2C2) might be found.

		// Using a limited search strategy would cut the runtime to 62% for
		// uncorrupted floppies, but I don't know if it's worth the complexity,
		// so I won't do it (yet).

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

			// If the match was right at the start, then the drive started
			// reading right at the edge of either a preamble or something
			// that looks a lot like it. Because we don't have the previous
			// byte, we can't tell, so skip.
			if (matches[j].match_location == 0) { continue; }

			size_t start_idx = matches[j].match_location - 1,
				end_idx = fluxes.size();

			if (j < matches.size()-1) {
				end_idx = matches[j+1].match_location;
			}

			// HACK: If the matched area is too short for a preamble, then
			// it's a false positive. Signal as such; this can happen with
			// a preamble that's slightly too long, e.g. A1A1A1A1F8 (four
			// instead of three). I have an idea for how to reduce this
			// but we should in any case make a sanity check. Probably in a
			// better way than this, though... TODO
			// I need some way to thread this through the filter_matches...
			// that it should be able to say "no, this must be truncated",
			// so I don't have to rely on every preamble being the same length
			// as I'm doing here.
			if (end_idx - start_idx < preamble_info.get_preamble_by_ID(0).size()) {
				std::cout << "Too short!\n";
				next.status = TS_TRUNCATED;
				continue;
			}
			std::cout << "Found " << m.match_location << " with clock " <<
				m.estimated_clock << " (interval " << start_idx << "-" <<
				end_idx << ")" << std::endl;

			double error;

			next.flux_data = std::vector<int>(f.fluxes.begin() + start_idx,
				f.fluxes.begin() + end_idx);
			next.mfm_train = get_MFM_train(m.estimated_clock,
					next.flux_data, error);

			// HACK HACK: Find the offset from the start of the MFM train
			// to the start of the preamble. (TODO: lots of optimization
			// can be done here. Not yet though.)
			std::vector<search_result> preamble_locations =
				preamble_search.find_matches(next.mfm_train.data, 1);
			if (preamble_locations.empty()) {
				// This happens when the clock estimate is wrong, which
				// shouldn't happen. (We should be signaling errors earlier
				// in ordinal_search.)
				throw std::logic_error("Found preamble but then couldn't!"
					" What's going on?");
			}
			next.status = TS_PREAMBLE_FOUND;
			next.preamble_offset = preamble_locations[0].idx;

			next.sec_data = decode_MFM_train(next.mfm_train,
				next.preamble_offset, next.mfm_train.data.size());

			floppy_line.insert(next);
		}

		IBM_decoder.decode(floppy_line, decoded);
	}

	// Very quick and dirty hard-coded values, fix later. TODO
	decoder IBM_decoder_final;
	IBM_decoder_final.dump_image(decoded, "output",
		std::max((size_t)80, last_track+1), 2, 18, 512);

	return 0;
}