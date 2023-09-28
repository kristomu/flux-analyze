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

// Perform brute-force adaptive (PLL/dewarping) decoding of flux streams.
// This returns the decoded tracks structure for the timeline.

decoded_tracks decode_brute_dewarp(timeline & floppy_line,
	double min_alpha, double max_alpha, double stepsize,
	bool verbose) {

	// Check for any faulty AMs and try to encode them with PLL mode
	// enabled instead. TODO: Refine this into a proper strategy approach.
	// will probably need to involve splitting decode() into something
	// that individually handles timeslices and something that sticks
	// them together, but then how will we keep the sector chunks
	// synchronized with the timeline? Might need a redesign. XXX.

	// ff_track18 shows that dewarping can sometimes make it worse
	// (going from DECODED_BAD to UNKNOWN). We need to keep the old one
	// in that case, which complicates matters considerably. We do that
	// by storing timeslices in a map; if a DECODED_BAD timeslice turns
	// into UNKNOWN, then we can just place the old data back afterwards.

	std::map<size_t, timeslice> best_timeslice_by_begin;

	decoder IBM_decoder;
	decoded_tracks decoded;

	size_t sectors = 0;

	for (double alpha = min_alpha; alpha <= max_alpha; alpha += stepsize) {
		std::cout << "Dewarping: alpha = " << alpha << std::endl;
		bool did_recode = false;
		for (timeslice & ts: floppy_line.timeslices) {
			if (ts.status != TS_DECODED_BAD) { continue; }

			// Store our backup.
			best_timeslice_by_begin[ts.flux_data_begin] = ts;

			double error;
			if (verbose) {
				std::cout << "Trying dewarp: indices: flux: " <<
					ts.flux_data_begin << " - " << ts.flux_data_end()
					<< ", sd: " << ts.sector_data_begin
					<< " - " << ts.sector_data_end() << ", alpha = "
					<< alpha;
				}

			// I'd really like get_mfm_train to automatically
			// handle the translation between levels on demand...
			// something more lazy perhaps.
			ts.mfm_train = get_MFM_train_dewarp(ts.clock_value,
					ts.flux_data, alpha, error);
			ts.sec_data = decode_MFM_train(ts.mfm_train,
				ts.preamble_offset, ts.mfm_train.data.size());

			if (verbose) {
				std::cout << " -- error: " << error << std::endl;
			}

			did_recode = true;
			ts.status = TS_UNKNOWN;
		}

		if (!did_recode) {
			continue;
		}

		// TODO: Don't be so chatty both times about what we've
		// decoded. Probably needs a decoder redesign.

		IBM_decoder.decode(floppy_line, decoded);

		// Restore any bad decodes that got turned into unknowns.
		// This is very iffy and needs a proper redesign; for that
		// matter, the timeline concept as such does. XXX

		for (timeslice & ts: floppy_line.timeslices) {
			if (ts.status != TS_UNKNOWN) {
				continue;
			}

			// Check if this timeslice exists in our backup
			// records, in case the preamble got corrupted by
			// an earlier decode.
			if (best_timeslice_by_begin.find(ts.flux_data_begin) !=
				best_timeslice_by_begin.end()) {

				// this must necessarily have status TS_DECODED_BAD;
				// otherwise it wouldn't have been added to the backup
				// record. So restore it.

				ts = best_timeslice_by_begin.find(
					ts.flux_data_begin)->second;
			}
		}
		sectors = decoded.last_decoded_sector;
		std::cout << sectors << "\n";
	}

	// Do a final decode so we've got something to return.
	IBM_decoder.decode(floppy_line, decoded);

	return decoded;
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
		preamble_info.ordinal_A1_sequence.needle, PREAMBLE_ID_A1);
	ordinal_preamble_search.add(
		preamble_info.ordinal_C2_sequence.needle, PREAMBLE_ID_C2);

	rabin_karp preamble_search(preamble_info.A1_sequence,
		PREAMBLE_ID_A1);
	preamble_search.add(preamble_info.C2_sequence,
		PREAMBLE_ID_C2);

	decoded_tracks all_decoded_tracks;
	size_t last_track = 0;

	if (flux_records.empty()) {
		std::cout << "Error: flux database contains no data!" << std::endl;
		return -1;
	}

	// TODO: Handle multiple tracks, perhaps by implementing a plus
	// operator for decoded_tracks.

	for (const flux_record & f: flux_records) {

		decoded_tracks decoded;
		std::vector<int> fluxes = f.fluxes;
		decoder IBM_decoder;

		std::cout << "Checking track " << f.track << ", side " << f.side << std::endl;
		last_track = std::max(last_track, (size_t)f.track);

		// Get locations/offsets into the flux where a magic preamble (A1A1A1
		// or C2C2C2) might be found.

		// Using a limited search strategy would cut the runtime to 62% for
		// uncorrupted floppies, but I don't know if it's worth the complexity,
		// so I won't do it (yet).

		std::vector<char> ordinal_flux = get_delta_coding(fluxes, true);
		std::vector<search_result> ordinal_locations =
			ordinal_preamble_search.find_matches(ordinal_flux);
		std::vector<match_with_clock> matches =
			get_flux_matches(fluxes, ordinal_locations,
				preamble_info);

		timeline floppy_line;
		// A linear sequence made up of each decoded chunk concatenated
		// in order.
		timeslice next;

		// We need to include a first timeslice covering everything
		// from the start of the flux record to the first match, even
		// if we don't know what's in there. This is necessary to make
		// flux indices in the timeline line up with those in the flux
		// record we're constructing it from.
		timeslice first;

		if (matches.size() > 0 && matches[0].match_location > 0) {
			first.flux_data = std::vector<int>(f.fluxes.begin(),
				f.fluxes.begin() + matches[0].match_location);
			floppy_line.insert(first);
		}

		for (size_t j = 0; j < matches.size(); ++j) {
			// We want to decode everything from this preamble to the
			// next one.
			match_with_clock m = matches[j];

			size_t start_idx = matches[j].match_location,
				end_idx = fluxes.size();

			if (j < matches.size()-1) {
				end_idx = matches[j+1].match_location;
			}

			// HACK: If the matched area is too short for a preamble, then
			// it's a false positive. Signal as such; it's mostly irrelevant
			// now due to non-overlapping search terms. However, it could
			// still theoretically happen.
			// I need some way to thread this through the filter_matches...
			// that it should be able to say "no, this must be truncated",
			// so I don't have to rely on every preamble being the same length
			// as I'm doing here.
			if (end_idx - start_idx < preamble_info.get_preamble_by_ID(0).size()) {
				std::cout << "Too short!\n";
				next.status = TS_TRUNCATED;
				floppy_line.insert(next);
				continue;
			}
			std::cout << "Found " << m.match_location << " with clock " <<
				m.estimated_clock << " (interval " << start_idx << "-" <<
				end_idx << ")";

			double error;

			next.clock_value = m.estimated_clock;
			next.flux_data = std::vector<int>(f.fluxes.begin() + start_idx,
				f.fluxes.begin() + end_idx);
			next.mfm_train = get_MFM_train(next.clock_value,
					next.flux_data, error);

			 std::cout << " -- Error: " << error << std::endl;

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
			next.preamble_offset = m.offset;

			next.sec_data = decode_MFM_train(next.mfm_train,
				next.preamble_offset, next.mfm_train.data.size());

			floppy_line.insert(next);

			if (floppy_line.timeslices.rbegin()->flux_data_begin != start_idx) {
				throw std::logic_error("Flux index start index mismatch. "
					"Was " + 
					itos(floppy_line.timeslices.rbegin()->flux_data_begin) +
					" but should be " + itos(start_idx));
			}
		}

		IBM_decoder.decode(floppy_line, decoded);

		// Do a brute force attempt to decode bad chunks that we haven't
		// got yet. This is very slow and could be sped up considerably with
		// some redesign, but the point is that we get a surprising amount
		// of corrupted chunks decoded "for free" this way.

		// Usually high alpha (less smoothing) is more promising, therefore
		// we only do a cursory check of the lower alpha parameters.
		// Known bug: checking alphas in reverse (from high to low) doesn't work.
		// I have no idea why.
		double stepsize = 0.005;
		decoded = decode_brute_dewarp(floppy_line, stepsize, 0.5, 0.1, false);
		decoded = decode_brute_dewarp(floppy_line, 0.5, 1-stepsize,
			stepsize, false);

		IBM_decoder.dump_sector_files(floppy_line,
			"check_sectors/t" + itos(f.track) + "h" + itos(f.side));

		all_decoded_tracks += decoded;
	}

	// Very quick and dirty hard-coded values, fix later. TODO
	decoder IBM_decoder_final;
	IBM_decoder_final.dump_image(all_decoded_tracks, "output",
		std::max((size_t)80, last_track+1), 2, 18, 512);

	return 0;
}
