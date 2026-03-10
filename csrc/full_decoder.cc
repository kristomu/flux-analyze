#include "full_decoder.h"

#include "pulse_train/all.h"

// ORDINAL
#include "ordinal_search.h"

// ONCE-THROUGH
// (to add includes here)

decoded_tracks full_decoder::decode_track(const flux_record & f) const {
	std::cout << "Checking track " << f.track << ", side " << f.side << std::endl;

	return decode_track(f.fluxes);
}

decoded_tracks full_decoder::decode_floppy(
	const std::vector<flux_record> & flux_records) const {

	decoded_tracks all_decoded_tracks;

	if (flux_records.empty()) {
		std::cout << "Error: flux database contains no data!" << std::endl;
		return all_decoded_tracks;
	}

	for (const flux_record & f: flux_records) {
		all_decoded_tracks += decode_track(f);
	}

	return all_decoded_tracks;
}

decoded_tracks full_decoder::decode_floppy(std::string flux_filename) const {
	std::vector<flux_record> flux_records = 
		get_flux_record(flux_filename, true);

	return decode_floppy(flux_records);
}


// ordinal below

decoded_tracks ordinal_full_decoder::decode_track(
	const std::vector<int> & fluxes) const {

	rabin_karp ordinal_preamble_search(
		preamble_info.ordinal_A1_sequence.needle, PREAMBLE_ID_A1);
	ordinal_preamble_search.add(
		preamble_info.ordinal_C2_sequence.needle, PREAMBLE_ID_C2);

	decoded_tracks decoded;
	decoder IBM_decoder;

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
			preamble_info, verbose);

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
		first.flux_data = std::vector<int>(fluxes.begin(),
			fluxes.begin() + matches[0].match_location);
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

		if (verbose) {
			std::cout << "[" << j << "/" << matches.size() << "]\tStart idx: " << matches[j].match_location << ", end idx: "
				<< end_idx << " offset " << matches[j].offset << "\n";
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
			if (verbose) {
				std::cout << "Too short!\n";
			}
			next.status = TS_TRUNCATED;
			floppy_line.insert(next);
			continue;
		}

		if (verbose) {
			std::cout << "Found " << m.match_location << " with clock " <<
				m.estimated_clock << " (interval " << start_idx << "-" <<
				end_idx << ")";
			}

		double error;

		// TODO: Allow any kind of clock decoder as an input -- or
		// even something like a "brute-force" band-based decoder
		// like in Python.

		constant_clock_decoder c_MFM_decoder;
		causal_EWMA_clock_decoder variable_decoder; // e.g.

		next.clock_value = m.estimated_clock;
		next.flux_data = std::vector<int>(fluxes.begin() + start_idx,
			fluxes.begin() + end_idx);
		c_MFM_decoder.set_clock(next.clock_value);
		variable_decoder.set_initial_clock(next.clock_value);
		variable_decoder.set_alpha(0.001);
		next.mfm_train = c_MFM_decoder.get_MFM_train(
			next.flux_data, error);

		if (verbose) {
			std::cout << " -- Error: " << error << std::endl;
		}

		// HACK HACK: Find the offset from the start of the MFM train
		// to the start of the preamble. (TODO: lots of optimization
		// can be done here; e.g. trying variable decoder and then
		// falling back onto c_MFM_decoder on failure. Not yet though.)
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

	bool show_stats = true;//verbose;
	IBM_decoder.decode(floppy_line, decoded, verbose, show_stats);

	// TODO: Get stats in a more unified manner than this...
	// I do have a stats structure now; use it for that, perhaps.

	if (verbose) {
		std::vector<size_t> covered = floppy_line.get_covered_flux_sizes();
		std::cout << "DEBUG: Timeline occupancy: ";
		std::copy(covered.begin(), covered.end(),
			std::ostream_iterator<size_t>(std::cout, " "));
		std::cout << "\n";
		std::cout << "Fraction of flux transitions decoded as good: " << floppy_line.get_good_fraction() << "\n";
	}

	return decoded;
}

// once-through:

decoded_tracks once_through_decoder::decode_track(
	const std::vector<int> & fluxes) const {

	double error = 0;

	MFM_train_data full_decoding =
		flux_to_MFM_decoder->get_MFM_train(fluxes, error);

	std::cout << "[Simpler decode] Error: " << error << "\n";

	std::vector<search_result> preamble_positions =
		preamble_search.find_matches(full_decoding.data);

	if (verbose) {
		for (const search_result & sr: preamble_positions) {
			std::cout << "Found preamble at " << sr.idx << "\t";
			switch(sr.ID) {
				case PREAMBLE_ID_A1: std::cout << "A1A1A1\n"; break;
				case PREAMBLE_ID_C2: std::cout << "C2C2C2\n"; break;
				default: std::cout << "???\n"; break;
			}
		}
	}

	std::cout << "Total length: " << full_decoding.data.size() << "\n";
	std::cout << std::endl;

	// Piece together a timeline and do a decoding from pulse train
	// to obtain more structured sector data.

	timeline floppy_line;

	decoder IBM_decoder;
	decoded_tracks decoded;

	if (preamble_positions.size() == 0) {
		return decoded; // nothing to do
	}

	// We need to include a first timeslice covering everything
	// from the start of the flux record to the first match, even
	// if we don't know what's in there. This is necessary to make
	// flux indices in the timeline line up with those in the flux
	// record we're constructing it from.
	timeslice first;

	// We'll be using two kinds of indices: indices into the flux
	// transition data and into the MFM pulse train. These must not
	// be confused, so I'll label them start/end_idx_flux and
	// start/end_idx_train. Search result indices are train
	// indices.

	if (preamble_positions[0].idx > 0) {
		size_t start_idx_train = preamble_positions[0].idx;
		size_t start_idx_flux = full_decoding.flux_indices[start_idx_train];

		first.flux_data = std::vector<int>(fluxes.begin(),
			fluxes.begin() + start_idx_flux);
		size_t offset = 0;

		first.mfm_train = full_decoding.split(0, start_idx_train, offset);

		floppy_line.insert(first);
	}

	// A linear sequence made up of each decoded chunk concatenated
	// in order.
	timeslice next;

	// This should be a stat of some sort...
	std::cout << preamble_positions.size() << " address mark markers found\n";

	for (size_t i = 0; i < preamble_positions.size(); ++i) {
		if (verbose) {
			std::cout << i << " of " << preamble_positions.size() << "\n";
		}
		size_t start_idx_train = preamble_positions[i].idx,
			end_idx_train = full_decoding.data.size();

		size_t start_idx_flux = full_decoding.flux_indices[start_idx_train],
			end_idx_flux;

		if (i+1 < preamble_positions.size()) {
			end_idx_train = preamble_positions[i+1].idx;
			end_idx_flux = full_decoding.flux_indices[end_idx_train];
		} else {
			end_idx_flux = fluxes.size();
		}

		// Skip TS_TRUNCATED check - should be OK as long as our needles
		// are non-overlapping (which they are, since IBM_preamble
		// includes the first all high bit of the F nibble that every
		// address mark contains).

		size_t offset = 0;

		while(full_decoding.flux_indices[start_idx_train-1] == start_idx_flux) {
			++offset;
			--start_idx_train;
		}

		if (verbose) {
			std::cout << "(Flux) Start idx: " << start_idx_flux<< ", end idx: "
				<< end_idx_flux << " offset " << offset << "\n";
			std::cout << "(Train) Start idx: " << start_idx_train << ", end idx: "
				<< end_idx_train << " offset " << offset << "\n";
		}

		// This seems to work, but I feel kinda iffy about it because
		// all the old code - particularly the timeslice stuff - is so
		// uneven.

		// Also whether offset is even used is kinda unclear. I need
		// a better timeline/timeslice structure.

		next.status = TS_PREAMBLE_FOUND;
		next.preamble_offset = offset; // TBD - this will be a mess...

		next.flux_data = std::vector<int>(fluxes.begin() +
			start_idx_flux, fluxes.begin() + end_idx_flux);
		next.mfm_train = full_decoding.split(start_idx_train, end_idx_train,
			offset);

		next.sec_data = decode_MFM_train(next.mfm_train,
			next.preamble_offset, next.mfm_train.data.size());

		floppy_line.insert(next);
	}

	bool show_stats = verbose;

	IBM_decoder.decode(floppy_line, decoded, verbose, show_stats);

	std::vector<size_t> covered = floppy_line.get_covered_flux_sizes();
	if (verbose) {
		std::cout << "DEBUG: Timeline occupancy: ";
		std::copy(covered.begin(), covered.end(),
			std::ostream_iterator<size_t>(std::cout, " "));
		std::cout << "\n";
	}
	std::cout << "Fraction of flux transitions decoded as good: " << floppy_line.get_good_fraction() << "\n";

	//throw std::logic_error("Job's done!");

	return decoded;
}