
// Grid search based on the standalone decoder.
// Enabling the DP-based MFM decoder is very slow and
// thus has been commented out by default.

#include <algorithm>
#include <iostream>

#include <fstream>
#include <string>
#include <stdexcept>
#include <vector>
#include <cmath>

#include <zlib.h>

#include "flux_record.h"

#include "pulse_train/all.h"
#include "pulse_train/dp/all.h"

#include "ordinal_search.h"
#include "sector_data.h"
#include "tools.h"

#include "timeline.h"
#include "decoder.h"
#include "full_decoder.h"

#include "crc16.h"

void do_grid_search(
	const std::vector<flux_record> & flux_records,
	const once_through_decoder & full_decoder,
	std::shared_ptr<pulse_decoder> & pulse_decoder,
	const std::vector<param_t> & param_type,
	const std::vector<double> & minimum,
	const std::vector<double> & maximum,
	const std::vector<int> & stepsize,
	std::vector<double> & current_params,
	size_t pos) {

	if (pos == current_params.size()) {
		// We've set all parameters; do a test.
		pulse_decoder->set_params(current_params);
		decoded_tracks all_decoded_images = full_decoder.decode_floppy(flux_records);

		bool revealed = false;

		// Output the performance figures.

		for (const std::pair<int, decoder_stats> stats_pairs: all_decoded_images.stats_per_track) {
			std::cout << pulse_decoder->name() << ":\t"
				<< "Track " << stats_pairs.first << ": "
				<< stats_pairs.second.num_recovered_sectors << " of "
				<< stats_pairs.second.num_sectors << " sectors recovered.\n";

			// Semicolon-separated values.
			// First the outcome (so we can sort regardless of parameter space)...
			std::cout << pulse_decoder->name() << ";";
			std::cout << stats_pairs.second.num_recovered_sectors
				<< ";" << stats_pairs.second.num_sectors
				<< ";" << stats_pairs.second.flux_fraction_good << ";";
			// Then the parameters.
			std::copy(current_params.begin(), current_params.end(),
				std::ostream_iterator<double>(std::cout, ";"));
			std::cout << "\n";

			revealed = true;
		}

		if (!revealed) {
			std::cout << pulse_decoder->name() << ";";
			std::cout << 0 << ";" << 0
				<< ";" << "N/A" << ";";
			std::copy(current_params.begin(), current_params.end(),
				std::ostream_iterator<double>(std::cout, ";"));
			std::cout << "\n";
		}

		return; // TODO: return some objective value
	}

	// Recurse on the given position.
	double last_step_val = minimum[pos]-1, current_step_val = minimum[pos]-1;

	// TODO: Fix off-by-one
	for (int current_step = 0; current_step < stepsize[pos]; ++current_step) {
		switch(param_type[pos]) {
			case PARAM_REAL:
				current_step_val = minimum[pos] +
					current_step/(double)(stepsize[pos]-1) * (maximum[pos]-minimum[pos]);
				break;
			case PARAM_INTEGER:
				current_step_val = minimum[pos] +
					current_step/(double)(stepsize[pos]-1) * (maximum[pos]-minimum[pos]);
				current_step_val = round(current_step_val);
				break;
			case PARAM_LOG_REAL:
				//current_step_val = (exp(current_step/(double)stepsize[pos])-exp(0))/(exp(1) - exp(0));
				current_step_val = exp(current_step-(stepsize[pos]-1));
				current_step_val = minimum[pos] +
					current_step_val * (maximum[pos]-minimum[pos]);
		}

		// Have we already evaluated this option? (E.g. for integer or booleans)
		if (current_step_val == last_step_val) {
			continue; // If so, don't do it again; there's no point.
		}

		current_params[pos] = current_step_val;

		do_grid_search(flux_records, full_decoder, pulse_decoder, param_type,
			minimum, maximum, stepsize, current_params, pos+1);

		last_step_val = current_step_val;
	}
}

void do_grid_search(
	const std::vector<flux_record> & flux_records,
	std::shared_ptr<pulse_decoder> pulse_decoder) {

	once_through_decoder full_decoder(pulse_decoder);

	size_t num_params = pulse_decoder->get_parameter_types().size();

	std::vector<double> current_params(num_params, 0);

	int desired_total_numiters = 1600;
	int stepsize_guess = 2;

	std::vector<int> stepsizes(num_params);
	bool found_match = false;

	auto param_types = pulse_decoder->get_parameter_types();
	auto minimum = pulse_decoder->get_parameter_min();
	auto maximum = pulse_decoder->get_parameter_max();

	std::cout << "adapting stepsizes to get a total of " << desired_total_numiters << "\n";

	// Find the right stepsize, Lagrange multiplier style.
	for (stepsize_guess = 2; stepsize_guess < 100 && !found_match; ++stepsize_guess) {
		int total_numiters = 1;
		for (size_t i = 0; i < param_types.size(); ++i) {
			stepsizes[i] = stepsize_guess;

			if (param_types[i] != PARAM_INTEGER) {
				// If the param type is real, we can subdivide arbitrarily fine, so
				// subdivide to stepsize_guess.
				total_numiters *= stepsize_guess;
			} else {
				// Otherwise, there's a finite number of choices, and going higher
				// won't cost us more.
				total_numiters *= std::min(maximum[i] - minimum[i] + 1,
					(double)stepsize_guess);
			}
		}
		found_match = total_numiters >= desired_total_numiters;
		std::cout << "stepsize_guess = " << stepsize_guess
			<< " gave " << total_numiters << "\n";
		++stepsize_guess;
	}

	do_grid_search(flux_records,
		full_decoder, pulse_decoder,
		pulse_decoder->get_parameter_types(),
		pulse_decoder->get_parameter_min(),
		pulse_decoder->get_parameter_max(),
		stepsizes, current_params, 0);
}

int main(int argc, char ** argv) {
	test_rabin_karp();
	test_crc16();

	ordinal_full_decoder ofd;

	// Example use of the once-through full decoder.
	std::shared_ptr<alternating_EWMA_decoder> AEC =
		std::make_shared<alternating_EWMA_decoder>();
	std::shared_ptr<causal_EWMA_clock_decoder> cacd =
		std::make_shared<causal_EWMA_clock_decoder>();
	std::shared_ptr<causal_EWMA_clock_decoder> acacd =
		std::make_shared<causal_EWMA_clock_decoder>();
	std::shared_ptr<orig_EWMA_causal_clock_decoder> orig_causal_clock =
		std::make_shared<orig_EWMA_causal_clock_decoder>();
	std::shared_ptr<constant_clock_decoder> constant_clock =
		std::make_shared<constant_clock_decoder>();
	std::shared_ptr<historical_EWMA_decoder> historical =
		std::make_shared<historical_EWMA_decoder>();
	cacd->set_initial_clock(24);
	acacd->set_initial_clock(24);
	//p_decoder->set_alpha(0.01);

	std::shared_ptr<QND_one_mode_DP> omd =
		std::make_shared<QND_one_mode_DP>();

	std::vector<std::shared_ptr<pulse_decoder> > pulse_decoders;
	//pulse_decoders.push_back(omd);
	pulse_decoders.push_back(AEC);
	pulse_decoders.push_back(cacd);
	pulse_decoders.push_back(acacd);
	pulse_decoders.push_back(constant_clock);
	pulse_decoders.push_back(orig_causal_clock);
	pulse_decoders.push_back(historical);
	std::reverse(pulse_decoders.begin(), pulse_decoders.end());

	// TODO: Grid search only tests one sector. We should only accept one
	// image, or splice together the sector images, or output grid search
	// results for all of them (though the latter will be very expensive
	// and the user would really be better off just launching a separate
	// process for each search).
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " <flux image> ... <flux image>\n";
		std::cout << "Specifying multiple images will decode each as a different\n";
		std::cout << "image of the same floppy; any sector that's valid in at\n";
		std::cout << "least one of them will be preserved in the output.\n";
		return -1;
	}

	std::string flux_filename = argv[1];
	std::cout << "Analyzing flux file " << flux_filename << "\n";

	std::vector<flux_record> flux_records =
		get_flux_record(flux_filename, true);

	for (std::shared_ptr<pulse_decoder> p_decoder: pulse_decoders) {
		do_grid_search(flux_records, p_decoder);
	}

	return 0;
}
