#include <algorithm>

#include "ewma.h"

std::vector<double> causal_get_approximate_clock(
	const std::vector<int> & fluxes, double alpha,
	double clock_prior, bool use_prior) {

	std::vector<double> appx_clock;

	double mean_clock = 1; // or insert your prior here
	size_t warmup_period = 300;

	if (use_prior) {
		mean_clock = clock_prior;
		warmup_period = 0;
	}

	for (size_t i = 0; i < fluxes.size(); ++i) {
		size_t observed_pulse_delay = fluxes[i];
		
		size_t zero_bits = round(
			2.0 * observed_pulse_delay / mean_clock) - 1;
		zero_bits = std::min(3, std::max(1, (int)zero_bits));

		double new_clock = 2 * observed_pulse_delay / (zero_bits + 1);

		mean_clock = (1-alpha) * mean_clock + alpha * new_clock;

		if (i < warmup_period) {
			appx_clock.push_back(-1);
		} else {
			appx_clock.push_back(mean_clock);
		}
	}

	return appx_clock;
}

std::vector<double> causal_get_approximate_clock(
	const std::vector<int> & fluxes, double alpha) {

	return causal_get_approximate_clock(fluxes, alpha,
		1, false);
}


std::vector<double> acausal_get_approximate_clock(
	std::vector<int> fluxes, double alpha) {

	// Forward pass
	std::vector<double> fwd = causal_get_approximate_clock(
		fluxes, 1-(1-alpha)/2.0);

	// Backward pass
	std::reverse(fluxes.begin(), fluxes.end());

	double prior = *fwd.rbegin();
	std::vector<double> backward = causal_get_approximate_clock(
		fluxes, 1-(1-alpha)/2.0, prior, true);

	std::reverse(backward.begin(), backward.end());

	// Average the two passes to cancel out lag (except at the beginning and
	// end). Strictly speaking, the -1 regions should be done with the actual
	// alpha (not halved towards one) with a starting value taken from the
	// first point where both forward and backward datapoints exist. But
	// that's overkill.

	std::vector<double> appx_clock(fluxes.size());

	for (size_t i = 0; i < appx_clock.size(); ++i) {
		if (fwd[i] == -1) {
			appx_clock[i] = backward[i];
			continue;
		}

		if (backward[i] == -1) {
			appx_clock[i] = fwd[i];
			continue;
		}

		appx_clock[i] = (backward[i] + fwd[i]) / 2.0; 
	}

	return appx_clock;
}

MFM_train_data decode_by_clock(const std::vector<int> & fluxes,
	const std::vector<double> & clock_values,
	int starting_offset) {

	MFM_train_data train;

	// Beware: this doesn't teka parity into account.
	// It probably works most of the time here because
	// the outer method starts decoding at an OOB marker.

	train.data.push_back(1);
	train.flux_indices.push_back(starting_offset);

	std::vector<size_t> pulse_train;

	for (size_t i = 0; i < fluxes.size(); ++i) {
		int observed_pulse_delay = fluxes[i];
                
		size_t zero_bits = round(
			2.0 * observed_pulse_delay / clock_values[i]) - 1;

		zero_bits = std::min(3, std::max(1, (int)zero_bits));

		for (size_t j = 0; j < zero_bits; ++j) {
			train.data.push_back(0);
			train.flux_indices.push_back(i + starting_offset);
		}
	}

	return train;
}

// OOP glue below.

MFM_train_data causal_EWMA_clock_decoder::get_MFM_train(
	const std::vector<int> & fluxes, size_t start_pos,
	size_t end_pos, double & error_out) const {

	// Error out is not supported, as optimizing
	// based on it obviously overfits. I may entirely
	// remove it later, or make it report the actual
	// objective value (according to a smoothing factor
	// of alpha).
	error_out = std::numeric_limits<double>::infinity();

	// It would be better to make the lower functions handle
	// start and end positions. Do that later.
	std::vector<int> cropped(fluxes.begin()+start_pos,
			fluxes.begin()+end_pos);

	std::vector<double> clock_values = causal_get_approximate_clock(
		cropped, alpha);

	// Use the last value as a prior to fill in the first values.
	// (This is kind of dirty.)

	clock_values = causal_get_approximate_clock(cropped, alpha,
		*clock_values.rbegin(), true);
	
	return decode_by_clock(cropped, clock_values, start_pos);
}

MFM_train_data acausal_EWMA_clock_decoder::get_MFM_train(
	const std::vector<int> & fluxes, size_t start_pos,
	size_t end_pos, double & error_out) const {

	error_out = std::numeric_limits<double>::infinity();

	std::vector<int> cropped(fluxes.begin()+start_pos,
			fluxes.begin()+end_pos);

	std::vector<double> clock_values = acausal_get_approximate_clock(
		cropped, alpha);
	
	return decode_by_clock(cropped, clock_values, start_pos);
}