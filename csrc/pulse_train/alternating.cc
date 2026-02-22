
#include "alternating.h"

std::vector<double> alternating_EWMA_decoder::get_new_clocks_filtfilt(
	const std::vector<int> & fluxes,
	const std::vector<double> & steps) const {

	// TODO: Explain this based on an approximation of
	// the EWMA frequency response with low alpha, as
	// H(w) = 1/(1 + lambda w)
	double alpha_mark = (2 + sqrt(2)) * alpha / (alpha + sqrt(2) + 1);

	// Forward pass.

	double mean = fluxes[0]/steps[0];

	double instant_clock, current;

	std::vector<double> clocks(fluxes.size());

	for (size_t i = 0; i < clocks.size(); ++i) {
		instant_clock = fluxes[i]/steps[i];
		current = alpha_mark * instant_clock + (1-alpha_mark) * mean;
		mean = current;

		clocks[i] = mean;
	}

	// Backward pass.

	std::vector<double> back_clocks(clocks.size());

	for (size_t i = 0; i < clocks.size(); ++i) {
		size_t j = clocks.size() - (i + 1);

		instant_clock = clocks[j];

		mean = alpha_mark * instant_clock + (1-alpha_mark) * mean;
		back_clocks[j] = mean;
	}

	return back_clocks;
}

std::vector<double> alternating_EWMA_decoder::get_new_clocks_avg(
	const std::vector<int> & fluxes,
	const std::vector<double> & steps) const {

	// This should have roughly the same frequency response
	// for a given alpha, no tweaking required.
	// (I think? I'm not an EE.)

	// Forward pass.

	double mean = fluxes[0]/steps[0];

	double instant_clock, current;

	std::vector<double> clocks(fluxes.size());

	for (size_t i = 0; i < clocks.size(); ++i) {
		instant_clock = fluxes[i]/steps[i];
		current = alpha * instant_clock + (1-alpha) * mean;
		mean = current;

		clocks[i] = mean;
	}

	// Backward pass.

	std::vector<double> back_clocks(clocks.size());

	for (size_t i = 0; i < clocks.size(); ++i) {
		size_t j = clocks.size() - (i + 1);

		instant_clock = fluxes[i]/steps[i];

		mean = alpha * instant_clock + (1-alpha) * mean;
		back_clocks[j] = mean;
	}

	// Combine the two.

	for (size_t i = 0; i < clocks.size(); ++i) {
		clocks[i] = (clocks[i] + back_clocks[i]) / 2.0;
	}

	return clocks;
}

std::vector<double> alternating_EWMA_decoder::get_new_clocks(
	const std::vector<int> & fluxes,
	const std::vector<double> & steps) const {

	switch (filter_type) {
		case AFT_FILTFILT:
			return get_new_clocks_filtfilt(fluxes, steps);
		case AFT_AVERAGE_FILT:
			return get_new_clocks_avg(fluxes, steps);
		default:
			throw std::logic_error("get_new_clocks: Unknown filter type!");
	}
}

std::vector<double> alternating_EWMA_decoder::get_new_steps(
	const std::vector<int> & fluxes,
	const std::vector<double> & clocks,
	int iteration) const {

	std::vector<double> steps;

	for (size_t i = 0; i < fluxes.size(); ++i) {
		// This is the optimal value of the step
		// variable if the values it can take are
		// completely unconstrained.
		double linear_step = fluxes[i]/clocks[i];

		// But if we want a perfect match to a valid
		// MFM number of zero bits, we need something
		// more like this. (Simplify later)

		int zero_bits = round(2.0 * fluxes[i] / clocks[i]) - 1;
		double new_clock = 2 * fluxes[i] / (zero_bits + 1.0);
		double quantized_step = fluxes[i] / new_clock;

		// Then do some linear combination of the two. I'm just
		// going to be mostly linear here. Experiment with a
		// changing temperature or whatever later.

		double linear_mix = linear_mix_intercept - linear_mix_slope * iteration;
		linear_mix = std::min(1.0, std::max(0.0, linear_mix));

		steps.push_back(linear_mix * linear_step + (1-linear_mix) * quantized_step);
	}

	return steps;
}

// Cut and paste from ewma.cc. TODO, fix later

MFM_train_data decode_by_clock_copy(const std::vector<int> & fluxes,
	const std::vector<double> & clock_values,
	int starting_offset) {

	MFM_train_data train;

	// Beware: this doesn't take parity into account.
	// It probably works most of the time here because
	// the outer method starts decoding at an OOB marker.

	std::vector<size_t> pulse_train;

	train.data.push_back(1);
	train.flux_indices.push_back(starting_offset);

	for (size_t i = 0; i < fluxes.size(); ++i) {
		int observed_pulse_delay = fluxes[i];

		/*if (clock_values[i] == CLOCK_OOB_DO_SKIP) {
			continue;
		}*/
                
		int zero_bits = round(
			2.0 * observed_pulse_delay / clock_values[i]) - 1;

		zero_bits = std::max(1, zero_bits);

		if (zero_bits < 0) {
			throw std::invalid_argument("decode_by_clock: invalid clock stream "
				"(negative number of zero bits!)");
		}

		// Trust the input to otherwise have done sensible clamping, so
		// we don't need to do any clamping on our own end. Hence
		// there's no std::min(3...) here.

		for (int j = 0; j < zero_bits; ++j) {
			train.data.push_back(0);
			train.flux_indices.push_back(i + starting_offset);
		}

		train.data.push_back(1);
		train.flux_indices.push_back(i + starting_offset);
	}

	return train;
}

MFM_train_data alternating_EWMA_decoder::get_MFM_train(
	const std::vector<int> & fluxes, size_t start_pos,
	size_t end_pos, double & error_out) const {

	error_out = std::numeric_limits<double>::infinity();

	std::vector<int> cropped_fluxes(fluxes.begin() + start_pos,
		fluxes.begin() + end_pos);

	// Default prior, make something uninformative later.
	std::vector<double> clocks(cropped_fluxes.size(), 22);
	std::vector<double> steps(cropped_fluxes.size(), 1);

	//clocks = get_new_clocks(fluxes, steps);

	for (size_t iter = 0; iter < num_iterations; ++iter) {
		steps = get_new_steps(fluxes, clocks, iter);
		clocks = get_new_clocks(fluxes, steps);
	}

	return decode_by_clock_copy(fluxes, clocks, start_pos);
}