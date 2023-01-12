#pragma once
#include <vector>

// For information about how these algorithms work, see
// the .cc file.

// Level one:

// Turn a flux record vector into an MFM pulse train by comparing the distance
// between flux reversals. The error_out double is set to a badness-of-fit
// value: the higher the worse.

std::vector<char> get_MFM_train(double clock,
		const std::vector<int> & fluxes, double & error_out);

double get_MFM_train_error(double clock, const std::vector<int> & fluxes);

// Baseline inference for clocks. This should be very quick and work for
// most non-corrupted floppies.

double infer_clock(const std::vector<int> & fluxes);