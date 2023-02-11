#pragma once
#include <vector>
#include <cmath>

// For information about how these algorithms work, see
// the .cc file.

typedef std::vector<char> MFM_train_t;

// Level one:

// Turn a flux record vector into an MFM pulse train by comparing the distance
// between flux reversals. The error_out double is set to a badness-of-fit
// value: the higher the worse.

MFM_train_t get_MFM_train(double clock,
		const std::vector<int> & fluxes, size_t start_pos,
		size_t end_pos, double & error_out);

MFM_train_t get_MFM_train(double clock,
		const std::vector<int> & fluxes, double & error_out);

double get_MFM_train_error(double clock, const std::vector<int> & fluxes);

// Baseline inference for clocks. This should be very quick and work for
// most non-corrupted floppies.

double infer_clock(std::vector<int>::const_iterator fluxes_start,
	std::vector<int>::const_iterator fluxes_end);

double infer_clock(const std::vector<int> & fluxes);