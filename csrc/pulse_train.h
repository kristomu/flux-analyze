#pragma once

#include <stdexcept>
#include <vector>
#include <cmath>

// For information about how these algorithms work, see
// the .cc file.

typedef std::vector<char> MFM_data_t;

class MFM_train_data {
	public:
		// This is the MFM train data itself.
		MFM_data_t data;

		// This is what flux delay a given bit of the MFM train
		// belongs to. It's used to keep track of what pieces of
		// flux data we've covered.
		std::vector<size_t> indices;

		// Append some other data encoded by another process.
		// This may be turned into something more sophisticated
		// later, but my indecision shows that at this point I just
		// need to iterate and then fix.
		MFM_train_data & operator+= (const MFM_train_data & after) {
			if (!indices.empty() && (
				after.indices[0] - *indices.rbegin() < 0 ||
				after.indices[0] - *indices.rbegin() > 1)) {

				throw std::logic_error("Can't append nonconsecutive data!");
			}
			std::copy(after.data.begin(), after.data.end(),
				std::back_inserter(data));
			std::copy(after.indices.begin(), after.indices.end(),
				std::back_inserter(indices));

			return *this;
		}
};


// Level one:

// Turn a flux record vector into an MFM pulse train by comparing the distance
// between flux reversals. The error_out double is set to a badness-of-fit
// value: the higher the worse.

MFM_train_data get_MFM_train(double clock,
		const std::vector<int> & fluxes, size_t start_pos,
		size_t end_pos, double & error_out);

MFM_train_data get_MFM_train(double clock,
		const std::vector<int> & fluxes, double & error_out);

double get_MFM_train_error(double clock, const std::vector<int> & fluxes);

// Baseline inference for clocks. This should be very quick and work for
// most non-corrupted floppies.

double infer_clock(std::vector<int>::const_iterator fluxes_start,
	std::vector<int>::const_iterator fluxes_end);

double infer_clock(const std::vector<int> & fluxes);