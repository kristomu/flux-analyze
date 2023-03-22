#pragma once

#include <stdexcept>
#include <vector>
#include <cmath>

// For information about how these algorithms work, see
// the .cc file.

// Level one:

typedef std::vector<char> MFM_data_t;

class MFM_train_data {
	public:
		// This is the MFM train data itself.
		MFM_data_t data;

		// This is what flux delay a given bit of the MFM train
		// belongs to. It's used to keep track of what pieces of
		// flux data we've covered.
		std::vector<size_t> flux_indices;

		// Append some other data encoded by another process.
		// This may be turned into something more sophisticated
		// later, but my indecision shows that at this point I just
		// need to iterate and then fix.
		MFM_train_data & operator+= (const MFM_train_data & after) {
			if (!flux_indices.empty() && (
				after.flux_indices[0] < *flux_indices.rbegin() ||
				after.flux_indices[0] - *flux_indices.rbegin() > 1)) {

				throw std::logic_error("Can't append nonconsecutive data!");
			}
			// If we end on a 1 and the other data start on a 1, chop
			// off the last index. This happens when concatenating
			// MFM trains that have been produced from flux data.
			if (!data.empty() && after.data[0] == 1 && *data.rbegin() == 1) {
				data.pop_back();
				flux_indices.pop_back();
			}

			std::copy(after.data.begin(), after.data.end(),
				std::back_inserter(data));
			std::copy(after.flux_indices.begin(), after.flux_indices.end(),
				std::back_inserter(flux_indices));

			return *this;
		}
};

// Turn a flux record vector into an MFM pulse train by comparing the distance
// between flux reversals. The error_out double is set to a badness-of-fit
// value: the higher the worse.

MFM_train_data get_MFM_train(double clock,
		const std::vector<int> & fluxes, size_t start_pos,
		size_t end_pos, double & error_out);

MFM_train_data get_MFM_train(double clock,
		const std::vector<int> & fluxes, double & error_out);

double get_MFM_train_error(double clock, const std::vector<int> & fluxes);

MFM_train_data get_MFM_train_dewarp(double clock,
	const std::vector<int> & fluxes, double & error_out);
