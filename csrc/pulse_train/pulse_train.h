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

		// Split a subsection of the MFM train data into its own MFM train data.
		// It will split on whole flux transitions because that's what decode_MFM_
		// train and timelines expect (there's a lot of ugly code here), and
		// return the offset to the actual start that was specified.
		MFM_train_data split(size_t start_idx_train, size_t end_idx_train,
			size_t & out_offset) {
			MFM_train_data out;

			out_offset = 0;

			size_t start_idx_flux = flux_indices[start_idx_train];

			// Is this ever used???

			/*while(start_idx_train != 0 &&
				flux_indices[start_idx_train-1] == start_idx_flux) {
				++out_offset;
				--start_idx_train;
			}*/

			out.data = MFM_data_t(data.begin() + start_idx_train,
				data.begin() + end_idx_train);

			for (size_t i = start_idx_train; i < end_idx_train; ++i) {
				out.flux_indices.push_back(flux_indices[i]-start_idx_flux);
			}

			return out;
		}
};