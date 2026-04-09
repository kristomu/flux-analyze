#include <iostream>
#include <stdexcept>

#include "preambles.h"
#include "stepwise_ewma.h"

#include "ewma_search/ewma_search.h"

MFM_train_data stepwise_EWMA_decoder::get_MFM_train(
	const std::vector<int> & fluxes, size_t start_pos,
	size_t end_pos, double & error_out) const {

	// Get the A1A1A1 preamble (fix this later)

	/* std::vector<char> short_A1 = {
       	0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1};
		needs to be repeated thrice. In addition, the first
		0 1 should be discarded because we don't know if
		that's a 01, 001, or 0001. Anything will match it.

		So we get the following nonzero sequence:
			3 2 3 2 1
			3 2 3 2 1
			3 2 3 2
		hence the following half-clocks:

	*/

	std::vector<int> half_clock_needle = {
		4, 3, 4, 3, 2,
		4, 3, 4, 3, 2,
		4, 3, 4, 3 };

	error_out = std::numeric_limits<double>::infinity();

	std::vector<int> one_way_truncated_fluxes(
		fluxes.begin() + start_pos,
		fluxes.begin() + end_pos);

	// Look for EWMA matches.

	std::vector<std::pair<int, interval> > matches =
		ewma_search_all(one_way_truncated_fluxes,
			half_clock_needle, alpha, tolerance);

	std::cout << matches.size() << " matches.\n";

	throw std::runtime_error("Missing functionality"); // fix later
}