#pragma once

#include "flux_decoder.h"
#include "stats/kmedian/kmedian.h"

class kmedian_decoder: public pulse_decoder {
	
	public:
		// Do I want to do this everywhere, or just do a lazy
		// two-step function approach instead? Hmm...
		using pulse_decoder::get_MFM_train;

		MFM_train_data get_MFM_train(
			const std::vector<int> & fluxes, size_t start_pos,
			size_t end_pos, double & error_out) const;

		std::string name() const {
			return "K-median constant clock decoder";
		}
};