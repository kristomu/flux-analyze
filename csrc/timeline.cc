#include "pulse_train.h"
#include "sector_data.h"
#include "timeline.h"

void timeline::insert(timeslice & next) {
	if (!timeslices.empty()) {
		auto lastptr = timeslices.rbegin();
		next.flux_data_begin = lastptr->flux_data_end();
		next.mfm_train_begin = lastptr->mfm_train_end();
		next.sector_data_begin = lastptr->sector_data_end();
	}

	if (next.flux_data.empty()) {
		throw std::invalid_argument("timeslice must have flux data");
	}

	if (next.mfm_train.data.empty()) {
		throw std::invalid_argument("timeslice must have flux data");
	}

	timeslices.push_back(next);
}