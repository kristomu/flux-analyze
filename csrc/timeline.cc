#pragma once
#include "pulse_train.h"
#include "sector_data.h"

class timeslice {
	public:
		// TODO, disambiguate from flux_record or refactor it.
		std::vector<int> flux_data;
		MFM_train_data mfm_train;
		sector_data sector_data;

		// These offsets give the relative position of the beginning
		// of each chunk, relative to the first chunk in the timeline.
		size_t flux_data_begin;
		size_t mfm_train_begin;
		size_t sector_data_begin;
};

class timeline {
	public:
		std::vector<timeslice> timeslices;
};