#include "pulse_train.h"
#include "sector_data.h"
#include "timeline.h"

#include <assert.h>

void timeline::insert(timeslice & next, int ID) {
	// Enforce the uniqueness constraint: we can't have two
	// with the same ID.
	if (used_IDs.find(ID) != used_IDs.end()) {
		throw std::invalid_argument("timeline: Can't insert timeslice"
			" with non-unique ID");
	}

	if (!timeslices.empty()) {
		auto lastptr = timeslices.rbegin();
		next.flux_data_begin = lastptr->flux_data_end();
		next.mfm_train_begin = lastptr->mfm_train_end();
		next.sector_data_begin = lastptr->sector_data_end();
	}

	if (next.flux_data.empty()) {
		throw std::invalid_argument("timeslice must have flux data");
	}

	next.ID = ID;
	used_IDs.insert(ID);

	timeslices.push_back(next);
}

int timeline::get_unused_ID() const {
	// Find an unused ID (warning, may hang if you have 2^32
	// timeslices, good luck!)
	int proposed_ID;
	do {
		proposed_ID = random();
	} while (used_IDs.find(proposed_ID) != used_IDs.end());

	return proposed_ID;
}

void timeline::insert(timeslice & next) {
	insert(next, get_unused_ID());
}

template<typename T> std::vector<T> after_idx(
	const std::vector<T> & in, size_t idx) {

	std::vector<T> out(in.begin() + idx, in.end());
	return out;
}

template<typename T> std::vector<T> resize_and_subtract(
	const std::vector<T> & in, size_t idx) {

	std::vector<T> out(in.begin() + idx, in.end());
	T start_offset = out[0];

	for (T & cur_entry: out) {
		cur_entry -= start_offset;
	}

	return out;
}

void timeline::split(std::list<timeslice>::iterator & to_split,
	size_t first_byte_after, split_strategy strategy) {

	size_t MFM_end =
		to_split->sec_data.MFM_train_indices[first_byte_after];

	size_t flux_end =
		to_split->mfm_train.flux_indices[MFM_end];

	timeslice ts_before = *to_split, ts_after;
	ts_before.flux_data = std::vector(to_split->flux_data.begin(),
		to_split->flux_data.begin() + flux_end);
	ts_after.flux_data = after_idx(to_split->flux_data, flux_end);

	// This is such a hack... propagate into the respective methods
	// please. TODO. It ought also to be verified against something
	// that's decoded from scratch.

	ts_before.mfm_train = to_split->mfm_train;
	ts_before.mfm_train.data.resize(MFM_end);
	ts_before.mfm_train.flux_indices.resize(MFM_end);

	ts_before.sec_data = to_split->sec_data;
	ts_before.sec_data.decoded_data.resize(first_byte_after);
	ts_before.sec_data.errors.resize(first_byte_after);
	ts_before.sec_data.MFM_train_indices.resize(first_byte_after);

	// TODO copy over all the metadata...
	ts_after.mfm_train = to_split->mfm_train;
	ts_after.mfm_train.data = after_idx(
		to_split->mfm_train.data, MFM_end);
	ts_after.mfm_train.flux_indices = resize_and_subtract(
		to_split->mfm_train.flux_indices, MFM_end);

	ts_after.sec_data = to_split->sec_data;
	ts_after.sec_data.decoded_data = after_idx(
		to_split->sec_data.decoded_data, first_byte_after);
	ts_after.sec_data.errors = after_idx(
		to_split->sec_data.errors, first_byte_after);
	ts_after.sec_data.MFM_train_indices = resize_and_subtract(
		to_split->sec_data.MFM_train_indices, first_byte_after);

	ts_after.sector_data_begin = ts_before.sector_data_begin + first_byte_after;
	ts_after.mfm_train_begin = ts_before.mfm_train_begin + MFM_end;
	ts_after.flux_data_begin = ts_before.flux_data_begin + flux_end;

	int new_ID = get_unused_ID();

	switch(strategy) {
		case PRESERVE_FIRST:
			ts_before.status = to_split->status;
			ts_before.ID = to_split->ID;
			ts_after.ID = new_ID;
			break;
		case PRESERVE_LAST:
			ts_after.status = to_split->status;
			ts_before.ID = new_ID;
			ts_after.ID = to_split->ID;
			break;
		case PRESERVE_NEITHER:
			break; // they're already both set to TS_UNKNOWN, do nothing
	}

	// The code above is pretty hairy, so do some sanity checks
	// before we overwrite anything.

	assert(ts_before.flux_data_begin == to_split->flux_data_begin);
	assert(ts_before.mfm_train_begin == to_split->mfm_train_begin);
	assert(ts_before.sector_data_begin == to_split->sector_data_begin);

	assert(ts_after.flux_data_end() == to_split->flux_data_end());
	assert(ts_after.mfm_train_end() == to_split->mfm_train_end());
	assert(ts_after.sector_data_end() == to_split->sector_data_end());

	assert(ts_after.flux_data_begin == ts_before.flux_data_end());
	assert(ts_after.mfm_train_begin == ts_before.mfm_train_end());
	assert(ts_after.sector_data_begin == ts_before.sector_data_end());

	// Finally, insert and update the indices.

	used_IDs.insert(new_ID);
	*to_split = ts_after;
	std::list<timeslice>::iterator new_first_part =
		timeslices.insert(to_split, ts_before);
	update(new_first_part);
}

void timeline::update(std::list<timeslice>::iterator & to_update) {

	auto next = to_update, last = to_update;
	++next;

	while (next != timeslices.end()) {
		next->flux_data_begin = last->flux_data_end();
		next->mfm_train_begin = last->mfm_train_end();
		next->sector_data_begin = last->sector_data_end();

		++next;
		++last;
	}
}