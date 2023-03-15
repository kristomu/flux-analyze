#include "pulse_train.h"
#include "sector_data.h"
#include "timeline.h"

#include <assert.h>

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

	timeslices.push_back(next);
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
	size_t first_byte_after, slice_status before_status,
	slice_status after_status) {

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

	ts_before.status = before_status;
	ts_after.status = after_status;

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

	// Finally, insert.

	*to_split = ts_after;
	timeslices.insert(to_split, ts_before);
}