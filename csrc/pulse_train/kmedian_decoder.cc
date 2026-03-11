#include "kmedian_decoder.h"

MFM_train_data kmedian_decoder::get_MFM_train(
	const std::vector<int> & fluxes, size_t start_pos,
	size_t end_pos, double & RMSE_out) const {

	KCluster kmedian;

	MFM_train_data train;

	train.data.push_back(1);
	train.flux_indices.push_back(start_pos);

	RMSE_out = 0;

	std::vector<double> clusters = kmedian.optimalKMedian(
		fluxes, 3);

	for (size_t i = start_pos; i < end_pos; ++i) {
		double distance_to_record = std::numeric_limits<double>::infinity();
		size_t closest_cluster = -1;

		for (size_t cluster = 0; cluster < clusters.size(); ++cluster) {
			double distance_to_cluster = fabs(
				fluxes[i]-clusters[cluster]);

			if (distance_to_cluster < distance_to_record) {
				distance_to_record = distance_to_cluster;
				closest_cluster = cluster;
			}
		}

		RMSE_out += distance_to_record * distance_to_record;

		size_t zeroes = closest_cluster + 1;

		for (size_t j = 0; j < zeroes; ++j) {
			train.data.push_back(0);
			train.flux_indices.push_back(i);
		}
		train.data.push_back(1);
		train.flux_indices.push_back(i);

	}

	RMSE_out = std::sqrt(RMSE_out / (end_pos - start_pos));

	return train;
}
