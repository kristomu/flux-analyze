#include "variable_offset.h"
#include <iostream>

double dot(const std::vector<double> & a,
	const std::vector<double> & b) {

	double out = 0;

	for (size_t i = 0; i < a.size(); ++i) {
		out += a[i] * b[i];
	}

	return out;
}

void matmul(const std::vector<std::vector<double> > & a,
	const std::vector<std::vector<double> > & b,
	std::vector<std::vector<double> > & out) {

	for (size_t row = 0; row < a.size(); ++row) {
		for (size_t column = 0; column < a[0].size(); ++column) {
			out[row][column] = 0;

			for (size_t element = 0; element < b.size(); ++element) {
				out[row][column] += a[row][element] * b[element][column]; 
			}
		}
	}
}

// Outputs another vector, out = matrix * vec^T. 
// This will be expensive, but I can optimize later.
std::vector<double> matvec(const std::vector<std::vector<double> > & matrix,
	std::vector<double> & vec) {

	std::vector<double> out(matrix.size(), 0);

	for (size_t row = 0; row < matrix.size(); ++row) {
		for (size_t col = 0; col < vec.size(); ++col) {
			out[row] += matrix[row][col] * vec[col];
		}
	}

	return out;
}

std::vector<double> matvec(std::vector<double> & vec,
	const std::vector<std::vector<double> > & matrix) {

	std::vector<double> out(matrix.size(), 0);

	for (size_t row = 0; row < vec.size(); ++row) {
		for (size_t col = 0; col < matrix.size(); ++col) {
			out[row] += matrix[col][row] * vec[col];
		}
	}

	return out;
}

double sq_norm(const std::vector<double> & a) {
	return dot(a, a);
}

MFM_train_data GD_variable_offset_decoder::get_MFM_train(
	const std::vector<int> & fluxes, size_t start_pos,
	size_t end_pos, double & error_out) const {

	MFM_train_data train;

	// state space multiplier
	std::vector<double> h = {1, 1};

	// state space
	// the 2/x is because we want half-clocks.
	std::vector<double> theta = {2/22.0, 0}; // prior, fix later

	train.data.push_back(1);
	train.flux_indices.push_back(start_pos);

	for (size_t i = start_pos; i < end_pos; ++i) {
		size_t observed_pulse_delay = fluxes[i];

		h[0] = observed_pulse_delay;

		double estimated_num_half_clocks = dot(h, theta);

		// we know that actual number of half-clocks must be
		// an integer, so round off and make that y(t).
		// This is the number of zero bits, plus one.

		double y = round(estimated_num_half_clocks);

		double error = y - estimated_num_half_clocks;

		// update theta
		for (int j = 0; j < 2; ++j) {
			theta[j] = theta[j] + stepsize/(1e-6 + sq_norm(h)) * error * h[j];
		}

		int zero_bits = y - 1;

		zero_bits = std::min(3, std::max(1, zero_bits));

		for (int j = 0; j < zero_bits; ++j) {
			train.data.push_back(0);
			train.flux_indices.push_back(i + start_pos);
		}

		train.data.push_back(1);
		train.flux_indices.push_back(i + start_pos);
	}

	error_out = std::numeric_limits<double>::infinity();

	return train;
}