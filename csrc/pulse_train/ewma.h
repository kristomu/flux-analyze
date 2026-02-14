#pragma once

#include "flux_decoder.h"

#include <vector>

// These should be "private", i.e. not exposed outside these classes.
// So there's no need to place them in a header accessible to all.
/*
std::vector<double> causal_get_approximate_clock(
	const std::vector<int> & fluxes, double alpha,
	double clock_prior, bool use_prior);

std::vector<double> causal_get_approximate_clock(
	const std::vector<int> & fluxes, double alpha);

std::vector<double> acausal_get_approximate_clock(
	const std::vector<int> & fluxes, double alpha);
*/

class causal_EWMA_clock_decoder : public pulse_decoder {
	private:
		double alpha;
		double initial_clock = -1;

		bool is_initial_clock_set() const { return initial_clock > 0; }

	public:
		using pulse_decoder::get_MFM_train;

		causal_EWMA_clock_decoder(double alpha_in,
			double initial_clock_in) {

			alpha = alpha_in;
			initial_clock = initial_clock_in;
		}

		void set_alpha(double alpha_in) {
			alpha = alpha_in;
		}

		void set_initial_clock(double initial_clock_in) {
			initial_clock = initial_clock_in;
		}

		void clear_initial_clock() {
			initial_clock = -1;
		}

		causal_EWMA_clock_decoder() : causal_EWMA_clock_decoder(0.5, 24) {}

		MFM_train_data get_MFM_train(
			const std::vector<int> & fluxes, size_t start_pos,
			size_t end_pos, double & error_out) const;

};

// This one doesn't work; the causal decoder is better on severely
// degraded disks. (It might be better to use something more fancy,
// like simultaneous equation solving with EM or something.. I may
// test that later.)

class acausal_EWMA_clock_decoder : public pulse_decoder {
	private:
		double alpha;

	public:
		using pulse_decoder::get_MFM_train;

		acausal_EWMA_clock_decoder(double alpha_in) {
			alpha = alpha_in;
		}

		void set_alpha(double alpha_in) {
			alpha = alpha_in;
		}

		acausal_EWMA_clock_decoder() : acausal_EWMA_clock_decoder(0.5) {}

		MFM_train_data get_MFM_train(
			const std::vector<int> & fluxes, size_t start_pos,
			size_t end_pos, double & error_out) const;
};
