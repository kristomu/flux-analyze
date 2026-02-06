#pragma once

#include "flux_decoder.h"

class constant_clock_decoder: public pulse_decoder {
	private:
		double clock;

	public:
		// Do I want to do this everywhere, or just do a lazy
		// two-step function approach instead? Hmm...
		using pulse_decoder::get_MFM_train;

		void set_clock(double clock_in) { clock = clock_in; }

		MFM_train_data get_MFM_train(
			const std::vector<int> & fluxes, size_t start_pos,
			size_t end_pos, double & error_out) const;
};

class orig_EWMA_causal_clock_decoder : public pulse_decoder {
	private:
		double alpha;
		double initial_clock;

	public:
		using pulse_decoder::get_MFM_train;

		orig_EWMA_causal_clock_decoder(double alpha_in,
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

		orig_EWMA_causal_clock_decoder() : orig_EWMA_causal_clock_decoder(0.5, 24) {}

		MFM_train_data get_MFM_train(
			const std::vector<int> & fluxes, size_t start_pos,
			size_t end_pos, double & error_out) const;
};

class historical_EWMA_decoder : public pulse_decoder {
	private:
		double alpha;
		double initial_clock;

	public:
		using pulse_decoder::get_MFM_train;

		historical_EWMA_decoder(double alpha_in,
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

		historical_EWMA_decoder() : historical_EWMA_decoder(0.5, 24) {}

		MFM_train_data get_MFM_train(
			const std::vector<int> & fluxes, size_t start_pos,
			size_t end_pos, double & error_out) const;
};
