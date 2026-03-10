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

		double max_relative_change = -1;

		bool ignore_lower = false;
		bool clamp = true;

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

		void set_ignore_lower(bool ignore_lower_in) {
			ignore_lower = ignore_lower_in;
		}

		void set_clamp(bool clamp_in) {
			clamp = clamp_in;
		}

		void set_max_relative_change(double max_relative_change_in) {
			max_relative_change = max_relative_change_in;
		}

		causal_EWMA_clock_decoder() : causal_EWMA_clock_decoder(0.5, 24) {}

		MFM_train_data get_MFM_train(
			const std::vector<int> & fluxes, size_t start_pos,
			size_t end_pos, double & error_out) const;

		// Initial clock is not considered a parameter here because
		// we can do without it, and fixing an initial clock wouldn't
		// generalize well. (TODO: Is this the right call?)

		void set_params(const std::vector<double> & params) {
			alpha = params[0];
			ignore_lower = (int)params[1];
			clamp = (int)params[2];
			max_relative_change = params[3];
		}

		std::vector<param_t> get_parameter_types() const {
			return {PARAM_LOG_REAL, PARAM_INTEGER, PARAM_INTEGER, PARAM_LOG_REAL};
		}

		std::vector<double> get_parameter_min() const {
			return {0, 0, 0, 1};
		}

		std::vector<double> get_parameter_max() const {
			return {1, 1, 1, 3};
		}

		std::vector<double> get_current_params() const {
			return {alpha, (double)((int)ignore_lower),
							(double)((int)clamp), max_relative_change};
		}

		std::string name() const {
			return "Causal EWMA decoder";
		}
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

		std::vector<param_t> get_parameter_types() const {
			return {PARAM_LOG_REAL};
		}

		std::vector<double> get_parameter_min() const {
			return {0};
		}

		std::vector<double> get_parameter_max() const {
			return {1}; // For reasonable floppy types?
		}

		std::vector<double> get_current_params() const {
			return {alpha};
		}

		void set_params(const std::vector<double> & params) {
			alpha = params[0];
		}

		std::string name() const {
			return "Acausal EWMA decoder";
		}
};
