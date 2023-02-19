#include <iostream>
#include <vector>

double sqr(double x) {
	return x*x;
}

// Returns either an error measure or denoises a flux signal, both according
// to the given clock estimate. See comments below about the algorithm used.
// Passing parameters like this is kind of ugly, but since there's so much
// repeated logic, doing it as two functions would be worse.

// For some reason, it's now 2x slower.

// And it doesn't actually dewarp the warped MS Plus floppy. I'll have to
// check if this can be salvaged, although I'm more likely to just steal
// the PLL dewarp from Python with various guesses at the clock estimate.

double test_dewarp(const std::vector<int> & fluxes, double clock_estimate,
	bool return_denoised_signal, std::vector<int> & denoised_signal) {
	// Try to remove low frequency variation in fluxes (due to uneven surface
	// or similar). This uses the model that

	// x[i] = p[i] * clock_estimate/2 + e[i]

	// where x[i] is the observed flux delay, p[i] is the half-clock multiple
	// (1, 2 or 3 for MFM), and e[i] is an error term. We then optimize the
	// choices of p so as to minimize
	
	// sum (e[i] - e[i-1])^2

	// by dynamic programming, thus minimzing the high frequency components
	// of the error term e.

	// But not bad for a first attempt. Kinda slow, though...
	// Maybe a better model is: x[i] * e[i] = p[i] * clock/2, because
	// then e[i] directly models the "stretching" of x[i] itself. But that
	// doesn't work in the other direction. It rewards small clocks, which is
	// not surprising because multiplying down costs less than multiplying up.
	// (compare e.g. ||1/2||^2 to 2^2.)

	// TODO: Find out a principled model that agrees with my idea of what
	// warping actually *is*.

	// The model needs to be
	//		- intuitive
	//		- unbiased between negative/too slow and positive/too fast
	//		- independent of clock scale.

	// Another problem is that I don't really know the frequency response
	// of (e[i] - e[i-1])^2. Maybe it's too heavy and adjusts clock_estimate
	// to filter out DC? One alternative is something like |e[i] - e[i-1]| if
	// < half clock, otherwise +inf, but that's hardly robust, is it?

	// We'll also need an e[i] term in there somewhere so that if there's no
	// warping at all, a more accurate fit is preferred. But it's not clear
	// how to do that without introducing a hyperparameter. I'll have to think
	// more about it.

	bool first = true;

	struct penalty_record {
		double penalty_so_far;
		double e_i;
		int p_i_choice;
		int last_p_i_choice;
	};

	const int min_p_i = 2, max_p_i = 6;

	// The dynamic programming mechanism needs a record of the
	// penalties (and associated error values) for the last round
	// as well as for the current round, so allocate two rows of penalty
	// records - one for the old and one for the new.

	// Actually getting the e_i values will require a full array
	// of max_p_i * fluxes.size(), so if return_denoised_signal is
	// true, that's what we do. Due to RAM usage, it's much slower!

	int num_rows = 2;
	if (return_denoised_signal) {
		num_rows = fluxes.size();
	}

	std::vector<std::vector<penalty_record> > penalties(num_rows,
		std::vector<penalty_record>(max_p_i-min_p_i));

	// These are indices into the penalties vector. We'll
	// swap them around each time so as not to have to
	// waste CPU power on a bunch of reallocations.

	int old_idx = 1, new_idx = 0;
	penalty_record our_best;

	for (size_t i = 0; i < fluxes.size(); ++i) {
		double x_i = fluxes[i];

		// TODO: Somehow be more smart about how many p_i values we
		// need to check, by making an observation that after some
		// point, (p[i] + e[i]) * clock_estimate/2 can only increase
		// e[i] away from e[i-1], hence worsening the penalty; the
		// problem is that we don't know what e[i+1] will be...
		for (int p_i_choice = min_p_i; p_i_choice < max_p_i; ++p_i_choice) {
			// Trying model x[i] - e[i] = p[i] * clock/2.
			double e_i = x_i - (clock_estimate * p_i_choice) / 2.0;
			// Evaluate the combined error by choosing the best
			// error for the last term. If it's the first flux delay,
			// everything is free (because we don't know what e_-1 is).
			double record_penalty = std::numeric_limits<double>::max();
			int recordholder = -1; // <-- this makes it slow!!!

			if (first) {
				record_penalty = 0;
			} else {
				for (const penalty_record & p: penalties[old_idx]) {
					double differential_penalty = sqr(p.e_i - e_i);
					double candidate_penalty = differential_penalty +
						sqr(e_i * 0.01) + p.penalty_so_far;

					if (candidate_penalty < record_penalty) {
						record_penalty = candidate_penalty;
						recordholder = p.p_i_choice;
					}
				}
			}

			our_best.penalty_so_far = record_penalty;
			our_best.e_i = e_i;
			our_best.p_i_choice = p_i_choice;
			our_best.last_p_i_choice = recordholder;
			penalties[new_idx][p_i_choice-min_p_i] = our_best;

			//std::cout << " i = " << i << " p_i choice = " << p_i_choice << " e_i = " << e_i << " r.p. " << record_penalty << std::endl;
		}

		if (first) {
			first = false;
		}

		if (return_denoised_signal) {
			old_idx = new_idx;
			new_idx++;
		} else {
			std::swap(old_idx, new_idx);
		}
	}

	penalty_record best = penalties[old_idx][0];
	for (penalty_record p: penalties[old_idx]) {
		if (p.penalty_so_far < best.penalty_so_far) {
			best = p;
		}
	}

	// Reconstruct the denoised signal if told to do so.
	if (return_denoised_signal) {
		// We need to start from the end and trace the path in
		// reverse, then reverse the output when done.
		denoised_signal.resize(0);

		for (size_t i = fluxes.size(); i > 0; --i) {
			denoised_signal.push_back(fluxes[i-1] - best.e_i);
			best = penalties[i-1][best.last_p_i_choice - min_p_i];
		}

		std::reverse(denoised_signal.begin(), denoised_signal.end());
	}

	return sqrt(best.penalty_so_far);
}

double test_dewarp(const std::vector<int> & fluxes, double clock_estimate) {
	std::vector<int> throwaway;

	return test_dewarp(fluxes, clock_estimate, false, throwaway);
}

std::vector<int> get_dewarped(const std::vector<int> & fluxes,
	double clock_estimate) {

	std::vector<int> denoised_signal;
	test_dewarp(fluxes, clock_estimate, true, denoised_signal);

	return denoised_signal;
}