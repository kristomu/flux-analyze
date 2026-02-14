#include "one_mode_dp.h"
#include "parity.h"
#include "tools.h"

#include <iostream>
#include <numeric>
#include <vector>

// Cut and paste code, welp.
// TODO: Fix that. Reimplement from scratch because this code isn't
// very old and already I think it's pretty ugly.

// perhaps some better class names?

// Inclusive parity is the parity of the sequence, i.e. number of
// bits emitted *including* the zeroes and trailing one due to the
// current pulse delay. Exclusive parity is the parity of the
// sequence up to, but not including, the current pulse delay.
class dp_idx {
	public:
		short clock;
		int idx;
		u_char inclusive_parity;
};

class dp_record {
	public:
		short zero_bits = -1;
		double penalty = 0;
		dp_idx previous_opt;
};

class dp_results {
	public:
		std::vector<int> clocks, zero_bits;
		std::vector<double> penalties;
		std::vector<int> valid;		// Does the current pules delay decode into
									// something that fails MFM rules?

		u_char starting_ex_parity; // starting exclusive parity
};

// Minimum and maximum clock index value.
const int MIN_CLOCK = 10, MAX_CLOCK = 50;
const int NUM_CLOCKS = MAX_CLOCK - MIN_CLOCK;

// How much increasing index by one increases the clock,
// e.g. CLOCK_GRANULARITY = 0.5 would make clock index 50
// represent a clock value of 50 * 0.5 = 25.
const double CLOCK_GRANULARITY = 1;

double fidelity_err(double x, size_t err_power) {

	if (err_power == 1) {
		return fabs(x)/fabs(100);
	}

	if (err_power == 2) {
		return sqr(x)/sqr(100);
	}

	return pow(fabs(x), err_power) / pow(100, err_power);
}

double smoothing_err(double x, size_t err_power) {

	if (err_power == 1) {
		return fabs(x)/fabs(100);
	}

	if (err_power == 2) {
		return sqr(x)/sqr(100);
	}

	return pow(fabs(x), err_power) / pow(100, err_power);
}

// State transition stuff.
// Is this the best way to do things? IDK
// Perhaps better to put it in a class later.

// ---- Hyperparameters ----

// Now I need some clean way to propose a state transition evolution.
// The simplest way to do so might be to divide this into two functions.
// The first function returns the penalty of doing the transition; and
// the second updates dp_record accordingly.
// Our inputs are the past and present dp_records (or past dp_record and
// present zero_bits???) as well as their dp_idx indices.

// This suggests that zero_bits is somehow special and that our
// abstraction is wrong at some point. But I don't know how to fix
// that at the moment.

// Since the fixed pattern is supposed to multiply the standard clock
// penalty by some factor, we also need to get the standard penalty
// (or rather, its current term) as an input.

/*class dp_record {
	public:
		short zero_bits = -1;
		short OOB_sequence_index = -1;
		double penalty = 0;
		dp_idx previous_opt;
};*/

// Everything about the current state should be set (i.e. copied from
// an admissible state, and index set to the proper bit number). Note that
// inclusive parity is specified ahead of time - we then work backwards to find
// every past state consistent with it.

// TBD.
/*
void dp_transition(dp_idx & current_state, int pulse_delay_length) {
	// Get the current pulse length based on the clock and delay,
	// as well as what the clock would be if the signal was
	// perfectly centered. The latter is used for penalty
	// calculations.

	double current_clock = current_state.clock * CLOCK_GRANULARITY;

	dp_record cur_record;

	size_t observed_pulse_delay = fluxes[i];
	cur_record.zero_bits = round(
		2.0 * observed_pulse_delay / current_clock) - 1;
	// Pulse lengths must be 1, 2, or 3.
	cur_record.zero_bits = std::min(3, std::max(1,
		(int)cur_record.zero_bits));

	double expected_pulse_delay =
		current_clock * (cur_record.zero_bits + 1)/ 2.0;

	// Measure the difference between what we observed and what we would
	// have expected if the clock was dead on.

	double pulse_delay_error = sqr(
		observed_pulse_delay - expected_pulse_delay);

	// If we're not the first flux delay, determine the least-penalty
	// previous state to use.

	if (i == 0) {
		cur_record.penalty = alpha * pulse_delay_error;

		dp[current_state.mode][current_state.clock]
			[current_state.idx][current_state.inclusive_parity] = cur_record;
		continue;
	}
}*/

int get_index(const dp_idx & x, int /*length*/) {
	return x.idx * NUM_CLOCKS * 2 + (x.clock-MIN_CLOCK) * 2 + x.inclusive_parity;
}

int get_max_index(int length) {
	return length * NUM_CLOCKS * 2 + (MAX_CLOCK - MIN_CLOCK) * 2 + 1;
}

dp_results do_dp(std::vector<int> fluxes, size_t len, double alpha,
	double in_band_error_penalty, size_t grace_period,
	size_t fidelity_err_power, size_t smoothing_err_power) {

	// Cut down the fluxes length so we can experiment with smaller
	// groups first.
	std::cout << "Flux length: " << fluxes.size() << " max length specified: " << len << "\n";
	if (fluxes.size() > len) {
		fluxes.resize(len);
	}

	/*std::vector<std::vector<std::vector<dp_record> > > dp(
			MAX_CLOCK+1, std::vector<std::vector<dp_record> >(
				len, std::vector<dp_record>(2)));
	*/
	std::vector<dp_record> dp_linear(get_max_index(len));

	// Get every possible index state so we can iterate over them
	// later - this saves a lot of nested for loops. Also keep
	// a track of possible states by parity, since valid past states
	// are limited to having the same parity.
	std::vector<dp_idx> possible_states;
	std::vector<std::vector<dp_idx> > possible_states_by_inc_parity(2);
	dp_idx current;

	for (current.clock = MIN_CLOCK; current.clock <= MAX_CLOCK;
		++current.clock) {
		for (current.inclusive_parity = 0; current.inclusive_parity < 2;
			++current.inclusive_parity) {

			possible_states.push_back(current);
			possible_states_by_inc_parity[
				current.inclusive_parity].push_back(current);
		}
	}

	for (size_t i = 0; i < fluxes.size(); ++i) {
		dp_record cur_record;
		for (dp_idx current_state: possible_states) {
			current_state.idx = i;

			// Get the current pulse length based on the clock and delay,
			// as well as what the clock would be if the signal was
			// perfectly centered. The latter is used for penalty
			// calculations.

			double current_clock = current_state.clock * CLOCK_GRANULARITY;

			//dp_record cur_record;

			size_t observed_pulse_delay = fluxes[i];
			cur_record.zero_bits = round(
				2.0 * observed_pulse_delay / current_clock) - 1;
			// Pulse lengths must be 1, 2, or 3.
			cur_record.zero_bits = std::min(3, std::max(1,
				(int)cur_record.zero_bits));

			double expected_pulse_delay =
				current_clock * (cur_record.zero_bits + 1)/ 2.0;

			// Measure the difference between what we observed and what we would
			// have expected if the clock was dead on.

			double pulse_delay_error = fidelity_err(
				observed_pulse_delay - expected_pulse_delay,
				fidelity_err_power);

			// If we're not the first flux delay, determine the least-penalty
			// previous state to use.

			if (i == 0) {
				cur_record.penalty = alpha * pulse_delay_error;

				dp_linear[get_index(current_state, len)] = cur_record;
				continue;
			}

			// Our current state has a particular inclusive parity, and
			// should thus be matched up with a past state that has the same
			// *inclusive* parity as our *exclusive* parity, so that adding
			// our emitted bits produces the given inclusive parity.

			u_char exclusive_parity = get_exclusive_parity(
				cur_record.zero_bits, current_state.inclusive_parity);

			bool primed = false;
			double best_penalty;
			dp_idx best_idx;
			dp_record last_record;

			for (dp_idx last_state: possible_states_by_inc_parity[exclusive_parity]) {

				last_state.idx = i-1;
				last_record = dp_linear[get_index(last_state, len)];

				// Calculate the penalty when transitioning from this
				// state. (I may move this to another function later.)

				double last_clock = last_state.clock * CLOCK_GRANULARITY;

				double candidate_penalty = alpha * pulse_delay_error +
					(1-alpha) * smoothing_err(current_clock - last_clock,
						smoothing_err_power);

				// XXX: Biasing very slightly away from MFM errors is
				// enough to break incorrectly classified/mass error
				// near-ties as seen on MS_Plus_disk3_OK_track with ~10k pulses,
				// but I want to do it in a more principled manner,
				// so that's commented out for now.

				if (!is_mfm_valid(cur_record.zero_bits,
					exclusive_parity) && i > grace_period) {
					candidate_penalty *= in_band_error_penalty;
				}

				// Add the accumulated penalty from the last state.
				candidate_penalty += last_record.penalty;

				if (!primed || candidate_penalty <= best_penalty) {
					best_penalty = candidate_penalty;
					best_idx = last_state;
					primed = true;
				}
			}

			if (!primed) {
				throw std::logic_error("DP: Could not find a prior state to build on!");
			} else {
				cur_record.penalty = best_penalty;
				cur_record.previous_opt = best_idx;

				dp_linear[get_index(current_state, len)] = cur_record;
			}
		}
	}

	// ============= Get the solution. ========

	// First determine the least penalty solution.

	double best_penalty = std::numeric_limits<double>::infinity();
	dp_idx best_idx;

	for (dp_idx current_state: possible_states) {
		current_state.idx = fluxes.size()-1;

		double candidate_penalty = dp_linear[
			get_index(current_state, len)].penalty;

		if (candidate_penalty <= best_penalty) {
			best_idx = current_state;
			best_penalty = candidate_penalty;
		}
	}

	std::cout << "Optimal penalty: " << best_penalty << "\n";

	std::vector<dp_idx> optimal_path;
	std::vector<dp_record> optimal_values;

	// Trace the solution back through the DP array

	for (size_t i = 0; i < fluxes.size(); ++i) {
		optimal_path.push_back(best_idx);

		dp_record at_current = dp_linear[get_index(best_idx, len)];

		optimal_values.push_back(at_current);

		best_idx = at_current.previous_opt;
	}

	// And reverse the order of that path to get beginning to end.

	std::reverse(optimal_path.begin(), optimal_path.end());
	std::reverse(optimal_values.begin(), optimal_values.end());

	dp_results out;
	out.starting_ex_parity = get_exclusive_parity(
		optimal_values[0].zero_bits,
		optimal_path[0].inclusive_parity);

	for (size_t i = 0; i < optimal_values.size(); ++i) {
		out.clocks.push_back(optimal_path[i].clock);
		out.zero_bits.push_back(optimal_values[i].zero_bits);
		out.penalties.push_back(optimal_values[i].penalty);

		int exclusive_parity = get_exclusive_parity(
			optimal_values[i].zero_bits,
			optimal_path[i].inclusive_parity);

		if (is_mfm_valid(optimal_values[i].zero_bits,
			exclusive_parity)) {
			out.valid.push_back(1);
		} else {
			out.valid.push_back(0);
		}
	}

	return out;
}

/*void print_dp_results(const dp_results & results) {

	// This shows a vast field of ones interspersed with
	// 1 0 1 1 1 1 0 1 1 1 1 0 1, which I think is the intentional MFM
	// errors of the OOB A1A1A1/C2C2C2 preambles.

	// But try the MS Plus disk3 OK track (track 4 side 0) with 10000
	// entries (or even just 500) and it flips out.

	std::cout << "Invalid and valid MFM bits (0 and 1 resp.):\n";
	print_vector(results.valid);
	
	std::cout <<"\n\nEstimated clock rate:\n";
	print_vector(results.clocks);

	std::cout << "Clock variance: " << variance(results.clocks) << "\n";

	std::cout << "Starting exclusive parity: ";
	switch(results.starting_ex_parity) {
		case 0: std::cout << "even.\n"; break;
		case 1: std::cout << "odd.\n"; break;
		default:
			throw std::logic_error("print_dp_results: Invalid parity value!");
	}
}*/

// Decode a stream of pulse delays (zero bit counts) and output which
// pulses create an MFM violation.

// We need to output pulses that create MFM violations rather than the
// bits that are OK, because the on-the-fly calculation using parity
// also works on a pulse-by-pulse basis.
std::vector<int> valid_from_zero_bits(const std::vector<int> & zero_bits,
	int starting_ex_parity) {

	std::vector<int> expanded_bit_train, responsible_pulse;

	// If the initial parity is odd, then pull in the preceding
	// trailing edge. (I.e. we assume the number of bits clocked is 1,
	// not -1, if the initial parity is odd.)

	if (starting_ex_parity == 1) {
		expanded_bit_train.push_back(1);
		responsible_pulse.push_back(-1);
	}

	size_t i;

	for (i = 0; i < zero_bits.size(); ++i) {
		// Add the leading zero bits
		for (int j = 0; j < zero_bits[i]; ++j) {
			expanded_bit_train.push_back(0);
			responsible_pulse.push_back(i);
		}
		// Then add the trailing edge.
		expanded_bit_train.push_back(1);
		responsible_pulse.push_back(i);
	}

	// Now decode.
	bool first_bit = true;
	int last_value = 0;

	std::vector<int> valid_bit, valid_pulse(zero_bits.size(), 1);

	for (i = 0; i < expanded_bit_train.size() - 1; i += 2) {
		// Get most and least significant bit of the pair.
		int a = expanded_bit_train[i],
			b = expanded_bit_train[i+1];

		int pulse_idx = std::max(0, responsible_pulse[i]);

		if (a == 0 && b == 0) {
			// Valid if last decoded bit was a one;
			// produces a zero.
			if (first_bit || last_value == 1) {
				valid_bit.push_back(1);
			} else {
				valid_bit.push_back(0);
				valid_pulse[pulse_idx] = 0;
			}
			last_value = 0;
		}

		if (a == 0 && b == 1) {
			// Always valid, and generates a 1
			valid_bit.push_back(1);
			last_value = 1;
		}

		if (a == 1 && b == 0) {
			// Valid if last decoded bit was a zero;
			// produces a zero.
			if (first_bit || last_value == 0) {
				valid_bit.push_back(1);
			} else {
				valid_bit.push_back(0);
				valid_pulse[pulse_idx] = 0;
			}
			last_value = 0;
		}

		if (a == 1 && b == 1) {
			throw std::invalid_argument("illegal MFM bit pair");
		}

		first_bit = false;
	}

	return valid_pulse;
}

// TODO: results.valid are by *pulse*, inferred_valid is by *bit pair*.
// Fix that.
void test_dp_validity_consistent(const dp_results & results) {

	std::vector<int> inferred_valid = valid_from_zero_bits(
		results.zero_bits, results.starting_ex_parity);

	if (inferred_valid.size() != results.valid.size()) {
		std::cout << "Inferred valid data is over " <<
			inferred_valid.size() << " bit pairs, but on-the-fly" <<
			" is over " << results.valid.size() << "\n";
		return;
	}

	for (size_t i = 0; i < inferred_valid.size(); ++i) {
		if (inferred_valid[i] != results.valid[i]) {
			std::cout << "MFM validity test: mismatch between on-the-fly "
				"and manual calc at " << i << "\n";
			return;
		}
	}

	std::cout << "Validity is consistent.\n";
}

MFM_train_data QND_one_mode_DP::get_MFM_train(
	const std::vector<int> & fluxes, size_t start_pos,
	size_t end_pos, double & error_out) const {

	error_out = std::numeric_limits<double>::infinity();

	size_t period_length = 400000;

	MFM_train_data whole_train;

	// Do a piecewise DP calculation because it's so
	// memory-heavy that we can't fit everything in one
	// go. NOTE: There may be off-by-one errors here.

	for (size_t offset = start_pos; offset < end_pos; offset += period_length) {

		std::cout << "Progress: " << (offset-start_pos)/(double)(end_pos-start_pos) << std::endl;

		std::vector<int> cut_fluxes(fluxes.begin() + offset,
			fluxes.begin() + std::min(end_pos, offset + period_length));

		dp_results res = do_dp(cut_fluxes, period_length,
			alpha, beta, grace_period, fidelity_err_power, smoothing_err_power);

		MFM_train_data out_train;

		if (offset == start_pos) {
			out_train.data.push_back(1);
			out_train.flux_indices.push_back(offset);
		}

		for (size_t i = 0; i < res.zero_bits.size(); ++i) {

			for (int j = 0; j < res.zero_bits[i]; ++j) {
				out_train.data.push_back(0);
				out_train.flux_indices.push_back(offset+i);
			}

			out_train.data.push_back(1);
			out_train.flux_indices.push_back(offset+i);
		}

		if (do_only_one) {
			return out_train;
		}

		whole_train += out_train;
	}

	return whole_train;
}
