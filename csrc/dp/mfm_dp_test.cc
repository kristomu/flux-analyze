// Experimenting with dynamic programming as a way to recover the data stream.

/*	My previous plan was a bit ambitious. Now I'll try to implement two different
	DP algorithms: one slow but reasonably sure, and one quick and approximate.
	They'll both assume MFM data (i.e. there's only one mode) with a penalty to
	MFM errors. (I can add fixed OOB sequences later; in a perhaps funny instance
	of "what's easy is hard, what's hard is easy", doing fixed sequence search
	naively is very slow; we need another dimension to the DP array with length
	equal to the longest needle/search string).

	Old comments below, adapted:

	The dynamic programming approach would read the whole flux stream for every
	possible clock at once, slowly varying the allowed clock according to a penalty
	function that punishes the pulse delay from being far from what's expected and
	for the current clock to change too much from the previous pulse's clock.

	This, as such, isn't too different from just doing a moving or exponential
	average. But I've then intended to make the DP system prioritize that all
	data must obey MFM rules.

	However, I think braiding all of this into a DP algorithm could get very messy.
	So beware, all who read this. I may need multiple iterations to get it to be
	somewhat presentable.
*/

/*	The current dynamic programming idea is as follows, where * indicates something
	that hasn't been implemented yet.

	The mode is:
		  0 - in-band

	The state for each step is:
		INDEX:	Current mode
		INDEX:	Current estimated clock value
		INDEX:	Current index into list of flux transitions (ith transition)
		INDEX:	MFM decoding parity

		STATE:	Current flux delay period (number of zero bits)
		STATE:	accumulated penalty/cost
		STATE:	previous optimal entry

	where INDEX is a dimension and STATE is a value at a concrete array cell.

	MFM decoding parity is pretty tricky and needs a more detailed explanation.
	I'll give that explanation later.
*/

/*	We're assuming the floppy is formatted with a standard MFM encoding, where 0
	indicates a lack of flux reversal and 1 is its presence:

		MFM sequence	last bit	current bit
		01				N/A			1
		00				1			0
		10				0			0

	and that flux transitions are based off a clock as follows:

		clock multiple		sequence
		1					01
		1.5					001
		2					0001
*/

#include "flux_record.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>

/*	Detecting what MFM transitions are allowed is pretty easy from the pairs of sequence
	bits alone. Assuming the first pair has already been verified as being part of an
	allowed sequence, we have

		XX 01	is valid (01 encodes for a 1 bit and can be preceded by anything)
		01 00	is valid, but no other prefix is valid
		00 10	is valid, and
		10 10	is valid,

	but no other combination is valid.

	However, in order to make in-band give a penalty on MFM transitions that aren't
	allowed, we need to keep track of when something disallowed happens even as each
	decoded pulse sequence contributes a different number of bits. Consider, for
	instance, an incoming 0001 sequence:

			prior bits    | new bits
			... 01 00 01  | 00 01
			even parity	  | even parity

	This is allowed, but

			prior bits    | new bits
			.... 01 00 1	0 00 1?
			odd parity    | odd parity

	here our bitstream stops at an odd number of bits, so we can't tell from the
	prior bits alone if everything is allowed. There is a trick, though; because
	every (legal) pulse length is at least one long, it produces at least one
	zero bit before the trailing edge. So we must have

			.... 01 00 1	0 00 10

	I'll put more information in a separate parity document, but the upshot
	is that, as long as we have no much-too-long or much-too-short flux
	transitions, then an MFM error occurs *only* at odd parity with a pulse
	length of three.
*/

double sqr(double x) { return x*x; }

// Debug and diagnostic functions.

template<typename T> void print_vector(const std::vector<T> & vec) {
	std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(std::cout, " " ));
	std::cout << "\n";
}

double variance(const std::vector<int> & vec) {
	double mean = std::accumulate(vec.begin(), vec.end(), 0) / (double)vec.size();
	double var = 0;

	for (int k: vec) {
		var += sqr(k-mean);
	}

	return var;
}

// Checks if the sequence 1, a zeroes, 1, b zeroes, breaks any rules, when
// the sequence 1, a zeroes has already been tested and can be assumed to
// not have broken any rules. This is done by padding with a one on the left
// and a zero on the right if required. Parity is from the perspective of the
// bits so far, including the trailing edge but not any of the zero bits
// counted by zero_bits nor the present pulse's trailing edge.

// See parity.md for more info.
bool is_mfm_valid(short zero_bits, u_char exclusive_parity) {
	// (Trivia: it actually doesn't matter if we use exclusive or
	// inclusive parity here because any pulse with three zero bits
	// emits four in total, and 4 mod 2 == 0, so parity doesn't
	// change. That said, doing so would be supremely confusing.)

	return !(exclusive_parity == 1 && zero_bits == 3);
}

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

// Given an inclusive parity, get the corresponding exclusive
// parity.
u_char get_exclusive_parity(short zero_bits,
	u_char inclusive_parity) {

	// The +1 is for the trailing edge's one bit.
	int signed_parity = inclusive_parity - (zero_bits + 1);

	while (signed_parity < 0) {
		signed_parity += 2;
	}

	return signed_parity % 2;
}

u_char get_inclusive_parity(short zero_bits,
	u_char exclusive_parity) {

	return (exclusive_parity + zero_bits + 1) % 2;
}

class dp_results {
	public:
		std::vector<int> clocks, zero_bits;
		std::vector<double> penalties;
		std::vector<int> valid;		// Does the current pules delay decode into
									// something that fails MFM rules?

		u_char starting_ex_parity; // starting exclusive parity
};

// Minimum and maximum clock index value.
const int MIN_CLOCK = 1, MAX_CLOCK = 100;
const int NUM_CLOCKS = MAX_CLOCK - MIN_CLOCK;

// How much increasing index by one increases the clock,
// e.g. CLOCK_GRANULARITY = 0.5 would make clock index 50
// represent a clock value of 50 * 0.5 = 25.
const double CLOCK_GRANULARITY = 1;

// State transition stuff.
// Is this the best way to do things? IDK
// Perhaps better to put it in a class later.

// ---- Hyperparameters ----

// How much does an MFM error cost when in the in-band mode?
const double IN_BAND_ERROR_PENALTY = 1;

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

dp_results do_dp(const flux_record & f, size_t len) {

	// Noise sensitivity hyperparameter.
	// Higher values of this allows greater changes in clock
	// from one pulse to the next, while lower values prioritize fitting
	// the set clock to the current pulse delay.
	// The range is 0..1.
	double alpha = 0.5;

	// Cut down the fluxes length so we can experiment with smaller
	// groups first.
	std::vector<int> fluxes = f.fluxes;
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

			double pulse_delay_error = sqr(
				observed_pulse_delay - expected_pulse_delay);

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
					(1-alpha) * sqr(current_clock - last_clock);

				// XXX: Biasing very slightly away from MFM errors is
				// enough to break incorrectly classified/mass error
				// near-ties as seen on MS_Plus_disk3_OK_track with ~10k pulses,
				// but I want to do it in a more principled manner,
				// so that's commented out for now.

				if (!is_mfm_valid(cur_record.zero_bits,
					exclusive_parity)) {
					candidate_penalty += IN_BAND_ERROR_PENALTY;
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

void print_dp_results(const dp_results & results) {

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
}

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

/*
std::vector<int> fluxes = f.fluxes;
	std::cout << "Flux length: " << fluxes.size() << " max length specified: " << len << "\n";
	if (fluxes.size() > len) {
		fluxes.resize(len);
	}
*/

void periodicity_test(const flux_record & /*f*/) {
	// Quick and dirty periodicity test. Take the differences between consecutive
	// flux delays, round them to nearest 10, and then look for rare patterns
	// that recur to find a possible period for the floppy signal.

	// With 300 RPM (0.2 secs per revolution), and with a clock rate of 500 kHz,
	// and one FE unit being 83 ns, we'd expect a running sum of about
	// 2.4 million between each loop. Suppose the normal clock value is 24.
	// Then that's 100k values (transitions).

	// (My DP can do that on an i7 in about three seconds now. That's a little
	// slow, but ideas...)

	// TODO: Actually implement that! :-P
}

// Get an approximate clock rate by using exponential weighted averaging
// with alpha=0.5. I think that's what a real FDD does; I might also be able
// to use this to narrow down the clock range that has to be investigated
// in the DP.
// The int value is -1 for the first 300 or so delays to signify that the
// EWMA has to "warm up" first.

std::vector<int> get_approximate_clock(const flux_record & f) {
	std::vector<int> fluxes = f.fluxes;
	std::vector<int> appx_clock;
	std::cout << "[Appx clock] Flux length: " << fluxes.size() << "\n";	

	double mean_clock = 1; // or insert your prior here
	size_t warmup_period = 300;

	for (size_t i = 0; i < fluxes.size(); ++i) {
		size_t observed_pulse_delay = fluxes[i];
		
		size_t zero_bits = round(
			2.0 * observed_pulse_delay / mean_clock) - 1;
		zero_bits = std::min(3, std::max(1, (int)zero_bits));

		double new_clock = 2 * observed_pulse_delay / (zero_bits + 1);

		mean_clock = 0.5 * mean_clock + 0.5 * new_clock;

		if (i < warmup_period) {
			appx_clock.push_back(-1);
		} else {
			appx_clock.push_back(round(mean_clock));
		}
	}

	return appx_clock;
}

int main(int argc, char ** argv) {

	std::string flux_filename = "../tracks/MS_Plus_disk3_OK_track.flux";

	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " <flux image>\n";
		std::cout << "Assuming " << flux_filename << "\n";
	} else {
		flux_filename = argv[1];
	}

	std::cout << "Analyzing flux file " << flux_filename << "\n";
	std::vector<flux_record> flux_records = 
		get_flux_record(flux_filename, true);

	// Just do something with the first record.
	dp_results res = do_dp(flux_records[0], 100000); // 1234567 causes OOM
	print_dp_results(res);
	test_dp_validity_consistent(res);

	return 0;
}
