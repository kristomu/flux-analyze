// Experimenting with dynamic programming as a way to recover the data stream.

/*	The overarching idea is that any data on a floppy may fall within one of three
	categories:

		- In-band data (obeys MFM state transition rules)
		- Out-of-band data (A1A1A1 or C2C2C2 preambles)
		- Illegible data (wrong format, corrupted, etc).

	The dynamic programming approach would read the whole flux stream for every
	possible clock at once, slowly varying the allowed clock according to a penalty
	function that punishes the pulse delay from being far from what's expected and
	for the current clock to change too much from the previous pulse's clock.

	This, as such, isn't too different from just doing a moving or exponential
	average. But I've then intended to make the DP system enforce that all in-band
	data must obey MFM rules. This would then make it switch modes exactly where
	the out-of-band data appears, while trying to recover legible MFM data where
	severe warping/clock skew occurs.

	Then I might add support for missing bits, where the system can add split up
	pulses where doing so is beneficial. In combination with known-length constraints
	(e.g. that a sector that metadata says should have 4096 bits has exactly 4096
	bits), this should make the process even better at recovering corrupted data.

	However, I think braiding all of this into a DP algorithm could get very messy.
	So beware, all who read this. I may need multiple iterations to get it to be
	somewhat presentable.
*/

/*	The current dynamic programming idea is as follows, where * indicates something
	that hasn't been implemented yet.

	There are four modes:
		  0 - in-band
		* 1 - out-of-band, A1A1A1 preamble
		* 2 - out-of-band, C2C2C2 preamble
		* 3 - everything else that doesn't fit the model ("unknown"/"illegible")

	The state for each step is:
		INDEX:	Current mode
		INDEX:	Current estimated clock value
		INDEX:	Current index into list of flux transitions (ith transition)
		INDEX:	MFM decoding parity				*

		STATE:	Current flux delay period (number of zero bits)
		STATE:	OOB sequence bit index			*
		STATE:	accumulated penalty/cost
		STATE:	previous optimal entry

	where INDEX is a dimension and STATE is a value at a concrete array cell.

	(When/if I add missing bit support, the current index INDEX should probably
	be "number of flux transitions processed", while the actual array index is
	(i, value left) so that we can insert flux transitions without disrupting
	which real flux transition is being checked.)

	MFM decoding parity is pretty tricky and needs a more detailed explanation.
	I'll give that explanation later.
*/

/*	For the first three modes, we're assuming the floppy is formatted with a
	standard MFM encoding, where 0 indicates a lack of flux reversal and 1 is
	its presence:

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

enum f_mode_t {
	IN_BAND = 0,
	PREAMBLE_A = 1,
	PREAMBLE_C = 2,
	ILLEGIBLE = 3,
	NUM_MODES = 4};

// perhaps some better class names?

// Inclusive parity is the parity of the sequence, i.e. number of
// bits emitted *including* the zeroes and trailing one due to the
// current pulse delay. Exclusive parity is the parity of the
// sequence up to, but not including, the current pulse delay.
class dp_idx {
	public:
		f_mode_t mode;
		short clock;
		int idx;
		u_char inclusive_parity;
};

class dp_record {
	public:
		short zero_bits = -1;
		short OOB_sequence_index = -1;
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
		std::vector<int> valid;
		u_char starting_ex_parity; // starting exclusive parity
};

// Minimum and maximum clock index value.
const int MIN_CLOCK = 1, MAX_CLOCK = 100;

// How much increasing index by one increases the clock,
// e.g. CLOCK_GRANULARITY = 0.5 would make clock index 50
// represent a clock value of 50 * 0.5 = 25.
const double CLOCK_GRANULARITY = 1;

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

	// We can probably cut down on the space complexity here by having
	// only two time instances (old and new) and recording the whole
	// trace... maybe do that later if necessary.
	// Or by doing it in chunks, e.g. 1000 a pop; get the trace from the
	// last 1000; then do the previous 1000 and stitch together from the
	// trace we got... O(kn) time, k space instead of O(n) time and space,
	// where k is the chunk length.
	// We'll see.

	std::vector<std::vector<std::vector<std::vector<dp_record> > > > dp(
		NUM_MODES, std::vector<std::vector<std::vector<dp_record> > >(
			MAX_CLOCK+1, std::vector<std::vector<dp_record> >(
				len, std::vector<dp_record>(2))));

	// Get every possible index state so we can iterate over them
	// later - this saves a lot of nested for loops. Also keep
	// a track of possible states by parity, since valid past states
	// are limited to having the same parity.
	std::vector<dp_idx> possible_states;
	std::vector<std::vector<dp_idx> > possible_states_by_inc_parity(2);
	dp_idx current;

	for (int cur_mode = 0; cur_mode < (int)NUM_MODES; ++cur_mode) {
		current.mode = (f_mode_t) cur_mode;
		for (current.clock = MIN_CLOCK; current.clock <= MAX_CLOCK;
			++current.clock) {
			for (current.inclusive_parity = 0; current.inclusive_parity < 2;
				++current.inclusive_parity) {

				possible_states.push_back(current);
				possible_states_by_inc_parity[
					current.inclusive_parity].push_back(current);
			}
		}
	}

	for (size_t i = 0; i < fluxes.size(); ++i) {
		for (dp_idx current_state: possible_states) {
			current_state.idx = i;

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
					[current_state.idx][current_state.parity] = cur_record;
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

			for (dp_idx last_state: possible_states_by_inc_parity[exclusive_parity]) {

				last_state.idx = i-1;
				dp_record last_record = dp[last_state.mode][last_state.clock]
					[last_state.idx][last_state.inclusive_parity];

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

				/*
				if (!is_mfm_valid(cur_record.zero_bits,
					exclusive_parity)) {
					candidate_penalty += 1;
				}
				*/

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

				dp[current_state.mode][current_state.clock][current_state.idx]
					[current_state.inclusive_parity] = cur_record;
			}
		}
	}

	// ============= Get the solution. ========

	// First determine the least penalty solution.

	double best_penalty = std::numeric_limits<double>::infinity();
	dp_idx best_idx;

	for (dp_idx current_state: possible_states) {
		current_state.idx = fluxes.size()-1;

		double candidate_penalty = dp[current_state.mode]
			[current_state.clock][current_state.idx]
			[current_state.inclusive_parity].penalty;

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

		dp_record at_current = dp[best_idx.mode]
			[best_idx.clock][best_idx.idx]
			[best_idx.inclusive_parity];

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
	dp_results res = do_dp(flux_records[0], 1000); // 1234567 causes OOM
	print_dp_results(res);
	test_dp_validity_consistent(res);

	return 0;
}
