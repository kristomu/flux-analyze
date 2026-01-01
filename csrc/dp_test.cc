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
#include <iostream>
#include <limits>
#include <cmath>

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

// TODO: Decide what the parity is from the perspective of - the end of LPZB
// or the end of it plus our own zero bits?

// Checks if the sequence 1, a zeroes, 1, b zeroes, breaks any rules, when
// the sequence 1, a zeroes has already been tested and can be assumed to
// not have broken any rules. This is done by padding with a one on the left
// and a zero on the right if required. It can probably be offloaded into a
// lookup table if it proves to be too slow.
bool is_mfm_valid(short zero_bits, u_char parity) {
	return parity == 0 || zero_bits != 3;
}

enum f_mode_t {
	IN_BAND = 0,
	PREAMBLE_A = 1,
	PREAMBLE_C = 2,
	ILLEGIBLE = 3,
	NUM_MODES = 4};

// perhaps some better class names?
class dp_idx {
	public:
		f_mode_t mode;
		short clock;
		int idx;
		u_char parity;
};

class dp_record {
	public:
		short zero_bits = -1;
		short OOB_sequence_index = -1;
		double penalty = 0;
		dp_idx previous_opt;
};

const int MIN_CLOCK = 1, MAX_CLOCK = 100;

double sqr(double x) { return x*x; }

std::vector<int> do_dp(const flux_record & f, size_t len) {

	// Noise sensitivity hyperparameter.
	// Higher values of this allows greater changes in clock
	// from one pulse to the next, while lower values prioritize fitting
	// the set clock to the current pulse delay.
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
	std::vector<std::vector<dp_idx> > possible_states_by_parity(2);
	dp_idx current;

	for (int cur_mode = 0; cur_mode < (int)NUM_MODES; ++cur_mode) {
		current.mode = (f_mode_t) cur_mode;
		for (current.clock = MIN_CLOCK; current.clock <= MAX_CLOCK;
			++current.clock) {
			for (current.parity = 0; current.parity < 2; ++current.parity) {
				possible_states.push_back(current);
				possible_states_by_parity[current.parity].push_back(current);
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

			dp_record cur_record;

			size_t observed_pulse_delay = fluxes[i];
			cur_record.zero_bits = round(
				2.0 * observed_pulse_delay / current_state.clock) - 1;
			// Pulse lengths must be 1, 2, or 3.
			cur_record.zero_bits = std::min(3, std::max(1,
				(int)cur_record.zero_bits));

			double expected_pulse_delay =
				current_state.clock * (cur_record.zero_bits + 1)/ 2.0;

			// If we're not the first flux delay, determine the least-penalty
			// previous state to use.

			if (i == 0) {
				cur_record.penalty =
					alpha * sqr(observed_pulse_delay - expected_pulse_delay);

				dp[current_state.mode][current_state.clock]
					[current_state.idx][current_state.parity] = cur_record;
				continue;
			}

			// Our current state has a particular parity, and the number of
			// zero bits forces all our prior states to be of a particular
			// parity. Determine which.

			// all zeroes plus trailing edge
			int bits_generated = cur_record.zero_bits + 1;
			u_char actual_parity = bits_generated % 2;

			// If actual parity is the same as the current state's claimed
			// parity, then the previous state must be even, otherwise odd
			// (because 0 parity is even).
			u_char required_past_parity = 0;
			if (actual_parity != current_state.parity) {
				required_past_parity = 1;
			}

			bool primed = false;
			double best_penalty;
			dp_idx best_idx;

			for (dp_idx last_state: possible_states_by_parity[required_past_parity]) {

				last_state.idx = i-1;
				dp_record last_record = dp[last_state.mode][last_state.clock]
					[last_state.idx][last_state.parity];

				// Calculate the penalty when transitioning from this
				// state. (I may move this to another function later.)

				double candidate_penalty = sqr(current_state.clock - last_state.clock) +
					alpha * sqr(observed_pulse_delay - expected_pulse_delay);

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

				dp[current_state.mode][current_state.clock]
					[current_state.idx][current_state.parity] = cur_record;
			}
		}
	}

	// ============= Get the solution. ========

	// First determine the least penalty solution.

	double best_penalty = std::numeric_limits<double>::infinity();
	dp_idx recordholder;

	for (dp_idx current_state: possible_states) {
		current_state.idx = fluxes.size()-1;

		double candidate_penalty = dp[current_state.mode]
			[current_state.clock][current_state.idx]
			[current_state.parity].penalty;

		if (candidate_penalty <= best_penalty) {
			recordholder = current_state;
			best_penalty = candidate_penalty;
		}
	}

	std::cout << "Optimal penalty: " << best_penalty << "\n";

	std::vector<dp_idx> optimal_path;
	std::vector<dp_record> optimal_values;

	// Trace the solution back through the DP array

	for (size_t i = 0; i < fluxes.size(); ++i) {
		optimal_path.push_back(recordholder);

		dp_record at_current = dp[recordholder.mode]
			[recordholder.clock][recordholder.idx]
			[recordholder.parity];

		optimal_values.push_back(at_current);

		recordholder = at_current.previous_opt;
	}

	// And reverse the order of that path to get beginning to end.

	std::reverse(optimal_path.begin(), optimal_path.end());
	std::reverse(optimal_values.begin(), optimal_values.end());

	std::vector<int> clocks, zero_bits;
	std::vector<double> penalties;
	std::vector<int> valid;

	for (size_t i = 0; i < optimal_values.size(); ++i) {
		clocks.push_back(optimal_path[i].clock);
		zero_bits.push_back(optimal_values[i].zero_bits);
		penalties.push_back(optimal_values[i].penalty);
		if (is_mfm_valid(optimal_values[i].zero_bits,
			optimal_path[i].parity)) {
			valid.push_back(1);
		} else {
			valid.push_back(0);
		}
	}

	// This shows a vast field of ones interspersed with
	// 1 0 1 1 1 1 0 1 1 1 1 0 1, which I think is the intentional MFM
	// errors of the OOB A1A1A1/C2C2C2 preambles.

	// But try the MS Plus disk3 OK track (track 4 side 0) with 10000
	// entries and it flips out. This suggests that my choice of the alpha
	// hyperparameter is not a good one - or the penalty function itself
	// isn't.

	std::copy(valid.begin(), valid.end(), std::ostream_iterator<int>(std::cout, " " ));
	std::cout << "\n";

	// TODO: Manually verify that the valid/invalid values are correct
	// by doing an old-fashioned MFM decode on the zero bits output and the
	// initial parity.

	return zero_bits;
}

int main(int argc, char ** argv) {

	std::string flux_filename = "tracks/MS_Plus_disk3_OK_track.flux";

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
	do_dp(flux_records[0], 1000); // 1234567 causes OOM

	return 0;
}
