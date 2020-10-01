#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>
#include <limits>

#include <math.h>

const int MAX_DELAY = 100;
const int BANDS = 3;
const int P = 1; // Lp norm penalization. Must be >= 1.

// Since floating point addition can have numerical imprecision problems,
// I use a sort of fixed point by having integer penalties and multiplying
// the results from distance() by this.
const int FIXEDPOINT = 1000;

// TODO: Fix off-by-one problem that leaves the most recent bit
// unclassified.


// All the MFM decoder cares about at time t is:
//		What are the upper bounds to categorize the flux delay
//			into a given band
//		Did the last flux delay produce a surplus (residual) bit, e.g
//			001 becoming 00 and then 1 to be prefixed to the next
//			flux delay?
//		Did this flux delay decode to something that produced a zero
//			immediately before the residual bit? This determines whether
//			the valid encoding for 0 is 10 (if the previous bit was 0) or
//			00 (if the previous bit was 1).

// In addition, we keep track of the penalty so far and the number of bits
// decoded, even though that isn't strictly a floppy state; as well as a
// pointer to the best previous state.

// ....

// Problem: classifying everything into band 1 will give a penalty of zero
// and no errors (try it). Obviously, I don't want that. (Even if I set the
// beginning to be fixed, it'll just gradually grow the 1 band after that
// constraint. But perhaps doing so will be too costly and so force it to
// behave. See test later on (marked TEST).)

// TODO: Force the first k band classifications to be of a particular type.
// When forcing, it should also ignore any valid/not valid distinctions.

// Also TODO: The triangle inequality trick. Should make it *very* fast.

class floppy_classifier_state {
	public:
		int upper_bounds[BANDS-1];
		bool has_residual_bit;
		bool decoded_a_zero;

		unsigned int penalty_so_far;
		unsigned int bits_decoded;

		const floppy_classifier_state * prev;

		floppy_classifier_state() {
			// init the array
			for (int i = 0; i < BANDS-2; ++i) { upper_bounds[0] = 0; }
			upper_bounds[BANDS-2] = MAX_DELAY-1;

			// Whether a zero was decoded or not has to be decided by the
			// decoder, so the default value is undefined. Let it be false.
			decoded_a_zero = false;

			// An empty stream must start with a residual one. Consider a
			// stream with only one flux delay, say 20. That means that
			// there's 20 time units since last flux reversal, which means
			// there must have been a flux reversal before the one
			// the flux delay number gives the time to. (Fencepost logic.)
			has_residual_bit = true;

			// This way all the "incorrect" floppy_classifier_state get a
			// penalty of infinity (or as close as we can get).
			penalty_so_far = std::numeric_limits<unsigned int>::max();

			bits_decoded = 0;

			prev = NULL;
		}

};

typedef
	std::vector<		// t-th flux transition
	std::vector<		// upper bound[0]
	std::vector<		// upper bound[1]
	std::vector<		// decoded a zero
	std::vector<		// has residual bit after decoding
	floppy_classifier_state
>>>>> penalty_arr;

double ipow(int base, int exponent) {
	if (exponent == 1) {
		return abs(base);
	}

	if (exponent % 2 == 0) {
		double x = ipow(base, exponent/2);
		return x * x;
	} else {
		return base * ipow(base, exponent-1);
	}
}

// This function is called tons of times, thus the optimization.
double distance(const int a[BANDS-1], const int b[BANDS-1], int p) {

	double sum = 0;
	for (size_t i = 0; i < BANDS-1; ++i) {
		sum += ipow((a[i] - b[i]), p);
	}

	return pow(sum, 1.0/p);
}

double distance(const floppy_classifier_state & a,
	const floppy_classifier_state & b, int p) {

	return distance(a.upper_bounds, b.upper_bounds, p);
}

unsigned int distance_fixed(const floppy_classifier_state & a,
	const floppy_classifier_state & b) {

	return distance(a, b, P) * FIXEDPOINT;
}


// Hold on a sec, we just need past residual, present band, past had zero,
// and we want present residual, present had zero. That can be done by a
// simple array lookup and precomputation! So I wasn't off the mark before.
// Do the precomputation here.

class present_data {
	public:
		bool has_residual_bit;
		bool decoded_a_zero;
		bool valid;
		int bits_decoded;

		present_data() {
			has_residual_bit = false;
			decoded_a_zero = false;
			valid = true;
			bits_decoded = 0;
		}
};

// TODO: Use this unmemoized. Then once we have the thing running,
// memoize it.
present_data get_present_data(bool past_has_residual_bit, bool past_had_zero,
	int present_band) {

	// Manually clock the bits. We can afford to be slow because this is
	// going into a precomputed array anyway.

	std::vector<bool> flux_bits;

	if (past_has_residual_bit) {
		flux_bits.push_back(true);
	}

	// Decode current flux band: as many zeroes as in the band
	// (one-indexed), then a one. Since our input is zero-indexed, it
	// must provide an additional zero, hence the <=.
	for (int i = 0; i <= present_band; ++i) {
		flux_bits.push_back(false);
	}
	flux_bits.push_back(true);

	present_data out;

	// MFM decodes groups of two bits each; if there's any trailing
	// bit, set the output info structure to reflect it.
	if (flux_bits.size() % 2 != 0) {
		out.has_residual_bit = true;
		flux_bits.resize(flux_bits.size()-1);
	}

	// Start by assuming the present comes just after the past, which
	// means the "decoded a zero from the previous pair" part must be
	// set to what the past one has.
	out.decoded_a_zero = past_had_zero;
	out.valid = true; // assume the best.
	for (int i = 0; i < flux_bits.size(); i += 2) {
		// Do an MFM decode!
		bool first = flux_bits[i], sec = flux_bits[i+1];

		if (!first && !sec) {	// 0 when preceded by 1, otherwise invalid
			if (out.decoded_a_zero) {
				out.valid = false;
				return out;
			}
			out.decoded_a_zero = true; // just decoded a zero
		}

		if (!first && sec) {	// 1
			out.decoded_a_zero = false; // just decoded a one
		}

		if (first && !sec) {	// 0 when preceded by 0, otherwise invalid
			if (!out.decoded_a_zero) {
				out.valid = false;
				return out;
			}
			out.decoded_a_zero = true;
		}
		++out.bits_decoded;
	}

	return out;
}

void grow_penalties_array(int t, const floppy_classifier_state & reference,
	penalty_arr & penalties) {
	// If the penalties array is too small, update it.
	if (penalties.size() < t) {
		penalties.resize(t);
	}
	if (penalties[t].size() < MAX_DELAY) {
		penalties[t].resize(MAX_DELAY);
	}
	if (penalties[t][reference.upper_bounds[0]].size() < MAX_DELAY) {
		penalties[t][reference.upper_bounds[0]].resize(MAX_DELAY);
	}
	if (penalties[t][reference.upper_bounds[0]][reference.upper_bounds[1]].size() < 2) {
		penalties[t][reference.upper_bounds[0]][reference.upper_bounds[1]].resize(2);
	}
	if (penalties[t][reference.upper_bounds[0]][reference.upper_bounds[1]][reference.decoded_a_zero].size() < 2) {
		penalties[t][reference.upper_bounds[0]][reference.upper_bounds[1]][reference.decoded_a_zero].resize(2);
	}
}

int get_band(int flux_delay, const floppy_classifier_state & state) {
	int band = 0;			// zero-indexed, i.e. first band is zero
	for (int i = 0; i < BANDS-1; ++i) {
		if (flux_delay >= state.upper_bounds[i]) {
			band = i + 1;
		}
	}
	return band;
}

// Display some info about the solution that has been found.
// Returns the bands the stream decodes into.
std::vector<int> backtrack(const floppy_classifier_state & last,
	const std::vector<int> & pulses) {

	std::vector<const floppy_classifier_state *> chain;
	std::vector<int> bands;

	// Push the whole sequence onto the list.
	chain.push_back(&last);
	while ((*chain.rbegin())->prev != NULL) {
		chain.push_back((*chain.rbegin())->prev);
	}

	std::cout << "CS:" << chain.size() << std::endl;

	// Reverse so we look at the first bit first.
	std::reverse(chain.begin(), chain.end());

	int t = 0;
	for (const floppy_classifier_state * cur: chain) {
		int pulse = pulses[t++];
		std::cout << "Pulse " << pulse << " (band " << get_band(pulse, *cur) << ")\t";
		std::cout << "Upper bounds: " << cur->upper_bounds[0] << ", " <<
			cur->upper_bounds[1] << "\t" << "penalty: " <<
			cur->penalty_so_far << std::endl;
		bands.push_back(get_band(pulse, *cur));
	}

	return bands;
}

// Update the present floppy classifier state to point at the valid past
// floppy classifier that minimizes the total penalty in p norm. We use
// the same trick as below to recurse through every possible candidate,
// but going from closest to furthest away to take advantage of a triangle
// inequality trick.
// note that t represents the timestep of *present*, not past, so we will
// be accessing penalties[t-1] to find the optimum past.
void update_optimal_past(int i, int t, bool past_has_residual_bit,
	bool past_had_zero, int bits_decoded_in_present,
	floppy_classifier_state & present,
	floppy_classifier_state & past, penalty_arr & penalties) {

	if (t == 0) {
		// There must be a 1 before anything happens.
		if (!past_has_residual_bit) {
			return;
		}
		// Everything else is allowed (assuming prior zero or not assuming
		// prior zero, everything goes).
		// TEST: Prior assumption of where the bounds reside.
		present.penalty_so_far = 0;
		past.penalty_so_far = 0;
		/*floppy_classifier_state test;
		test.upper_bounds[0] = 30;
		test.upper_bounds[1] = 40;
		present.penalty_so_far = distance_fixed(test, present);*/
		return;
	}

	// Set past upper bound
	// but consider the early pruning failure that might happen!
	// No, I don't think that would be a problem.
	// Pruning failure: our triangle inequality heuristic uses the fact
	// that if a point closer to the present is valid, then it must
	// have a lesser penalty than any point further away. But suppose that
	// the "sorted order" invariant forces some of the points to change.
	// E.g. we start at ub[0] = 50 and that forces ub[1] to 50, but the
	// point (43, 49) would've been closer.
	// But that's not a problem because we start with past being equal to
	// present (because that has distance 0); so no such clamping can happen
	// when going from present to past.
	if (i < 2) {
		// The upper bound values must be in sorted order.
		int min_value = 0, cur_delta = 0;
		if (i > 0) { min_value = past.upper_bounds[i-1]; }

		for (int delta = 0; delta < MAX_DELAY; ++delta) {
			// Walk towards 0 from the present's value
			int prospective_ub = present.upper_bounds[i] - delta;
			if (prospective_ub >= min_value) {
				past.upper_bounds[i] = prospective_ub;
				update_optimal_past(i+1, t, past_has_residual_bit,
					past_had_zero, bits_decoded_in_present, present, past,
					penalties);
			}

			// Walk towards max
			prospective_ub = present.upper_bounds[i] + delta;
			if (prospective_ub >= min_value && prospective_ub < MAX_DELAY) {
				past.upper_bounds[i] = prospective_ub;
				update_optimal_past(i+1, t, past_has_residual_bit,
					past_had_zero, bits_decoded_in_present, present, past,
					penalties);
			}

			// Cleanup
			past.upper_bounds[i] = present.upper_bounds[i];
		}
		return;
	}

	// Redimension penalties array in case we're accessing something
	// inaccesible.
	grow_penalties_array(t-1, past, penalties);

	// If we get here, all the state for the past has been set.
	past = penalties[t-1][past.upper_bounds[0]][past.upper_bounds[1]]
		[past.decoded_a_zero][past.has_residual_bit];

	unsigned int past_present_dist,
		past_present_bitcount = past.bits_decoded + bits_decoded_in_present;

	// Avoid overflow
	if (past.penalty_so_far == std::numeric_limits<unsigned int>::max()) {
		past_present_dist = past.penalty_so_far;
	} else {
		past_present_dist = past.penalty_so_far +
			distance_fixed(past, present);
	}

	/*std::cout << "\t" << past_present_dist << " vs " << present.penalty_so_far << std::endl;
	std::cout << "Past UB: " << past.upper_bounds[0] << ", " << past.upper_bounds[1] << std::endl;
	std::cout << "Present UB: " << present.upper_bounds[0] << ", " << present.upper_bounds[1] << std::endl;

	if (past_present_dist < present.penalty_so_far) {
		std::cout << "HM" << std::endl;
	}*/

	// If choosing this as the past beats the record, update the record.
	if (past_present_dist < present.penalty_so_far ||
		(past_present_dist == present.penalty_so_far &&
			past_present_bitcount > present.bits_decoded &&
			past_present_dist != std::numeric_limits<unsigned int>::max())) {

		present.penalty_so_far = past_present_dist;
		present.bits_decoded = past_present_bitcount;
		present.prev = &penalties[t-1][past.upper_bounds[0]]
			[past.upper_bounds[1]][past.decoded_a_zero]
			[past.has_residual_bit];
		//std::cout << "OK" << std::endl;
	}
}

// The challenger is better than the record if its penalty is lower, or if
// the penalty is the same but the challenger decodes more bits.
bool challenger_better_than_record(const floppy_classifier_state & challenger,
	const floppy_classifier_state & old_record) {

	return challenger.penalty_so_far < old_record.penalty_so_far ||
		(challenger.penalty_so_far == old_record.penalty_so_far &&
			challenger.bits_decoded > old_record.bits_decoded);
}

// i is the state index (upper bound 0 and 1, has residual bit,
// decoded a zero). t is the flux delay index that we're generating all
// possible options for (i.e. t-1 is the last valid slice of the
// penalties array before entering walk_state).
// record is the best (lowest penalty) solution of them all, and
// forced_bands is a list of bands to force the first few values to attain.
// (E.g. the A1A1A1 preamble.) As forced bands may be in error (e.g. the OOB
// markers in A1A1A1), validity checks are disabled as long as the band is
// forced.
void walk_state(int i, int t, int num_bits_desired,
	const std::vector<int> & pulses, floppy_classifier_state & present,
	penalty_arr & penalties, floppy_classifier_state & record,
	const std::vector<int> & forced_bands) {

	// Set upper bound
	if (i < 2) {
		// The upper bound values must be in sorted order.
		int min_value = 0;
		if (i > 0) { min_value = present.upper_bounds[i-1]; }

		for (present.upper_bounds[i] = min_value;
			present.upper_bounds[i] < MAX_DELAY; ++present.upper_bounds[i]) {
			walk_state(i+1, t, num_bits_desired, pulses, present,
				penalties, record, forced_bands);
		}
		return;
	}

	// Decode the flux transition to determine whether we're odd or even.
	// If odd: a prior residual bit means no residual bit this time. If
	// even: a prior residual bit means residual bit this time.

	int band = get_band(pulses[t], present); // zero-indexed.

	// If we've forced the band, and the band we got doesn't match the
	// forced band, exit.
	if (t < forced_bands.size() && band != forced_bands[t]) {
		//std::cout << "Rejected " << present.upper_bounds[0] << ", " << present.upper_bounds[1] << " because force." << std::endl;
	 return; }

	// Do this the ugly way. Assume the past is a certain way,
	// and then determine if the present, as defined above, is consistent
	// with it.

	for (bool past_has_residual_bit: {false, true}) {
		for (bool past_had_zero: {false, true}) {
			present_data present_settings = get_present_data(
				past_has_residual_bit, past_had_zero, band);

			// No use continuing if that settings combination wasn't
			// valid -- unless we're forced.
			if (!present_settings.valid && t >= forced_bands.size()) {
				//std::cout << "Rejected " << present.upper_bounds[0] << ", " << present.upper_bounds[1] << " because invalid. " << std::endl;
				continue;
			}

			// if (t >= forced_bands.size()) { std::cout << "Out of forced domain. " << std::endl; }

			// Set the present to the settings provided...
			present.has_residual_bit = present_settings.has_residual_bit;
			present.decoded_a_zero = present_settings.decoded_a_zero;

			// Bits decoded is a cumulative count and thus depends on the
			// past (which we don't know). Penalty ditto, so set them both
			// to the worst possible values so they'll be overwritten.
			present.bits_decoded = 0;
			present.penalty_so_far = std::numeric_limits<unsigned int>::max();

			// And go looking for an optimal past.
			floppy_classifier_state past = present;
			// Perhaps refactor this so the for loop directly uses these?
			past.has_residual_bit = past_has_residual_bit;
			past.decoded_a_zero = past_had_zero;

			update_optimal_past(0, t, past_has_residual_bit, past_had_zero,
				present_settings.bits_decoded, present, past, penalties);

			// We now have the optimal penalty for present under the
			// constraints given. Update if it's better than what we have.

			// If the penalties array is too small, update it.
			grow_penalties_array(t, present, penalties);

			floppy_classifier_state * old_record = &penalties[t][present.upper_bounds[0]][present.upper_bounds[1]]
				[present.decoded_a_zero][present.has_residual_bit];

			// Update the record if what we have is better.
			if (challenger_better_than_record(present, *old_record)) {
				*old_record = present;
			}

			// std::cout << "with past: " << past.penalty_so_far << std::endl;
		}
	}

	// std::cout << "finding present: " << t << ", " << present.upper_bounds[0] << ", " << present.upper_bounds[1] << ": " << present.penalty_so_far << std::endl;

	// If we've found a decoding that gets the number of bits that's
	// desired, output it and (sloppily) exclude it from consideration.
	if (present.bits_decoded >= num_bits_desired) {
		std::cout << t << ": " << present.penalty_so_far << "\t" <<
			present.bits_decoded << std::endl;
		//backtrack(present, pulses);
		if (challenger_better_than_record(present, record)) {
			std::cout << t << std::endl;
			record = present;
		}
		present.penalty_so_far = std::numeric_limits<unsigned int>::max();
	}
}

floppy_classifier_state walk_state(int t, int num_bits_desired,
	const std::vector<int> & pulses, penalty_arr & penalties,
	const std::vector<int> & forced) {

	floppy_classifier_state present, old_record;

	walk_state(0, t, num_bits_desired, pulses, present, penalties,
		old_record, forced);

	return old_record;
}

floppy_classifier_state walk_state(int t, int num_bits_desired,
	const std::vector<int> & pulses, penalty_arr & penalties) {

	return walk_state(t, num_bits_desired, pulses, penalties);
}

// A valid IDAM.
std::vector<int> get_test_pulses() {
	return {
		23, 22, 23, 22, 22, 22, 23, 22, 23, 22, 22, 22, 23, 22, 22, 22, 23,
		22, 22, 22, 22, 22, 23, 22, 23, 22, 23, 22, 23, 22, 23, 22, 23, 22,
		22, 22, 23, 22, 23, 22, 23, 22, 23, 22, 23, 22, 23, 22, 23, 22, 23,
		22, 23, 22, 23, 34, 47, 34, 47, 34, 22, 46, 36, 45, 35, 21, 47, 34,
		47, 34, 23, 22, 23, 22, 23, 22, 23, 34, 23, 22, 22, 22, 35, 34, 23,
		22, 23, 22, 22, 22, 23, 34, 35, 22, 23, 34, 47, 21, 35, 22, 23, 22,
		23, 34, 47, 22, 47, 21, 35, 22, 22, 35, 22, 34, 35, 46, 47, 34};
}

// The right bands for an A1A1A1 preamble with zeroes in front.
// (Zero-indexed.)
std::vector<int> get_preamble_assignment() {
	return {
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1, 2, 1, 0, 2, 1, 2, 1, 0,
		2, 1, 2, 1};
}

// What the IDAM should decode to if everything goes well.
std::vector<int> desired_band_assignment() {
	return {
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1, 2, 1, 0, 2, 1, 2, 1, 0,
		2, 1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0,
		0, 0, 0, 1, 1, 0, 0, 1, 2, 0, 1, 0, 0, 0, 0, 1, 2, 0, 2, 0, 1, 0,
		0, 1, 0, 1, 1, 2, 2, 1};
}

std::vector<int> mult_warp(const std::vector<int> & pulses,
	double amplitude, double freq) {

	std::vector<int> warped;

	for (int i = 0; i < pulses.size(); ++i) {
		// Scale a sine from -1..1 to (1-amplitude)...(1+amplitude)
		double warp_factor = (sin(i * freq) + 1)/2; // now 0..1
		warp_factor = 1 + 2 * amplitude - amplitude; // now 1-amp..1+amp

		warped.push_back(pulses[i] * warp_factor);
	}
	return warped;
}

std::vector<int> additive_warp(const std::vector<int> & pulses,
	double minimum, double maximum, double freq) {

	std::vector<int> warped;

	for (int i = 0; i < pulses.size(); ++i) {
		// Scale a sine from -1..1 to
		double offset = (sin(i * freq) + 1)/2;
		offset = offset * (maximum - minimum) + minimum;

		warped.push_back(pulses[i] + offset);
	}
	return warped;
}

int main() {
	/*int T = 75; // e.g.
	int num_bits_desired = 100; // e.g.
	std::vector<int> pulses; // pulses, e.g
	for (int i = 0; i < T; ++i) {
		pulses.push_back(random() % MAX_DELAY);
	}
	*/

	std::vector<int> pulses = get_test_pulses(), forced = get_preamble_assignment();
	int T = pulses.size();
	int num_bits_desired = 130;

	// Let's make it a bit more challenging
	// And it breaks due to the draw towards penalty=0
	pulses = additive_warp(pulses, 33, 50, 0.05);

	std::copy(pulses.begin(), pulses.end(), std::ostream_iterator<int>(std::cout, " "));
	std::cout << std::endl;
	penalty_arr penalties;

	penalties.resize(T);

	floppy_classifier_state best;

	for (int t = 0; t < T; ++t) {
		std::cout << t << std::endl;
		floppy_classifier_state best_this_t = walk_state(t, num_bits_desired,
			pulses, penalties, forced);

		if (challenger_better_than_record(best_this_t, best)) {
			best = best_this_t;
		}
	}

	std::cout << "\nLeast penalty: " << std::endl;
	std::vector<int> observed = backtrack(best, pulses);
	std::vector<int> expected = desired_band_assignment();

	std::cout << "Num pulses: " << pulses.size() << std::endl;
	std::cout << "Observed " << observed.size() << " assignments, expected " << expected.size() << std::endl;

	std::vector<int> expected_subset = expected;
	expected_subset.resize(observed.size());

	if (expected_subset == observed) {
		std::cout << "All observed bands are correct." << std::endl;
		if(expected.size() != observed.size()) {
			std::cout << "But the number differs!" << std::endl;
		}
	} else {
		std::cout << "Discrepancy between observed and expected!" << std::endl;
	}
}