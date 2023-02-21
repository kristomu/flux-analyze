// From https://pastebin.com/w0S40FV3

void pll(int delay) {
	// Store pulse delay in a variable.
	int old_delay = delay;

	// Calculate limits
	// This is definitely too ugly, fix later.
	// As far as I can see, these correspond to the limits
	// given by treating a delay as a clock and then rounding to
	// the nearest half-clock:
	//			2 half-clocks range from 1.5 to 2.5 half-clocks
	//				i.e. 3/2 to 5/2 half-clocks, i.e. 3/4 to 5/4 clocks.
	//			3 half-clocks range from 2.5 to 3.5 half-clocks, i.e.
	//				to 7/2 half-clocks, or 7/4 clocks.
	//			4 half-clocks range from 3.5 to 4.5 half-clocks, i.e.
	//				to 9/2 half-clocks, or 9/4 clocks. Alright.
	int one_impulse_low = 3 * old_delay/4;
	int one_impulse_nom = old_delay;
	int one_impulse_high = 5 * old_delay / 4;

	int two_impulse_low = one_impulse_high;
	int two_impulse_nom = 3 * old_delay / 2;
	int two_impulse_high = 7 * old_delay / 4;

	int three_impulse_low = two_impulse_high;
	int three_impulse_nom = 2 * old_delay;
	int three_impulse_high = 9 * old_delay / 4;

	int length_snap = length;
	int new_delay;

	// But these intervals are *weird*!
	// Everything below 1/2 half-clock is too small
	//	everything between 1/2 and 5/2 half-clocks count as 1 us
	//  everything between 5/2 and 7/2 count as a 1.5 us
	//  everything btween 7/2 and 12/2 count as 2
	//  and everything higher is out of scope. (????)
	// Maybe https://oeis.org/A280291 ???

	// The phase (delay) updates are also very weird.
	// But since the filter and update only updates based on a
	// pass/fail criterion, it would strictly speaking be possible
	// to propagate these over.

	// ...

	// I'll imagine that all these updates are really of the
	// form: alpha * old_delay + (1-alpha) * c_hat, where c_hat is what the
	// clock would be for the observed delay to be exactly nom.

	// E.g. for 1 us between pulses,
	// c_hat = length_snap, so
	// new_delay = 1/2 * old_delay + 1/2 * c_hat.
	// For 1.5 us between two pulses
	// length_snap = 3/2 c_hat, so
	// new_delay = 16/32 old_delay + 11/32 * 3/2 c_hat =
	//	32/64 old_delay + 33/64 c_hat		almost!
	// and for 2 us between pulses
	// length_snap = 2 * c_hat, so
	// new_delay = (2 * old_delay + 2 * c_hat) / 4
	
	// So yep. With the exception of the 33/64 part, which is off by 1/64,
	// it's all
	// new_delay = old_delay/2 + c_hat/2.
	// I guess I could handwave away the 1/64 by saying clocks
	// usually are below that and so it would be rounded off anyway.

	// ... which helps explain why my sliding window solution in Python
	// is so good.

	// Another clue to get from this is to not update on very long
	// intervals. From one POV, it makes sense because we don't know
	// how the interval was divided up - e.g. is it four intervals of
	// two half-clocks or two of four? On the other, we could just go
	// with a reasonable approximation.

	if (length_snap < delay/4) {
		// too small
	} else if (length_snap < one_imp_high) {
		// 1 us between two pulses
		new_delay = (old_delay + length_snap)/2; // ???
	} else if (length_snap < two_imp_high) {
		// 1.5 us between two pulses
		new_delay = (16 * old_delay + 11 * length_snap) / 32;
	} else if (length_snap < 3 * delay) {
		// 2 us between two pulses
		new_delay = (2 * old_delay + length_snap) / 4;
	} else {
		// too long
	}

	// I like this robust approach: if the new clock estimate is high,
	// then multiply by a fixed amount; if it's low, then divide, and 
	// if it's neither then do nothing. I might crib it.

	// This is still a two-hyperparameter solution, though: the alpha
	// (above), and the increment/decrement step below.

	// Filter and update delay.
	if (new_delay > 9 * old_delay / 8) {
		delay = 9 * old_delay / 8;
	} else if (new_delay < 7 * old_delay / 8) {
		delay = 7 * old_delay / 8;
	}
}