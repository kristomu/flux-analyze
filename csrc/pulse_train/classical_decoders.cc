#include <algorithm>
#include <limits>
#include <vector>

#include <iostream>

#include "classical_decoders.h"
#include "tools.h"

// Turn a flux record vector into an MFM pulse train by comparing the distance
// between flux reversals. Our job is made more difficult because apparently
// the mapping is nonlinear: one clock length corresponds to 01, 1.5 clock
// lengths to 001, and 2 clock lengths to 0001.

// For this reason, longer delays are undefined. I may handle these later,
// but note that they are out of spec and should never occur on uncorrupted
// floppies with the right clock.

// The error_out double is set to a badness-of-fit value: the higher the worse.
// It is used as a quick and dirty clock inference method (see below), but
// is more useful for debugging (e.g. when I implement dewarping, it could
// show how well the signal has been dewarped).

// In addition, the nonlinearity might save us sometimes: suppose that we have an
// interval of 2.5 clocks. Then it's likely that this is either one clock followed
// by 1.5 clocks (01001) or the other way around (00101), with some central flux
// reversal having been erased. But I'm not going to do recovery before I've got
// this working on normal images.

// TODO: It might be useful to create an "explicit clock vector decoder" that
// takes a list of clock values and creates an MFM_train_data object based on
// it. Then constant-clock is trivial and the "modern" decoders are fairly
// easy too. Only problem is that it's creating and destroying a bunch of
// temporary vectors when not doing that would be faster. Optimize later,
// perhaps.

// Level one:

MFM_train_data constant_clock_decoder::get_MFM_train(
		const std::vector<int> & fluxes, size_t start_pos,
		size_t end_pos, double & error_out) const {

	MFM_train_data train;

	train.data.push_back(1);
	train.flux_indices.push_back(start_pos);

	error_out = 0;

	for (size_t i = start_pos; i < end_pos; ++i) {
		// We'll model the nonlinearity like this:
		//		- There's always half a clock's delay before anything happens.
		//		- Then one zero corresponds to half a clock more, two zeroes is
		//			two halves more, and three zeroes is three, and so on.

		int flux_delay = fluxes[i];
		int half_clocks = round((flux_delay*2)/clock);

		// Subtract the constant half-clock offset of one, then round the next
		// to get the number of zeroes.
		int zeroes = std::max(0, half_clocks-1);

		// Update Euclidean error.
		double error_term = flux_delay - half_clocks * clock/2.0;
		error_out += error_term * error_term;

		// Ignore RR -- treat it as noise (as the reference VHDL code
		// does). Better might be to just add the delay to the next
		// term. Maybe do that later.
		if (zeroes == 0) { continue; }

		for (int j = 0; j < zeroes; ++j) {
			train.data.push_back(0);
			train.flux_indices.push_back(i);
		}
		train.data.push_back(1);
		train.flux_indices.push_back(i);
	}

	error_out = std::sqrt(error_out / (end_pos - start_pos));

	return train;
}

// Horribly cut and paste code, but it'll work as a proof of concept.
// MFM decoder with dewarping based on exponential smoothing.
// Better (and less cut and paste) stuff later.
// The exponential smoothing has no right being this powerful, either
// (not that I'm really complaining).

// If different friction or warping of the surface were to bias the
// clock, we'd expect the bias to be proportional to the length between
// clock transitions as a longer delay would mean more time to accumulate.
// The dewarp function below instead assumes that the bias is independent
// of length, and bizarrely, that seems to work *better*, at least for
// certain chunks that are otherwise very tough to decode. I don't get it,
// but I can't dispute the results.

MFM_train_data orig_EWMA_causal_clock_decoder::get_MFM_train(
	const std::vector<int> & fluxes, size_t start_pos,
	size_t end_pos, double & RMSE_out) const {

	MFM_train_data train;

	train.data.push_back(1);
	train.flux_indices.push_back(0);

	RMSE_out = 0;
	double bias = 0;

	double clock = initial_clock;

	for (size_t i = start_pos; i < end_pos; ++i) {
		int corrected_flux_delay = fluxes[i] - bias;
		int half_clocks = round((corrected_flux_delay*2)/clock);

		// The error term is the difference between the
		// observed (corrected) flux delay and what we would
		// have seen if the delay had been an exact multiple of
		// the initial clock rate.

		double error_term = corrected_flux_delay - half_clocks * clock/2.0;
		RMSE_out += error_term * error_term;

		bias = bias * (1-alpha) + error_term * alpha;

		int zeroes = std::max(0, half_clocks-1);

		if (zeroes == 0) { continue; }

		for (int j = 0; j < zeroes; ++j) {
			train.data.push_back(0);
			train.flux_indices.push_back(i);
		}
		train.data.push_back(1);
		train.flux_indices.push_back(i);
	}

	RMSE_out = std::sqrt(RMSE_out / (end_pos - start_pos));

	return train;
}

// Trying a reasonably faithful reproduction of the VHDL PLL code.
// The only difference is: the center for three half-clocks is at 3x
// a half clock, not 32/11 x.

// Maybe not faithful enough?
MFM_train_data historical_EWMA_decoder::get_MFM_train(
		const std::vector<int> & fluxes, size_t start_pos,
		size_t end_pos, double & RMSE_out) const {

	MFM_train_data train;
	
	train.data.push_back(1);
	train.flux_indices.push_back(0);

	RMSE_out = 0;
	double estimated_half_clock = initial_clock / 2;

	/*
	-- length limits
	variable one_imp_low  : integer range 0 to 1000 := 0;
	variable one_imp_nom  : integer range 0 to 1000 := 0;
	variable one_imp_high : integer range 0 to 1000 := 0;
	 
	variable two_imp_low  : integer range 0 to 1000 := 0;
	variable two_imp_nom  : integer range 0 to 1000 := 0;
	variable two_imp_high : integer range 0 to 1000 := 0;
	 
	variable tre_imp_low  : integer range 0 to 1000 := 0;
	variable tre_imp_nom  : integer range 0 to 1000 := 0;
	variable tre_imp_high : integer range 0 to 1000 := 0;
	*/

	int one_imp_low, one_imp_nom, one_imp_high;
	int two_imp_low, two_imp_nom, two_imp_high;
	int tre_imp_low, tre_imp_nom, tre_imp_high;
	 
	/*-- delay:
	variable old_delay    : integer range 0 to 1000 := 0;
	variable new_delay    : integer range 0 to 1000 := 0;
	 
	variable length_snap  : integer range 0 to 1000 := 0;

	*/

	int old_delay, new_delay, length_snap;
	int delay = fluxes[start_pos]; // initial clock estimate ???

	for (size_t i = start_pos; i < end_pos; ++i) {
		/*
		begin
	 
	    	if (rdata_pulse = '1' and rdata_pulse'event ) then
	 
	        if (read_front /= read_back and read_cancel = '0') then
	    */

		// pulse edge, store delay in a variable
		int old_delay = delay;

		one_imp_low = old_delay - old_delay / 4;
		one_imp_nom = old_delay;
		one_imp_high = old_delay + old_delay / 4;

		two_imp_low = old_delay + old_delay / 4;
		two_imp_nom = old_delay + old_delay / 2;
		two_imp_high = old_delay + old_delay - old_delay / 4;

		tre_imp_low = old_delay + old_delay - old_delay / 4;
		tre_imp_nom = old_delay + old_delay;
		tre_imp_high = old_delay + old_delay + old_delay / 4;

		int length = fluxes[i]; // is this what length is ???

		length_snap = length; // <--- ????

		/*
	            -- process MFM bits
	            if (skip > 0) then
	                skip <= skip - 1;
	            elsif (length_snap < delay/4) then
	                -- too small
	            elsif (length_snap < one_imp_high) then
	                -- 1us between two pulses
	                new_delay := old_delay - (old_delay - length_snap)/2;
	                mfm_type <= 2;
	                mfm_front <= mfm_front + 1;
	            elsif (length_snap < two_imp_high) then
	                -- 1.5us between two pulses
	                mfm_type <= 3;
	                mfm_front <= mfm_front + 1;
	                new_delay := old_delay - (old_delay-(length_snap+length_snap)/4-
	                                (length_snap+length_snap)/8 +
	                                (length_snap+length_snap)/32)/2;
	            elsif (length_snap < delay+delay+delay) then
	                -- 2us between two pulses
	                mfm_type <= 4;
	                mfm_front <= mfm_front + 1;
	                new_delay := old_delay - (old_delay - length_snap/2)/2;
	            else
	                -- too long duration between pulses
	                skip <= 5;
	            end if;
	    */

	    // Not entirely sure if this is what's going on
	    // "skip" is just outputting a zero, I think

	    bool emit = false;
	    int mfm_type = -1;
	    int skip = 0;
	    new_delay = delay;
	    int nd_simplified;

	    if (length_snap < delay/4) {
	    	// nothing
	    } else if (length_snap < one_imp_high) { // 1 us
	    	new_delay = old_delay - (old_delay - length_snap) / 2;
	    	mfm_type = 2;
	    	emit = true;
	    } else if (length_snap < two_imp_high) { // 1.5 us
	    	emit = true;
	    	mfm_type = 3;
	    	// WTF?
	    	new_delay = old_delay - (old_delay-(length_snap+length_snap)/4-
	                                (length_snap+length_snap)/8 +
	                                (length_snap+length_snap)/32)/2;
	    	nd_simplified = old_delay - (old_delay-
    								length_snap/2-
	                                length_snap/4 +
	                                length_snap/16)/2;
	    	if (nd_simplified != new_delay) {
	    		std::cout << new_delay << "\t" << nd_simplified << "\n";
	    		throw std::logic_error("oops");
	    	}
	    } else if (length_snap < delay + delay + delay) { // 2 us
	    	emit = true;
	    	mfm_type = 4;
	    	new_delay = old_delay - (old_delay - length_snap/2)/2;
	    } else { // Too long a duration
	    	skip = 5; // ???
	    }

	    if (emit) {
	    	for (int j = 0; j < mfm_type-1; ++j) {
				train.data.push_back(0);
				train.flux_indices.push_back(i);
			}
			train.data.push_back(1);
			train.flux_indices.push_back(i);
	    }

	    /*
	 
	            -- filter and update delay:
	            if (new_delay > old_delay + old_delay/8) then
	                delay <= old_delay + old_delay/8;
	            elsif (new_delay < old_delay - old_delay/8) then
	                delay <= old_delay - old_delay/8;
	            else
	                delay <= new_delay;
	            end if;
	    */
		if (new_delay > old_delay + old_delay / 8) {
			delay = old_delay + old_delay / 8;
		} else if ( new_delay < old_delay - old_delay / 8) {
			delay = old_delay - old_delay / 8;
		} else {
			delay = new_delay;
		}
	}

	RMSE_out = std::sqrt(RMSE_out / (end_pos - start_pos));

	return train;

		/*
	            -- initialize length for next pulse
	            reset_front <= reset_front + 1;
	 
	        else
	 
	            -- initialize skip to 1
	            skip       <= 1;
	            -- initialize delay to 100
	            delay      <= 100;
	 
	        end if;
	 
	    end if;
	 
	end process;
	*/

}


/*
MFM_train_data historical_EWMA_decoder::get_MFM_train(
		const std::vector<int> & fluxes, size_t start_pos,
		size_t end_pos, double & RMSE_out) const {

	MFM_train_data train;
	
	train.data.push_back(1);
	train.flux_indices.push_back(0);

	RMSE_out = 0;
	double estimated_half_clock = initial_clock / 2;

	for (size_t i = start_pos; i < end_pos; ++i) {
		int half_clocks = round(fluxes[i]/estimated_half_clock);

		double error_term = fluxes[i] - half_clocks * estimated_half_clock;
		RMSE_out += error_term * error_term;

		if (half_clocks == 0 || half_clocks > 4) { continue; }

		// The VHDL code has alpha = 0.5, but we'll let it be
		// a parameter.
		double proposed = estimated_half_clock * (1-alpha) + 
			fluxes[i]/(double)half_clocks * alpha;

		// Outlier thresholding: if the deviation is too far, ignore.
		if (proposed > 2.25 * estimated_half_clock) {
			proposed = 2.25 * estimated_half_clock;
		}
		if (proposed < 0.75 * estimated_half_clock) {
			proposed = 0.75 * estimated_half_clock;	
		}

		estimated_half_clock = proposed;

		int zeroes = std::max(0, half_clocks-1);

		if (zeroes == 0) { continue; }

		for (int j = 0; j < zeroes; ++j) {
			train.data.push_back(0);
			train.flux_indices.push_back(i);
		}
		train.data.push_back(1);
		train.flux_indices.push_back(i);
	}

	RMSE_out = std::sqrt(RMSE_out / (end_pos - start_pos));

	return train;
}
*/