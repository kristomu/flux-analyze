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

		STATE:	Current flux delay period
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

typedef enum mode_t {
	IN_BAND = 0,
	PREAMBLE_A = 1,
	PREAMBLE_C = 2,
	ILLEGIBLE = 3,
	NUM_MODES = 4};

class dp_idx {
	public:
		mode_t mode;
		int clock;
		int idx;
		int parity;
};

class dp_record {
	public:
		int cur_period;
		int OOB_sequence_index;
		double penalty;
		dp_idx previous_opt;
};

const int MIN_CLOCK = 1, MAX_CLOCK = 100;

void do_dp(const & flux_record f) {

	// Noise sensitivity hyperparameter.
	// Higher values of this allows greater changes in clock
	// from one pulse to the next, while lower values prioritize fitting
	// the set clock to the current pulse delay.
	double alpha = 0.5;

	std::vector<std::vector<std::vector<std::vector<dp_record> > > > dp;
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

	// Just do something with the first track.

	return 0;
}
