"""

Do some prototyping of a dynamic programming solution for flux reconstruction, without
missing bits yet.

There are four modes:
	0 - in-band (obeys MFM rules)
	1 - out-of-band A1A1A1 preamble
	2 - out-of-band C2C2C2 preamble
	3 - everything else that doesn't fit the model ("unknown")

The state for each step is:
	INDEX: Mode
	INDEX: Current estimated clock value
	INDEX: Current flux transition (input index)

	STATE: Current flux delay period
	STATE: OOB sequence bit index (mod 16)
	STATE: accumulated penalty

I'll also keep the previous entry (the optimum we're building on) inside
the array; it makes the code easier to read. I can always remove it later
if I need the space.

For the first three modes, we're assuming the floppy is formated with a standard MFM
encoding, where 0 is no flux reversal and 1 is:

	MFM sequence		last bit	current bit
	01					N/A			1
	00					1			0
	10					0			0

and that flux transitions are based off a clock as follows:
	clock multiple		sequence
	1					01
	1.5					001
	2					0001

Doing separate bands (like my previous Python attempt) might be better, but it would
also cost a lot more.

Another concern is that the quadratic penalty makes the hyperparameters depend on the
scale of the clock, so good parameter choices for DD floppies may be different than
for HD floppies. I should find a way of handling that.
"""

import numpy as np
import pickle

MIN_CLOCK = 1
MAX_CLOCK = 100

# indices
MODE = 0
PERIOD = 1
PENALTY = 2
OPTIMAL_PREFIX = 3

# special values
NO_PRIOR_FLUX = -1
WRONG_MODE = -1

# modes
MODE_IN_BAND = 0
MODE_OUT_OF_BAND_A = 1
MODE_OUT_OF_BAND_C = 2
MODE_UNKNOWN = 3

MODE_DEFAULT = -2	# for debugging; passes everything

# noise sensitivity hyperparameter
# Higher values of this allows greater changes in clock from
# one pulse to the next, while lower values prioritize fitting
# the set clock to the current pulse delay.
alpha = 0.5

# Simple DP array would have the following indices: clock and index.
flux_delays = pickle.load(open("flux_delays_list.pickle", "rb"))[:2000]

num_fluxes = len(flux_delays)

# Dimension the DP list.
dp_list = []
for i in range(num_fluxes):
	dp_list.append([None] * (MAX_CLOCK+1))

# Use the previous DP state to repeatedly advance one step, and accumulate
# a penalty for changing the clock.

for i in range(num_fluxes):
	print(i/num_fluxes)
	for cur_clock in range(MIN_CLOCK, MAX_CLOCK+1):
		best_penalty = float("inf")
		recordholder = -1

		# Get the current pulse length based on the clock and delay,
		# as well as what the clock would be if the signal was perfectly
		# centered.
		observed_pulse_delay = flux_delays[i]
		pulse_length = np.round(2 * observed_pulse_delay / cur_clock) - 1
		# Pulse lengths must be 1, 2, or 3.
		pulse_length = int(min(3, max(1, pulse_length)))
		expected_pulse_delay = cur_clock * (pulse_length + 1)/ 2

		for previous_clock in range(MIN_CLOCK, MAX_CLOCK+1):
			# Calculate the penalty when transitioning from this
			# clock.
			candidate_penalty = (cur_clock - previous_clock)**2 + \
				alpha * (observed_pulse_delay - expected_pulse_delay)**2

			# Add the accumulated penalty up to this point.
			# If it's the first step, we incur no penalty from prior
			# steps; only from our own.

			if i != 0:
				candidate_penalty += dp_list[i-1][previous_clock][PENALTY]

			if candidate_penalty < best_penalty:
				best_penalty = candidate_penalty
				recordholder = previous_clock

		dp_list[i][cur_clock] = [MODE_DEFAULT, pulse_length, best_penalty,
			recordholder]

# Unspool the DP choices to find the pulse lengths and clocks.
last_time = np.array(dp_list[-1][MIN_CLOCK:MAX_CLOCK+1])
optimal_ending_clock = np.argmin(last_time[:,2]) + MIN_CLOCK
optimal_clock = optimal_ending_clock

clocks = []
pulse_lengths = []

for i in reversed(range(0, num_fluxes)):
	clocks.append(optimal_clock)
	pulse_lengths.append(dp_list[i][optimal_clock][PERIOD])

	optimal_clock = dp_list[i][optimal_clock][OPTIMAL_PREFIX]

clocks.reverse()
pulse_lengths.reverse()