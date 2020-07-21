# Dewarping functions.

# TODO: Incorporate into the new decoding process.

# warped.au
# (These are old values, fix later)
# No PLL
# Valid CRC chunks: 2072	Invalid: 50	Fraction valid: 0.976
# Number of sectors: 21	Bad sectors: 0
#
# Simplest PLL, alpha=0.10
# Valid CRC chunks: 2108  Invalid: 13     Fraction valid: 0.994
# Number of sectors: 21   Bad sectors: 0
#
# Comprehensive dewarping with Lloyd's
# Valid CRC chunks: 2122	Invalid: 0	Fraction valid: 1.000
# Number of sectors: 21	Bad sectors: 0
#
# Another_PLL
# Valid CRC chunks: 2121	Invalid: 0	Fraction valid: 1.000
# Number of sectors: 21	Bad sectors: 0
#
# So I know how to do this, now. Sort of.

import numpy as np
from collections import deque

def simplest_ever_pll(pulses, clustering, alpha=0.05):
	assignments = []

	# https://www.me.utexas.edu/~jensen/ORMM/omie/operation/unit/forecast/exp.html
	lag = int((1-alpha)/alpha)

	for i in range(len(pulses)):
		assignment = np.argmin(np.abs(pulses[i] - clustering))
		residual = pulses[i] - clustering[assignment]

		# Relative residual
		rel_residual = residual / clustering[assignment]

		# Adjust the pulse delays accordingly.
		clustering = clustering + clustering * rel_residual * alpha

		if i >= lag:
			assignment_earlier = np.argmin(np.abs(pulses[i-lag] - clustering))
			assignments.append(assignment_earlier)

	for i in range(len(pulses)-lag, len(pulses)):
		assignment = np.argmin(np.abs(pulses[i] - clustering))
		assignments.append(assignment)

	return assignments

def another_pll(pulses, x=None, alpha=0.05):
	if x is None:
		x = kmedian.wt_optimal_k_median(pulses[:3000], 3)
	assignments = []

	for i in range(31, len(pulses)-15):
		residuals = []
		for j in range(-15, 15):
			assignment = np.argmin(np.abs(pulses[i+j] - x))
			residual = pulses[i+j] - x[assignment]
			rel_residual = residual / x[assignment]

		# Relative residual
		rel_residual = np.median(rel_residual)

		# Adjust the pulse delays accordingly.
		x = x + x * rel_residual * alpha

		assignment = np.argmin(np.abs(pulses[i] - x))
		assignments.append(assignment)
		#residual = pulses[i] - x[assignment]

	return assignments

# Drop that alpha and use a much larger window instead?
def shifted_window_pll(pulses, start_pos, end_pos, x=None, alpha=0.005):
	if x is None:
		x = kmedian.wt_optimal_k_median(pulses[start_pos:start_pos+3000], 3)
	assignments = []

	window_size = 55
	left_of_center = window_size//2
	right_of_center = window_size - 1 - left_of_center

	rel_residuals = deque()

	for i in range(start_pos - left_of_center, start_pos - left_of_center + window_size):
		assignment = np.argmin(np.abs(pulses[i] - x))
		residual = pulses[i] - x[assignment]
		rel_residuals.append(residual / x[assignment])

	# Get the mean relative residual
	mean_rel_residual = np.mean(rel_residuals)

	for i in range(start_pos - left_of_center, end_pos - left_of_center):
		# Adjust the pulse delays accordingly.
		x = x + x * mean_rel_residual * alpha

		# Assign the point in the center of the window
		assignment = np.argmin(np.abs(pulses[i+left_of_center] - x))
		assignments.append(assignment)

		# Assign the point on the right, shift it in, remove
		# what's to the left
		assignment = np.argmin(np.abs(pulses[i+window_size] - x))
		residual = pulses[i+window_size] - x[assignment]
		new_rel_residual = residual / x[assignment]
		old_rel_residual = rel_residuals.popleft()
		rel_residuals.append(new_rel_residual)

		mean_rel_residual -= old_rel_residual / window_size
		mean_rel_residual += new_rel_residual / window_size

	return assignments