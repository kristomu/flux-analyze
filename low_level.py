import numpy as np

def get_pulse_deltas(au_filename):
	file = open(au_filename, "rb")
	byte = b'x'
	counter = 0
	pulse_deltas = []
	while byte != b"":
		byte = file.read(1)
		if byte == b"\x7f":
			pulse_deltas.append(counter)
			counter = 0
		else:
			counter = counter + 1

	return np.array(pulse_deltas)

# Assign each point to the cluster closest to it. Returns the assignment
# of points to clusters, as well as the residuals (signed distances).
# Works with local clusters if cluster_pts is the transpose of the local
# clusters list (as a matrix); see assign_to_local_clusters.
# TODO? "Soft assignments" indicating how far the point is from the next
# closest cluster, or from some reasonable distance -- to detect dropped
# flux reversals in the case of very long delays (e.g. 80 most likely two
# 40 pulses with one of them having been lost by magnetic decay).
def assign_with_residuals(pulse_deltas, cluster_pts):
	# First assign everything to cluster 0.
	residuals = pulse_deltas - cluster_pts[0]
	assignment = np.array([0]*pulse_deltas)

	# Then assign to other clusters.
	for i in range(1, len(cluster_pts)):
		residuals_this = pulse_deltas - cluster_pts[i]
		assignment[np.abs(residuals_this) < np.abs(residuals)] = i
		residuals[np.abs(residuals_this) < np.abs(residuals)] = \
			residuals_this[np.abs(residuals_this) < np.abs(residuals)]

	return assignment, residuals

def assign(pulse_deltas, cluster_pts):
	return assign_with_residuals(pulse_deltas, cluster_pts)[0]

def get_residual_abs_sum(pulse_deltas, cluster_pts):
	assignment, residuals = assign_to_clusters(pulse_deltas, cluster_pts)
	return np.sum(np.abs(residuals))

# The same, but with dynamic clusters as returned by local_k_median.
def assign_locally(pulse_deltas, local_clusters_per_point):
	return assign_with_residuals(pulse_deltas, local_clusters_per_point.T)

# From a band assignment, produce a binary flux stream where a flux reversal
# is represented by a 1, and no reversal is represented by a 0. Also return
# a positions index that, for each bit (reversal or not), gives the first
# assignment position responsible for this bit.
def get_flux_stream_and_pos(assignments):
	# TODO? Probabilistic judgment of where a possible extra reversal
	# might occur?

	flux_code = []
	position = []
	counter = 0

	for assignment in assignments:
		# Category X corresponds to 1+X no-flux-reversals, then a reversal.
		flux_code += [0] * (assignment + 1) + [1]
		position += [counter] * (assignment + 2)
		counter += 1

	return flux_code, position

# Do the reverse. Note that the flux stream *must* start with a reversal (R).
# Any out-of-spec oddities (consecutive reversals, or more than the maximum
# allowed length non-reversal) leads to an exception.

# The default max run length is set according to MFM specs. For FM, the
# max length would be 2, and other RLL encoding types have the max length
# as part of the spec.
def encode_flux_stream(flux_stream, max_run_length=3):
	reversal_idx = np.where(flux_stream == 1)[0]
	if reversal_idx[0] != 0:
		raise Exception("encode: Not a valid flux stream.")

	assignment = reversal_idx[1:] - reversal_idx[:-1] # get deltas
	assignment -= 1 # subtract the reversal itself.

	if np.max(assignment) > max_run_length or np.min(assignment) == 0:
		raise Exception("encode: Flux stream fails clock constraints.")

	return assignment

# The following functions relate to a strategy I discovered for finding
# fixed flux patterns without knowing what the clustering is. This strategy
# can also be used to determine the bands.

# The idea is that, although we may not know the exact pulse delays that
# go into each band, we know that if a flux stream's encoding implies that
# the first delay belongs to the first band and the second to the second,
# then the second delay point must be greater than the first.

# For instance, for detecting the out-of-band A1A1A1 preamble in IBM MFM
# encoding without knowing the clustering, observe:

# We have a bunch of zeroes, which are encoded as RNRNRN..
# then we have the pattern NRNNNRNNRNNNRNNR thrice.
# If we consider one zero before the actual sequence, this is equivalent
# to, in terms of delays, to

#  1   2   3    2   3    2
# RN RNN RNNN RNN RNNN RNN R...

# or, with two runs,

#  1  2    3    2   3    2   1   3   2    3    2
# RN RNN RNNN RNN RNNN RNN R N RNNN RNN RNNN RNN R...

# And this, in turn, implies that the pulse delays must go:
# increasing (1->2), increasing (2->3), decreasing (3->2), increasing (2->3)
# decreasing (3->2), decreasing (2->1), and so on.

# And this ordinal pattern is something we *can* search for without knowing
# any clustering. However, also note that we can't search for stretches of
# identical delays (e.g. the RNRNRN zeroes), because noise is very likely
# to impart a random increasing/decreasing pattern on the pulse delays. So
# any such sequence must be treated as a wildcard when searching for
# patterns.

# Returns a string where 1 means an increase and 0 means no increase.
# If wildcard is set, then any identical runs are replaced with the
# specified character.
def get_ordinal_pattern(flux_delays, wildcard=None):
	flux_changes = (flux_delays[1:] - flux_delays[:-1]).astype(np.int8)
	flux_changes[flux_changes>0] = ord('1')
	if wildcard is None:
		flux_changes[flux_changes<=0] = ord('0')
	else:
		flux_changes[flux_changes<0] = ord('0')
		flux_changes[flux_changes==0] = ord(wildcard)

	return flux_changes.tostring()

def search_for_pattern(haystack_flux_delays, needle_flux_delays):
	haystack_pattern = get_ordinal_pattern(haystack_flux_delays)
	needle_pattern = get_ordinal_pattern(needle_flux_delays, '.')

	# Performance hack: remove a long string of wildcards before the
	# variable part, then just subtract the length from the return
	# values, since a long string of wildcards matches anything.
	for i in range(len(needle_pattern)):
		if i != '.':
			# https://stackoverflow.com/questions/4664850
			matches = [m.start() for m in re.finditer(needle_pattern[i:],
					haystack_pattern)]
			return np.array(matches) + i

	return None # matches everything

# Should perhaps be in MFM.
def get_potential_preambles_new(pulses):
	preamble_locations = np.array([])
	for flux_stream in A1_prefix_code, C2_prefix_code:
		needle_delays = encode_flux_stream(np.array(flux_stream))
		preamble_locations = np.concatenate([preamble_locations,
			search_for_pattern(pulses, needle_delays)])

	return np.array(sorted(preamble_locations))