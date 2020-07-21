import scipy.stats
import numpy as np
import re

# I should probably change the default categorization from clusters to
# intervals, because suppose I have intervals like this:

# [1, 2] [10, 20] [100, 200]

# then the cluster centers won't indicate the means/median position in
# each band. This will mess up the residuals as well as anything based
# on soft assignment. It isn't absolutely critical yet, though.

# Now somewhat faster!
def get_pulse_deltas(au_filename, bufsize=16384):
	file = open(au_filename, "rb")
	in_bytes = b'x'
	all_pulse_locations = []
	counter = 0
	while in_bytes != b"":
		in_bytes = file.read(bufsize)
		in_array = np.frombuffer(in_bytes, dtype=np.int8)
		pulse_locations = np.where(in_array==127)[0] + counter
		all_pulse_locations.append(pulse_locations)
		counter += bufsize

	pulse_locations = np.concatenate(all_pulse_locations)

	# -1 because distances between pulses are exclusive, i.e. don't
	# count the \x7f itself. Shouldn't matter much, really...
	return pulse_locations[1:] - pulse_locations[:-1] - 1

# Assign each point to the cluster closest to it. Returns the assignment
# of points to clusters, as well as the residuals (signed distances).
# Works with local clusters if cluster_pts is the transpose of the local
# clusters list (as a matrix); see assign_to_local_clusters.
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

# I'm planning to use the same parameter format for assignment functions
# that also do dewarping.
def assign_partial(pulse_deltas, start_pos, stop_pos, cluster_pts):
	return assign(pulse_deltas[start_pos:stop_pos], cluster_pts)

def get_residual_abs_sum(pulse_deltas, cluster_pts):
	assignment, residuals = assign_to_clusters(pulse_deltas, cluster_pts)
	return np.sum(np.abs(residuals))

# The same, but with dynamic clusters as returned by local_k_median.
def assign_locally(pulse_deltas, local_clusters_per_point):
	return assign_with_residuals(pulse_deltas, local_clusters_per_point.T)

# Returns an array where the first value is the index of the center
# closest to the value, the second is the index of the center next to
# closest, and so on.
# Passing the transpose of a numpy array as value also works and will
# give the indices for each value in that array.
def get_full_assignment(value, centers):
	return np.argsort(np.abs(value - centers))

# From a band assignment, produce a binary flux stream where a flux reversal
# is represented by a 1, and no reversal is represented by a 0. Also return
# a positions index that, for each bit (reversal or not), gives the first
# assignment position responsible for this bit.
# TODO? Preallocate arrays. But that could be hard indeed.
def get_flux_stream_and_pos(assignments):

	# It's implied that there's a flux transition before the first delay,
	# because each delay counts the time from the last to the next
	# transition.
	flux_code = [1]
	position = [0]

	counter = 0

	for assignment in assignments:
		# Category X corresponds to 1+X no-flux-reversals, then a reversal.
		# I don't think this is right. TODO: Fix.
		flux_code += [0] * (assignment + 1) + [1]
		position += [counter] * (assignment + 2)
		counter += 1

	return flux_code, position

# Do the reverse. Note that the flux stream *must* start with a reversal (R),
# because otherwise we don't know what band the first flux delay belongs
# to (e.g. the pattern "0 1" could be "1 0 1" or "1 0 0 1", and we
# wouldn't know which).
# Any out-of-spec oddities (consecutive reversals, or more than the maximum
# allowed length non-reversal) will produce an exception.

# The default max run length is set according to MFM specs. For FM, the
# max length would be 2, and other RLL encoding types have the max length
# as part of the spec.

# Since the pulse delay output array has values 1, 2, etc. instead of
# actual timings, it's called an assignment here, because that's what it
# looks like.
def encode_to_assignment(flux_stream, max_run_length=3):
	reversal_idx = np.where(flux_stream == 1)[0]
	if reversal_idx[0] != 0:
		raise Exception("encode: Not a valid flux stream.")

	assignment = reversal_idx[1:] - reversal_idx[:-1] # get deltas
	# Subtract the reversal itself...
	assignment -= 1

	if np.max(assignment) > max_run_length or np.min(assignment) == 0:
		raise Exception("encode: Flux stream fails clock constraints.")

	# ... and once more, because assignments are zero-indexed.
	assignment -= 1

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

# Strictly speaking, if the noise overwhelms the signal or we're in a trend
# (e.g. a warping of the floppy surface), the above reasoning won't work,
# but then we're SOL anyway.

# Returns a byte rendering of an np array denoting either increase or
# no increase of the flux delay in question. The first element gives
# whether the the second delay is greater than the first, and so on.
def get_ordinal_pattern(flux_delays):
	flux_changes = (flux_delays[1:] - flux_delays[:-1]) > 0
	return flux_changes.tostring()

# Returns a tuple where the first element gives the number of wildcard
# elements at the start of the pattern, and the second gives the non-
# wildcard pattern that follows it.
# If there are runs of zero difference inside the pattern (not just
# at the beginning, this will throw an exception.
def get_ordinal_pattern_with_wildcards(flux_delays):
	flux_changes = (flux_delays[1:] - flux_delays[:-1])
	wildcard_count = np.where(flux_changes)
	for wildcard_count in range(len(flux_changes)):
		if flux_changes[wildcard_count] != 0:
			if np.any(flux_changes[wildcard_count:] == 0):
				raise Exception("GOPWW: Wildcard in middle of ordinal string!")

			return (wildcard_count,
				(flux_changes[wildcard_count:]>0).tostring())

def ordinal_search_for_pattern(haystack_pattern, needle_flux_delays):
	needle_wild_count, needle_pattern = get_ordinal_pattern_with_wildcards(needle_flux_delays)

	# https://stackoverflow.com/questions/4664850
	matches = np.array([m.start() for m in re.finditer(needle_pattern,
		haystack_pattern)])

	# Adjust all match locations to account for the initial run
	# of wildcards, if any. Cast to integer because positions are,
	# by their nature, integers.
	return (matches - needle_wild_count).astype(np.int)

# Pruning: suppose we've found some ordinal pattern matches derived from an
# assignment flux pattern of [1, 3, 2, 1]. Then we know that the first pulse
# delay in the match belongs to band 1, the second to band 3 and so on, which
# automatically gives us band clustering ranges. In addition, if the thus
# defined intervals overlap each other, we know that the ordinal match must
# be a false one - there's no way to tell the bands apart if the sequence is
# true, and since the bands should be distinct, that means the sequence is
# false.

# The functions below separate out bands and check if the intervals overlap.
# Note, severe noise/degradation may cause a sequence to appear overlapping
# when it really isn't (absent the noise).

def classify_bands(haystack_flux_delays, needle_match_pos,
	needle_assignments, num_bands=3):

	points_by_band = [np.array([]) for x in range(num_bands)]

	for i in range(num_bands):
		# Assignments are zero-indexed
		positions_with_this_band = np.where(needle_assignments==i)[0] + \
			needle_match_pos
		points_by_band[i] = haystack_flux_delays[
			positions_with_this_band]

	return points_by_band

# Get appropriate cluster intervals. It returns one (min, max) tuple per band,
# or None if it's impossible to assign intervals to tell the bands apart.

# If it doesn't return None, then row k of the matrix gives the inclusive
# interval that must be assigned to the kth band.
def get_band_intervals(points_by_band):
	if len(points_by_band) == 0:
		return None

	band_intervals = np.array([[min(x), max(x)] for x in points_by_band])
	for i in range(len(band_intervals)-1):
		# If our max is greater than the next one's min, there's an
		# overlap and no classification can possibly work.
		if band_intervals[i][1] >= band_intervals[i+1][0]:
			return None
	return band_intervals

def get_valid_ordinal_matches(haystack, needle, matches):
	filtered_matches = [ x for x in matches if get_band_intervals(
		classify_bands(haystack, x, needle)) is not None]

	return filtered_matches

def search_for_pattern(haystack_flux_delays, needle_assignments):
	haystack_pattern = get_ordinal_pattern(haystack_flux_delays)

	# First do an ordinal search...
	ordinal_matches = ordinal_search_for_pattern(haystack_pattern,
		needle_assignments)

	# Then filter out those that produce overlapping cluster intervals.
	filtered_matches = get_valid_ordinal_matches(haystack_flux_delays,
		needle_assignments, ordinal_matches)

	# Explicitly cast to int because empty np arrays default to float.
	return np.array(filtered_matches).astype(np.int)

def search_for_patterns(haystack_flux_delays, needles_assignments):
	return [ search_for_pattern(haystack_flux_delays, x) for x in \
		needles_assignments ]

def band_intervals_to_clusters(band_intervals):
	# The clusters have to be chosen so that the distance from the right
	# cluster to the end of the interval is smaller than the distance from
	# any other cluster to that point.

	# Suppose that the clusters, in sorted order, are x_1, x_2, x_3
	# and the closed intervals for each band is [min_1, max_1]
	# [min_2, max_2] [min_3, max_3].
	# Then one way of accomplishing the above is

	# max_1 - x_1 = x_2 - min_2
	# max_2 - x_2 = x_3 - min_3

	# or in general,

	# x_(i+1) = min_(i+1) + max_i - x_i
	# with x_1 free.

	# In addition, each of these max - ... terms must be nonnegative.
	# Otherwise, the following example fails:

	# [[17, 28],
	#  [33, 37],
	#  [38, 54]].

	# The overflow problem occurs when min_(i+1) + max_i - x_i > max_(i+1).
	# The mirror "underflow problem", min_(i+1) + max_i - x_i < min(i+1)
	# can never happen as long as max_i > x_i.
	# At first glance, that suggests maximizing x_i. However, maximizing
	# x_i will minimize x_(i+1) and thus may not be the answer.

	clusters = [[] for x in range(len(band_intervals))]

	# For now I just choose the lowest value of x_1 that passes
	# test, or the median if it's included.

	lower_bound = band_intervals[1][0] + band_intervals[0][1] - band_intervals[1][1]
	clusters[0] = max(np.median(band_intervals[0]), lower_bound)

	for i in range(1, len(band_intervals)):
		clusters[i] = band_intervals[i][0] + band_intervals[i-1][1] - clusters[i-1]

	if sorted(clusters) != clusters:
		raise Exception("Internal error: Output clusters are not sorted.")

	return np.array(clusters)

def random_cluster(band_intervals):
	# Let the three cluster points be x, y, z, and the min and max
	# for the intervals be min_i, max_i, indexed from 1.

	# Then x must be closer to max_1 than y is.
	# y must be closer to min_2 than x is and must be closer to max_2 than z is.
	# z must be closer to min_3 than y is.

	# In other words:
	# max_1 - x < y - max_1		==>		y > 2 * max_1 - x
	# y - min_2 < min_2 - x		==>		y < 2 * min_2 - x
	# max_2 - y < z - max_2		==>		z > 2 * max_2 - y
	# z - min_3 < min_3 - y		==>		z < 2 * min_3 - y

	eps = 1e-9

	x = 3
	y = 2
	z = 1

	while x >= y or y >= z:
		x = np.random.uniform(0, band_intervals[0][1])
		y = np.random.uniform(2 * band_intervals[0][1] - x + eps,
			2 * band_intervals[1][0] - x)
		z = np.random.uniform(2 * band_intervals[1][1] - y + eps,
			2 * band_intervals[2][0] - y)

	return np.array([x, y, z])

def valid_clustering(band_intervals, proposed_clustering):

	for i in range(1, len(proposed_clustering)):
		if proposed_clustering[i] <= 2 * band_intervals[i-1][1] - proposed_clustering[i-1]:
			return False
		if proposed_clustering[i] >= 2 * band_intervals[i][0] - proposed_clustering[i-1]:
			return False

	return True

# TODO: Handling A1 matches and C2 matches in the same list is turning out
# to be quite a pain. Split them up or somehow indicate *which* match it
# is.

# Gets a cluster for this particular chunk.
def get_cluster(haystack_flux_delays, needle_match_pos,
	needle_assignments, num_bands=3):

	classified_points = classify_bands(haystack_flux_delays,
		needle_match_pos, needle_assignments, num_bands)

	band_intervals = get_band_intervals(classified_points)

	return band_intervals_to_clusters(band_intervals)

# Surprisingly efficient clustering function: categorize every preamble,
# then take the median (or a trimmed mean) of each category as the cluster
# center.

def get_combined_band_categorization(haystack, needles, num_bands=3):
	bands = [np.array([]) for x in range(num_bands)]

	for needle_assignments in needles:
		for preamble_pos in search_for_pattern(haystack, needle_assignments):
			classified_points = classify_bands(haystack, preamble_pos,
				needle_assignments, num_bands)
			for j in range(len(classified_points)):
				bands[j] = np.concatenate([bands[j], classified_points[j]])

	# For ease of interpretation
	for j in range(num_bands):
		bands[j] = np.sort(bands[j])

	return bands

# Returns a trimmed midrange, used for estimating cluster centers.
# Nothing fancy, but it should do the trick.

def midsummary(x, trim=None):
	return (np.quantile(x, trim) + np.quantile(x, 1-trim)) / 2.0

def cluster_on_custom_preambles(pulses, needles, trim=None):
	if trim is None:
		return np.array([np.median(x) for x in \
			get_combined_band_categorization(pulses, needles)])

	return np.array([midsummary(x, trim) for x in \
		get_combined_band_categorization(pulses, needles)])