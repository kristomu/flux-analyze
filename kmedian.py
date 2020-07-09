import numpy as np
import sortedcontainers as sc

# Implement a k * n log n method for optimal 1D K-median.
# The K-median problem consists of finding k clusters so that the
# sum of absolute distances from each of the n points to the closest
# cluster is minimized.

# GRØNLUND, Allan, et al. Fast exact k-means, k-medians and Bregman
# divergence clustering in 1d. arXiv preprint arXiv:1701.07204, 2017.

# Rather deep magic happens here, and so I've not commented everything.
# It is possible to do this even quicker at the expense of even more
# impenetrable magic; I haven't done so.

# Get the prefix sum array x (for use in CC etc). Points are given in
# ascending order.
def get_cdf(points):
	return np.cumsum(points)

# Grønlund et al.'s CC function for the K-median problem. x is the sorted
# list of points in ascending order, and px is the prefix sum from the above
# function.
# CC(x,i,j) is the cost of grouping x_i...x_j into one cluster with the
# optimal cluster point (the median point). By programming convention, the
# interval is half-open and indexed from 0, unlike the paper's convention
# of closed intervals indexed from 1.
def CC(x, px, i, j):
	if i >= j: return 0

	m = (j+i-1)/2.0

	# Lots of nasty off-by-one traps here.
	mu = (x[int(np.floor(m))] + x[int(np.ceil(m))]) / 2		# median

	# Lower part: entry i to median point between i and j, median excluded.
	sum_below = px[int(np.floor(m))]
	if i > 0:
		sum_below -= px[i-1]
	# Upper part: everything from the median up, median included if it's a
	# real point.
	sum_above = px[j-1] - px[int(np.floor(m))]

	return np.floor(m - i + 1) * mu - sum_below + \
		   sum_above - np.ceil(j - m-1) * mu

# the kn log n algorithm.
# D_previous is the D vector for (i-1) clusters, or empty if i < 2.
# It's possible to do this even faster (and more incomprehensibly).
# See Grønlund et al. for that.

# Calculate C_i[m][j] given D_previous = D[i-1]. p. 4
def C_i(x, px, i, D_previous, m, j):
	if i == 1:
		return CC(x, px, 0, m)
	if i == 2:
		return CC(x, px, 0, min(j, m)) + CC(x, px, j, m)
	else:
		return D_previous[min(j, m)] + CC(x, px, j, m)

# Find the (location of the) minimum value of each row of an nxn matrix in
# n log n time, given that the matrix is monotone (by the definition of
# Grønlund et al.) and defined by a function M that takes row and column
# as parameters. All boundary values are closed, [min_row, max_row]
# p. 197 of
# AGGARWAL, Alok, et al. Geometric applications of a matrix-searching
# algorithm. Algorithmica, 1987, 2.1-4: 195-208.
def monotone_matrix_indices(M, min_row, max_row, min_col, max_col,
	Tout, Dout):
	if max_row == min_row: return

	cur_row = int(min_row + np.floor((max_row-min_row)/2.0))

	# Find the minimum column for this row.
	min_here = min_col + np.argmin([M(cur_row, j) for j in \
			range(min_col, max_col+1)])

	if Tout[cur_row] != 0:
		raise Exception("tried to set a variable already set.")

	Tout[cur_row] = min_here
	Dout[cur_row] = M(cur_row, min_here)

	# No lower row can have a minimum to the right of the current minimum.
	# Recurse on that assumption.
	monotone_matrix_indices(M, min_row, cur_row, min_col, min_here+1, Tout, Dout)
	# And no higher row can have a minimum to the left.
	monotone_matrix_indices(M, cur_row+1, max_row, min_here, max_col, Tout, Dout)

def optimal_k_median(points, num_clusters):
	pts_sorted = sorted(points)
	pts_cumulative = get_cdf(pts_sorted)
	pts_len = len(pts_sorted)

	T, D = [[]], [[]]

	for i in range(1, num_clusters+1):
		# If we have achieved an optimal clustering with fewer clusters,
		# our work is done.
		if len(D[-1]) > 0 and D[-1][-1] == 0:
			continue

		M = lambda m, j: C_i(pts_sorted, pts_cumulative, i, D[-1], m, j)

		Tout = np.array([0] * (pts_len+1))
		Dout = np.array([0.0] * (pts_len+1))

		monotone_matrix_indices(M, i, pts_len+1, i-1, pts_len+1, Tout, Dout)

		T.append(Tout)
		D.append(Dout)

	# Backtrack. The last cluster has to encompass everything, so the
	# rightmost boundary of the last cluster is the last point in the array,
	# hence given by the last position of T. Then the previous cluster
	# transition boundary is given by where that cluster starts, and so on.

	cur_clustering_range = len(points)
	cluster_boundaries = []
	centers = []

	for i in range(num_clusters, 0, -1):
		print(i)
		cluster_boundaries.append(T[i][cur_clustering_range])
		new_cluster_range = T[i][cur_clustering_range]
		centers.append(
			np.median(pts_sorted[new_cluster_range:cur_clustering_range]))
		cur_clustering_range = new_cluster_range

	return np.array(sorted(centers))

# --- WEIGHTED STUFF ---

# The points all have one of a very small range of values (compared to the
# number of points). Since k-medians clustering is either a point or the
# mean of two points, we can then do calculations with n being the number
# of distinct points, rather than the number of points, so that O(kn log n)
# becomes much quicker in concrete terms.

# The cdf (or prefix sum) is a 2D numpy array, each row being:
# [number of points seen so far, sum of point values up to and including
#  this point value, this point value].

# Finding the median takes O(log n) time instead of O(1) time, which
# increases the constant factor in front of the kn log n term.

# Construct the prefix sum by using a counting sort. Since the range
# is very small, this cuts out time down from the O(n log n) required
# to sort, to O(n) with a manageable constant factor.
def wt_get_prefix_sum(pts):
	count_per_value = np.bincount(pts)

	prefix_sum = []
	sum_values = 0
	sum_counts = 0
	for cur_value in range(len(count_per_value)):
		if count_per_value[cur_value] == 0:
			continue

		sum_counts += count_per_value[cur_value]
		sum_values += cur_value * count_per_value[cur_value]
		prefix_sum.append((sum_counts, sum_values, cur_value))

	return np.array(prefix_sum)

# Determines the median point of all points between the ith category
# in weighted_cdf (inclusive), and the jth (exclusive). It's used to
# reconstruct clusters later.
def wt_median(weighted_cdf, i, j):
	if i >= j: return 0

	lower_count = 0
	if i > 0:
		lower_count = weighted_cdf[i-1][0]
	upper_count = weighted_cdf[j-1][0]

	# Note, this is not the index, but number of points start inclusive,
	# which is why the -1 has to be turned into +1.
	median_pt = (lower_count+upper_count+1)/2.0

	lower_median_idx = weighted_cdf[0:,0].searchsorted(
		np.floor(median_pt))

	if (int(np.ceil(median_pt)) == int(np.floor(median_pt))) or \
		int(np.floor(median_pt)) != weighted_cdf[lower_median_idx][0]:
		return weighted_cdf[lower_median_idx][2]
	else:
		return (weighted_cdf[lower_median_idx][2] + \
			weighted_cdf[lower_median_idx+1][2]) / 2.0

def wt_median_test():
	return (wt_median(wt_get_cdf([0, 0, 1, 1])) == 0.5) \
		and (wt_median(wt_get_cdf([0, 1, 3])) == 1)

# If weighted_cdf[index] is the first weighted_cdf value with count at
# or above i, returns the cumulative sum up to the ith entry of the
# underlying sorted points list.
def wt_cumulative_at(weighted_cdf, index, i):
	sum_below = 0
	num_points_below = 0

	if index > 0:
		sum_below = weighted_cdf[index-1][1]
		num_points_below = weighted_cdf[index-1][0]

	return sum_below + weighted_cdf[index][2] * (i - num_points_below)

# Note: i and j are indices onto weighted_cdf. So e.g. if i = 0, j = 2 and
# weighted cdf is [[2, 0, 0], [4, 2, 1], [5, 2, 2]], then that is the cost
# of clustering all points between the one described by the zeroth weighted
# cdf entry, and up to (but not including) the last. In other words, it's
# CC([0, 0, 1, 1], 0, 5).
def wt_CC(weighted_cdf, i, j):
	if i >= j: return 0

	lower_count = 0
	if i > 0:
		lower_count = weighted_cdf[i-1][0]
	upper_count = weighted_cdf[j-1][0]

	median_pt = (lower_count+upper_count+1)/2.0

	lower_median_idx = weighted_cdf[0:,0].searchsorted(
		np.floor(median_pt))

	if (int(np.ceil(median_pt)) == int(np.floor(median_pt))) or \
		int(np.floor(median_pt)) != weighted_cdf[lower_median_idx][0]:
		mu = weighted_cdf[lower_median_idx][2]
	else:
		mu = (weighted_cdf[lower_median_idx][2] + \
			weighted_cdf[lower_median_idx+1][2]) / 2.0

	# Lower part: entry i to median point between i and j, median excluded.
	sum_below = wt_cumulative_at(weighted_cdf, lower_median_idx,
		np.floor(median_pt))
	if i > 0:
		sum_below -= weighted_cdf[i-1][1]
	# Upper part: everything from the median up, median included if it's a
	# real point.
	sum_above = weighted_cdf[j-1][1] - wt_cumulative_at(weighted_cdf,
		lower_median_idx, np.floor(median_pt))

	return np.floor(median_pt - lower_count) * mu - sum_below + \
		   sum_above - np.ceil(upper_count - median_pt) * mu

# A test that weighted CC is the same as unweighted.
def test_wt_CC():
	for i in range(1000):
		test = [np.random.randint(5) for j in range(9)]
		wtcdf = wt_get_cdf(sorted(test))
		wtCC = wt_CC(wtcdf, 1, len(wtcdf))
		ordCC = CC(sorted(test), get_cdf(sorted(test)), wtcdf[0][0],
			len(test))
		if wtCC != ordCC:
			print(test)
			return False

	return True

# These are pretty much cut and paste of the unweighted version.

def wt_C_i(weighted_cdf, i, D_previous, m, j):
	if i == 1:
		return wt_CC(weighted_cdf, 0, m)
	else:
		return D_previous[min(j, m)] + wt_CC(weighted_cdf, j, m)

def wt_optimal_k_median_from_cumul(wt_pts_cumulative, num_clusters):
	unique_points = len(wt_pts_cumulative)

	T, D = [[]], [[]]

	for i in range(1, num_clusters+1):
		# If we have achieved an optimal clustering with fewer clusters,
		# our work is done.
		if len(D[-1]) > 0 and D[-1][-1] == 0:
			continue

		# Otherwise, find the optimal T (cluster boundaries) and D
		# (cluster costs) using the Aggarwal algorithm.
		M = lambda m, j: wt_C_i(wt_pts_cumulative, i, D[-1], m, j)

		Tout = np.array([0] * (unique_points+1))
		Dout = np.array([0.0] * (unique_points+1))

		monotone_matrix_indices(M, i, unique_points+1, i-1, unique_points+1,
			Tout, Dout)

		T.append(Tout)
		D.append(Dout)

	# Backtrack. The last cluster has to encompass everything, so the
	# rightmost boundary of the last cluster is the last point in the array,
	# hence given by the last position of T. Then the previous cluster
	# transition boundary is given by where that cluster starts, and so on.

	cur_clustering_range = len(wt_pts_cumulative)
	cluster_boundaries = []
	centers = []

	for i in range(num_clusters, 0, -1):
		cluster_boundaries.append(T[i][cur_clustering_range])
		new_cluster_range = T[i][cur_clustering_range]
		# Reconstruct the cluster that's the median point between
		# new_cluster_range and cur_clustering_range.
		centers.append(wt_median(wt_pts_cumulative, new_cluster_range,
			cur_clustering_range))
		cur_clustering_range = new_cluster_range

	return np.array(sorted(centers))

def wt_optimal_k_median(points, num_clusters):
	return wt_optimal_k_median_from_cumul(wt_get_prefix_sum(points),
		num_clusters)

# Assign each point to the cluster closest to it. Returns the assignment
# of points to clusters, as well as the residuals (signed distances).
# Works with local clusters if cluster_pts is the transpose of the local
# clusters list (as a matrix); see assign_to_local_clusters.
# TODO? "Soft assignments" indicating how far the point is from the next
# closest cluster, or from some reasonable distance -- to detect dropped
# flux reversals in the case of very long delays (e.g. 80 most likely two
# 40 pulses with one of them having been lost by magnetic decay).
def assign_to_clusters(pulse_deltas, cluster_pts):
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

def sum_distances(pulse_deltas, cluster_pts):
	assignment, residuals = assign_to_clusters(pulse_deltas, cluster_pts)
	return np.sum(np.abs(residuals))

# Adaptive k-median (dewarping).

# If the disk has a wobble, is warped, or just a mechanical snag somewhere,
# then everything on a certain area of the plastic might register with
# consistently higher (or lower) pulse delays than elsewhere. Correcting
# such a warp seems to improve the chance of restoring data in edge cases.

# I've chosen a strategy of using a sliding window centered on the point to
# be corrected. If points in the neighborhood are all higher than usual,
# then we'll correct that point down, and vice versa. Just like the code
# above, I'll use k-median for this.

# There is a snag, however. It may happen that every pulse within the
# sliding window is the same multiple of the clock. If so, directly using
# k-median would place all three clusters near the band in question, and
# points that rightly belong to the same band would be classified as three
# different ones.

# To fix that, we'll first use an optimal k-median on the sliding window
# area around the first point. Then we'll compare the points to the global
# k-median; if more than one of the local clusters are assigned to the same
# global cluster, then we replace the one furthest away with the global
# cluster that's not represented. We'll repeat until done. That should
# represent the bands that aren't present inside the sliding window.

# Then for further points, we just have to keep the invariant that every
# band has a cluster assigned to it, even when that band's not physically
# present. So my strategy is to, when we slide the window one pulse to the
# right, pop off the old value, push in the new value, and run Lloyd's
# algorithm on the clusters from the previous point. That will leave any
# clusters with no representation alone, while moving clusters that are
# represented according to the general shift in delays.

# I've sped up Lloyd's algorithm by using a sorted structure with O(log n)
# insertion. In C++, that can be done with a map.

# Known bugs: The dewarper can sometimes get confused in areas with a lot
# of corruption. (I've added a quick and dirty fix.) Also quite sloow.

# Window size should be odd: the center gets one point, and then equal sized
# wings to either direction of the center.

def local_k_median(pts, start_center=0, window_size=301, test_interval=801,
	global_clusters=None, num_clusters=3):

	# Window sides
	left_of_center = window_size//2
	right_of_center = window_size - 1 - left_of_center

	# One for each point in the pts array.
	local_clusters_per_point = []

	# Get the global clusters in sorted order.
	if global_clusters is None:
		global_clusters = wt_optimal_k_median(pts, num_clusters)

	# Since there's nothing to the left of the leftmost point, start with
	# a window that includes the center plus points to the right.
	# If a different start point has been specified, also include points
	# to the left (if any).
	left_side = max(0, start_center - left_of_center)
	right_side = min(len(pts), start_center + right_of_center + 1)

	local_clusters = wt_optimal_k_median(pts[left_side:right_side],
		num_clusters)

	# Find out which local cluster is closest to the kth global cluster,
	# so we can augment and replace if some global clusters are missing.
	closest_local_clusters = [(np.inf, np.inf)] * num_clusters
	for i in range(num_clusters):
		distances = np.abs(local_clusters[i] - global_clusters)
		closest_global, distance = np.argmin(distances), np.min(distances)
		candidate_cluster_info = (i, distance)

		if candidate_cluster_info[1] < closest_local_clusters[closest_global][1]:
			closest_local_clusters[closest_global] = candidate_cluster_info

	local_clusters_augmented = []
	for i in range(num_clusters):
		if closest_local_clusters[i][0] == np.inf:
			local_clusters_augmented.append(global_clusters[i])
		else:
			local_clusters_augmented.append(local_clusters[closest_local_clusters[i][0]])

	local_clusters = np.array(local_clusters_augmented)
	local_clusters_per_point.append(local_clusters_augmented)

	# We now have a suitable local cluster group for the first point.
	# For our invariant, insert every point we just classified into
	# a sliding window, and then start from the second point on.

	classified_points = sc.SortedList(pts[:(right_of_center+1)])
	test_counter = 0

	for i in range(start_center+1, len(pts)):
		# If we're at a test interval, check if the global clustering fits
		# better than the last local cluster. This may happen if the
		# invariant is violated at some point, e.g. Lloyd's gets confused
		# due to noise. (What a hack.)

		test_counter += 1
		if test_counter >= test_interval:
			test_counter = 0
			current_section = pts[(i - left_of_center):(i + right_of_center + 1)]
			if sum_distances(current_section, global_clusters) < \
				sum_distances(current_section, local_clusters):
				local_clusters = global_clusters

		# If the point we're to remove is the same as the point we're
		# to add, there's no need to do anything at all: just push
		# the same cluster.
		if (i - left_of_center > 0) and (i + right_of_center < len(pts)):
			if pts[i - left_of_center - 1] == pts[i + right_of_center]:
				local_clusters_per_point.append(local_clusters)
				continue

		# Remove the point at the left edge of the window, if it's a
		# real point (i.e. not the first window_size/2 times)
		if i - left_of_center > 0:
			classified_points.discard(pts[i - left_of_center - 1])

		# Add the point at the right edge of the window, if it's a
		# real point (i.e. not the last window_size/2 times)
		if i + right_of_center < len(pts):
			classified_points.add(pts[i + right_of_center])

			# If what we add is equal to any cluster's center, then every
			# cluster stays the same. (Only test the lowest cluster because
			# that one happens most often.)
			if pts[i+right_of_center] == local_clusters[0]:
				local_clusters_per_point.append(local_clusters)
				continue

		# Now all sorted points are on [i - left_of_center,
		# i + right_of_center + 1) unless truncated due to the beginning
		# or end of the pts array.

		# Get lower bounds for the range covered by each cluster but the
		# first, whose lower bound is -infinity.
		# https://stackoverflow.com/questions/23855976/
		lower = (local_clusters[1:] + local_clusters[:-1]) / 2

		# Get the sorted list indices for these, and for the first cluster
		# (which is 0). Also add the upper bound for the final cluster,
		# which is the number of points.
		lower_idx = np.array([0] + [classified_points.bisect_left(x) for x in lower] \
			+ [len(classified_points)])

		# Get the midpoints and do Lloyd's algorithm by finding new medians.
		# We need to do this explicitly because one of the bins may be empty
		# (e.g. no pulses in the corresponding band), which must then be
		# detected.
		new_local_clusters = np.array([0]*num_clusters)
		for j in range(num_clusters):
			# If we have too few points to make a proper decision,
			# then don't. Use the cluster from the last iteration to
			# avoid being skewed by outliers.
			if lower_idx[j+1] - lower_idx[j] <= 10:
				new_local_clusters[j] = local_clusters[j]
			else:
				median_idx = (lower_idx[j+1] + lower_idx[j] - 1) / 2
				# TODO: Handle medians that fall between points, here.
				# Or maybe this is good enough??
				new_local_clusters[j] = \
					classified_points[int(np.floor(median_idx))]

		local_clusters = new_local_clusters
		local_clusters_per_point.append(local_clusters)

	return np.array(local_clusters_per_point)

# The same, but with dynamic clusters as returned by local_k_median.
def assign_to_local_clusters(pulse_deltas, local_clusters_per_point):
	return assign_to_clusters(pulse_deltas, local_clusters_per_point.T)

# Returns a corrected pulse delay list given assignments and distances,
# i.e. what it would be like with absolutely no variation in local
# clusters. Mostly useful for fancy plots.
def correct_warp(assignments, residuals, global_clusters):
	return global_clusters[assignments] + residuals
