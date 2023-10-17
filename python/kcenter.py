import numpy as np
import sortedcontainers as sc

# HACK HACK, attempt to implement k-center.
# The optimal point for a k-center cluster is still the median, since
# max norm is an Lp norm. However, the cost function is different, and
# the cost C_i is

# C_i(m, j) = max(D_previous..., CC(...))
# with CC being max(mu - x[i], x[j-1] - mu).

# With Lagrangian outlier censoring, every difference from the median
# greater than the specified radius gets replaced by that radius. Tune the
# radius to remove as many or few points as you'd like: lower radius means
# more points are removed.

# It's very slow, and the stuff below is very cut and paste.
# Now it's not as slow, but unweighted is slower than it used to be.
# There's a bug somewhere involving tie breaking, with weighted providing
# a different result than unweighted. However, both results have the same
# score (number of dropped points and maximum distance).

# Get the prefix sum array x (for use in CC etc). Points are given in
# ascending order.
def get_cdf(points):
	return np.cumsum(points)

# Grønlund et al.'s CC function, adapted to the k-center problem with
# outliers.
# CC(x,i,j) is the cost of grouping x_i...x_j into one cluster with the
# optimal cluster point. By programming convention, the interval is
# half-open and indexed from 0, unlike the paper's convention of closed
# intervals indexed from 1.
# Since the k-center problem with outliers allows for discarding some
# number of points, the optimal point is no longer the median, because it
# may be the case that another point allows for discarding fewer points,
# even if the maximum value of the residuals is worse. Thus we need an
# auxiliary function to determine the optimum center, which is done in p
# log p time for p points.
# It can be done quicker if required (definitely in p time, possibly in
# log(p) time). However, that would come at a cost of additional code
# complexity, so I'll leave it alone for the time being.
# x must be sorted in ascending order.
# Returns the location of the center, the maximum residual (distance from
# center), and the number of dropped points.
def find_optimal_center(x, i, j):
	if i >= j:
		raise Exception("find_optimal_center: empty area")

	recordholder = None
	record_distance = 0

	# Assume the left end of the interval captured by the center is k.
	# The other end is wherever the value is greater than the value at k,
	# by more than 2 * radius -- or the end of the range if it doesn't
	# happen.
	for k in range(i, j):
		# No future record can break the current one
		if j-k < record_distance:
			continue

		# First element exceeding this value.
		k_end = k + x[k:j].searchsorted(x[k] + 2 * radius, side='right')
		if (k_end - k) > record_distance:
			record_distance = k_end - k
			recordholder = k

	if recordholder is None:
		raise Exception("Could not find recordholder - bug in code")

	first, last = recordholder, recordholder + record_distance - 1
	center = (x[first] + x[last]) / 2.0
	max_residual = max(center-x[i], x[j-1]-center)
	max_residual = min(radius, max_residual)

	dropped_values = j-i-record_distance

	return center, max_residual, dropped_values

def CC(x, px, i, j):
	if i >= j: return (0,0)

	m = (j+i-1)/2.0

	# Lots of nasty off-by-one traps here.
	mu = (x[int(np.floor(m))] + x[int(np.ceil(m))]) / 2		# median

	# For k-center with outliers, there're two options: either some
	# of the points between i and j will be classified as outliers,
	# or none will be. If none, then the median is the optimal position
	# for the center and we're good to go. But if either end of the
	# range exceeds the radius (threshold for considering something an
	# outlier), with the median as the cluster, then there may be a
	# better point and we need to manually find it.

	if mu - x[i] > radius or x[j-1] - mu > radius:
		center, maxval, dropped_points = find_optimal_center(x, i, j)
	else:
		center, maxval, dropped_points = mu, \
			max(mu - x[i], x[j-1] - mu), 0

	return maxval, dropped_points

# the stuff from the kn log n algorithm. I no longer know its exact
# complexity, due to find_optimal_center.
# D_previous is the D vector for (i-1) clusters, or empty if i < 2.
# It's possible to do this even faster (and more incomprehensibly).
# See Grønlund et al. for that.

# Calculate C_i[m][j] given D_previous = D[i-1]. p. 4
radius = 30000

def C_i(x, px, i, D_previous, m, j):
	if i == 1:
		return CC(x, px, 0, m)
	else:
		prev = D_previous[min(j, m)]
		this = CC(x, px, j, m)
		return (max(prev[0], this[0]), (prev[1] + this[1]))

# Find the (location of the) minimum value of each row of an nxn matrix in
# n log n time, given that the matrix is monotone (by the definition of
# Grønlund et al.) and defined by a function M that takes row and column
# as parameters. All boundary values are closed, [min_row, max_row]
# p. 197 of
# AGGARWAL, Alok, et al. Geometric applications of a matrix-searching
# algorithm. Algorithmica, 1987, 2.1-4: 195-208.
def monotone_matrix_indices(M, min_row, max_row, min_col, max_col,
	Tout, Dout):
	if max_row <= min_row: return

	# PERF: Only runs M() once per test, thus speeding up the process
	# by a factor of two.
	def M_combined_score(m, j):
		out = M(m, j)
		return out[0] + 1e-3 * out[1]

	cur_row = int(min_row + np.floor((max_row-min_row)/2.0))

	# Find the minimum column for this row.
	# PERF: Improve performance by noting that j > m is redundant.
	min_here = min_col + np.argmin([(M_combined_score(cur_row, j)) for j in \
			range(min_col, min(max_col, cur_row)+1)])

	if Tout[cur_row] != 0:
		raise Exception("tried to set a variable already set.")

	Tout[cur_row] = min_here
	Dout[cur_row] = M(cur_row, min_here)

	# No lower row can have a minimum to the right of the current minimum.
	# Recurse on that assumption.
	monotone_matrix_indices(M, min_row, cur_row, min_col, min_here+1, Tout, Dout)
	# And no higher row can have a minimum to the left.
	monotone_matrix_indices(M, cur_row+1, max_row, min_here, max_col, Tout, Dout)

def optimal_k_center(points, num_clusters, radius_in):
	pts_sorted = np.array(sorted(points))
	pts_cumulative = get_cdf(pts_sorted)
	pts_len = len(pts_sorted)

	# Ugly hack since I can't yet be bothered to thread the radius
	# parameter through all the functions needed... Fix later!
	global radius
	radius = radius_in

	T, D = [[]], [[]]

	for i in range(1, num_clusters+1):
		# If we have achieved an optimal clustering with fewer clusters,
		# our work is done.
		if len(D[-1]) > 0 and D[-1][-1] == 0:
			continue

		M = lambda m, j: C_i(pts_sorted, pts_cumulative, i, D[-1], m, j)

		Tout = [0] * (pts_len+1)
		Dout = [0.0] * (pts_len+1)

		#print(len(Tout))

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
	num_dropped = D[num_clusters][cur_clustering_range][1]

	for i in range(num_clusters, 0, -1):
		cluster_boundaries.append(T[i][cur_clustering_range])
		new_cluster_range = T[i][cur_clustering_range]
		center, maxval, dropped_points = find_optimal_center(pts_sorted,
				new_cluster_range, cur_clustering_range)
		centers.append(center)
		cur_clustering_range = new_cluster_range

	return np.array(sorted(centers)), num_dropped

# --- WEIGHTED STUFF ---

# The points all have one of a very small range of values (compared to the
# number of points). Since k-center clustering is either a point or the
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

# K-center note: I'm still not entirely sure of the new complexity.
# It's probably not kn log n.
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

	lower_median_idx = i + weighted_cdf[i:,0].searchsorted(
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

# Returns the number of points represented by the weighted cdf interval
# [i, j)
def wt_num_points_between(weighted_cdf, i, j):
	j = min(j, len(weighted_cdf))
	i = max(i, 0)

	if i >= j:
		return 0

	points_below = 0
	if i > 0:
		points_below = weighted_cdf[i-1][0]

	return weighted_cdf[j-1][0] - points_below

# Note: i and j are indices onto weighted_cdf. So e.g. if i = 0, j = 2 and
# weighted cdf is [[2, 0, 0], [4, 2, 1], [5, 2, 2]], then that is the cost
# of clustering all points between the one described by the zeroth weighted
# cdf entry, and up to (but not including) the last. In other words, it's
# CC([0, 0, 1, 1], 0, 5).

def wt_find_optimal_center(weighted_cdf, i, j):
	if i >= j:
		raise Exception("wt_find_optimal_center: empty area")

	recordholder = None
	record_distance = 0
	record_dist_index = 0

	for k in range(i, j):
		# No future record can break the current one
		if wt_num_points_between(weighted_cdf, k, j) < record_distance:
			continue

		# First element exceeding this value.
		k_end = k + weighted_cdf[k:j,2].searchsorted(weighted_cdf[k][2] + 2 * radius, side='right')
		if wt_num_points_between(weighted_cdf, k, k_end) > record_distance:
			record_distance = wt_num_points_between(weighted_cdf, k, k_end)
			record_dist_index = k_end - k
			recordholder = k

	if recordholder is None:
		raise Exception("Could not find recordholder - bug in code")

	first, last = recordholder, recordholder + record_dist_index - 1
	center = (weighted_cdf[first][2] + weighted_cdf[last][2]) / 2.0
	max_residual = max(center-weighted_cdf[i][2], weighted_cdf[j-1][2]-center)
	max_residual = min(radius, max_residual)

	dropped_values = wt_num_points_between(weighted_cdf, i, j) - record_distance
	return center, max_residual, dropped_values

def wt_CC(weighted_cdf, i, j):
	if i >= j: return (0,0)

	# Before we go about calculating centers, determine if any points
	# will be dropped. If so, just go directly to the function that
	# determines which those are.

	if weighted_cdf[j-1][2] - weighted_cdf[i][2] > 2 * radius:
		center, maxval, dropped_points = wt_find_optimal_center(weighted_cdf,
			i, j)
		return maxval, dropped_points

	# Hold on, isn't the midrange better? Consider e.g. the set
	# [0, 0, 0, 100, 100]. The mid-range has max residual 50, but the
	# median has max residual 100. Whoops.

	# With leximin, it would be even trickier. Good thing we don't have
	# to deal with that here...

	center = (weighted_cdf[i][2] + weighted_cdf[j-1][2]) / 2.0
	maxval = weighted_cdf[j-1][2] - center
	dropped_points = 0

	return maxval, dropped_points

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
		prev = D_previous[min(j, m)]
		this = wt_CC(weighted_cdf, j, m)
		return (max(prev[0], this[0]), (prev[1] + this[1]))

# for the floppy images, 2 is a good choice for radius. 1 is too little,
# but even something like 10 will work.
def wt_optimal_k_center_from_cumul(wt_pts_cumulative, num_clusters, radius_in):

	unique_points = len(wt_pts_cumulative)

	# Ugly hack since I can't yet be bothered to thread the radius
	# parameter through all the functions needed... Fix later!
	global radius
	radius = radius_in

	T, D = [[]], [[]]

	for i in range(1, num_clusters+1):
		# If we have achieved an optimal clustering with fewer clusters,
		# our work is done.
		if len(D[-1]) > 0 and D[-1][-1] == 0:
			continue

		# Otherwise, find the optimal T (cluster boundaries) and D
		# (cluster costs) using the Aggarwal algorithm.
		M = lambda m, j: wt_C_i(wt_pts_cumulative, i, D[-1], m, j)

		Tout = [0] * (unique_points+1)
		Dout = [(0.0, 0.0)] * (unique_points+1)

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
	num_dropped = D[num_clusters][cur_clustering_range][1]

	for i in range(num_clusters, 0, -1):
		cluster_boundaries.append(T[i][cur_clustering_range])
		new_cluster_range = T[i][cur_clustering_range]
		# Reconstruct the cluster that's the median point between
		# new_cluster_range and cur_clustering_range.
		center, maxval, dropped_points = wt_find_optimal_center(
			wt_pts_cumulative, new_cluster_range, cur_clustering_range)
		centers.append(center)
		cur_clustering_range = new_cluster_range

	return np.array(sorted(centers)), num_dropped

def wt_optimal_k_center(points, num_clusters, radius_in):
	return wt_optimal_k_center_from_cumul(wt_get_prefix_sum(points),
		num_clusters, radius_in)

def wt_debug(points, wt_pts_cumulative=None):
	if wt_pts_cumulative is None:
		wt_pts_cumulative = wt_get_prefix_sum(points)

	for i in range(len(wt_pts_cumulative)+1):
		print("")
		for j in range(len(wt_pts_cumulative)+1):
			print("Debug: CC", i, ",", j, "=", wt_CC(wt_pts_cumulative, i, j))

def get_residuals(data, clustering):
	residuals = np.abs(data - clustering[0])
	for i in range(1, len(clustering)):
		residuals = np.minimum(residuals, np.abs(data - clustering[i]))

	return residuals
