import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict, deque
import struct, timeit, copy
import dewarp
import kmedian
import kcenter
import re

from low_level import *
from mfm import *

# USAGE: Create the AU file with
# fluxengine convert fluxtoau -s SOURCE/:t=TRK:s=SIDE -o somename.au
# call get_pulse_deltas on somename.au
# and call various functions (see demonstrate).

# As an example, call demonstrate("tracks/MS_Plus_warped_track.au")

# TODO? Detect an abrupt end to a sector.
# Also be able to clock backwards.

# Call the groupings (categories) "bands" to follow KF terminology.

# TODO? Do something clever with the error report when decoding MFM.

# A typical run now is:
# pulses = get_pulse_deltas("file.au")
# preambles_clusterings = assign_clusters_to_preambles(pulses)
# assn = assign_on_preambles_clusterings(pulses, preambles_clusterings,
#	cluster_on_preambles(pulses))
# print_stats(decode_from_assignments(assn, datalen=512))

# and then your score is the "Number of sectors" count.

# This is considerably less format-agnostic than k-medians, but in turn
# is considerably better. However, I still can't beat fluxengine on
# every sample I have. In particular, fluxengine gets 18 good sectors on
# "low level format with noise", while I'm only able to get 17.

# preambles_clusterings doesn't contain all preambles. But then again,
# I'm not searching for C2 codes yet.

# VIGOR has C2 codes.

# Ideas:
#	- Multiple runs of the same sector: find the sector with the largest
#		distance to the first error. Align this with some other sector
#		that has an error before this to remove that error. Repeat until
#		all errors are gone.
#	- Fuzzy search for every possible IDAM for this track.


def scatter_plot(pulse_deltas, len, start=0, step=1):
	plt.scatter(range(0, len-start, step), pulse_deltas[start:len:step], s=10)
	plt.show()

# If the chunk is not correctly decoded (exhibits errors), try all k-center
# radii up to 10, filter those that are correct clusterings according
# to the preamble band intervals, and try to decode using those until
# something returns zero errors.

# Possible speedup: remember what worked in the past, do that first.
# Perhaps a bit of a hack...
radius_test_order = []

def kcenter_semi_brute(pulse_deltas, chunk_pos_preamble, pulse_len,
	data_length):

	global radius_test_order

	min_radius, max_radius = 1, 10
	start_pos = chunk_pos_preamble[0]

	band_intervals = get_band_intervals(classify_bands(
		pulse_deltas, start_pos,
		low_level.encode_to_assignment(chunk_pos_preamble[2])))

	candidate_clusterings = []
	with_previous_best_radius = None

	if radius_test_order == []:
		radius_test_order = list(range(min_radius, max_radius+1))

	# Try k-center, looking for an approach that gives zero errors.
	for radius_idx in range(len(radius_test_order)):
		clustering = kcenter.wt_optimal_k_center(
			pulse_deltas[start_pos:start_pos+pulse_len],
			len(band_intervals), radius_test_order[radius_idx])[0]

		if not valid_clustering(band_intervals, clustering):
			continue

		try:
			next_struct = decode_chunk(pulse_deltas, chunk_pos_preamble,
				pulse_len, data_length, override_clustering=clustering)
			if next_struct["_num_errors"] == 0:
				# That particular radius worked. Move it to the front of
				# the list so that it can be tested earlier next time around.
				radius_test_order.insert(0,
					radius_test_order.pop(radius_idx))
				return next_struct
		except KeyError:
			pass

	return None

# It's probably better to directly search for the out-of-band A1A1A1 tokens.
# But then I would have to preprocess and classify the whole pulse stream
# at once, first. Perhaps do that later.
def decode_all(pulse_deltas, data_length=0):
	print("Getting preamble positions...")
	preamble_pos_clusters = get_preamble_positions_clusters(pulse_deltas)

	decoded_structures = []
	data_length_needed = False

	for i in range(len(preamble_pos_clusters)):
		brute_corrected = False

		if data_length_needed and data_length != 0:
			print("Got data length, starting from the beginning.")
			return decode_all(pulse_deltas, data_length)

		start_pos = preamble_pos_clusters[i][0]
		if i == len(preamble_pos_clusters)-1:
			end_upper_bound = len(pulse_deltas)
		else:
			end_upper_bound = preamble_pos_clusters[i+1][0]

		print ("Decoding %d/%d: pos %d ... %d\t" % (i,
			len(preamble_pos_clusters), start_pos, end_upper_bound), end="")

		try:
			next_struct = decode_chunk(pulse_deltas, preamble_pos_clusters[i],
				end_upper_bound-start_pos, data_length)

			if next_struct["_num_errors"] > 0:
				given_length = next_struct["_end_pos"] - next_struct["_start_pos"]
				# The +10 to length is here because errors may make a structure
				# seem a little shorter or longer than it really is.
				corrected_next_struct =  kcenter_semi_brute(pulse_deltas,
					preamble_pos_clusters[i], given_length + 10,
					data_length)

				if corrected_next_struct is not None:
					next_struct = corrected_next_struct
					brute_corrected = True

			if "datalen" in next_struct:
				# Don't trust the data length if CRC isn't OK.
				if "CRC_OK" in next_struct and next_struct["CRC_OK"]:
					data_length = next_struct["datalen"]

			print("%s %s errors" % (next_struct["header"].rjust(13),
				str(next_struct["_num_errors"]).rjust(4)), end="")
			if next_struct["_num_errors"] > 0:
				print(" x")
			else:
				if brute_corrected:
					print (" C")
				else:
			 		print("")

			decoded_structures.append(next_struct)
		except KeyError:
			print("Could not distinguish the floppy chunk header!")
			pass
		except IndexError:
			if data_length == 0:
				print("Found DAM header without knowledge of its data length. Skipping.")
				data_length_needed = True
				pass

	return decoded_structures

# The coverage level tends to be around 90% for uncorrupted IBM MFM
# floppies with 512-byte sectors.
def get_coverage(pulse_deltas, decoded_structures):
	covered = sum([x["_end_pos"] - x["_start_pos"] for x in \
		decoded_structures])

	return covered/len(pulse_deltas)

# Group DAM data chunks by (track, head, sector) of the immediately
# preceding IDAM. More tricks may be required to fail safe if there's
# corruption (e.g. the structure as intended was IDAM1 DAM1 IDAM2 DAM2,
# but DAM1 and IDAM2 got severely corrupted, then DAM2 will be attributed
# to the sector in IDAM1). A simple check to see that the headers are
# consecutive should suffice. DONE, sort of. The 200 number is picked
# empirically.

def categorize_decoded(decoded_structures):
	by_sector = defaultdict(list)
	by_sector_valid = defaultdict(list)

	for i in range(1, len(decoded_structures)):
		this, prev = decoded_structures[i], decoded_structures[i-1]
		if this["header"] == "DAM" and prev["header"] == "IDAM":
			if "CRC" in this and prev["CRC_OK"]:
				ths = (prev["track"], prev["head"], prev["sector"])
				by_sector[ths].append((prev, this))
				gap_to_prev = this["_start_pos"] - prev["_end_pos"]
				if this["CRC_OK"] and gap_to_prev < 200:
					by_sector_valid[ths].append((prev, this))

	return by_sector, by_sector_valid

def get_invalid_chunk_stats(decoded):
	valid = len([ x for x in decoded if "CRC_OK" in x and x["CRC_OK"]])
	invalid = len([ x for x in decoded if "CRC_OK" in x and not x["CRC_OK"]])

	return (valid, invalid)

# "RETRY-77-4-t69.0.au" is a nice new standard.
# First you have to find the right clustering. Best so far is 22, 33, 44.
# Then dewarping helps immensely. But fluxengine gets 12 good sectors
# *without* dewarping. I'm not there yet!
# Now I get 15 without dewarping.

# Returns the number of MFM errors with this cluster assignment. Slow.
# Possible future improvement: only count errors up to the end of whatever
# struct is being investigated. That would require that decode inform the
# function about how many MFM sequence bits it consumed.
def get_tentative_errors(pulses, window_start, window_end, cluster_guess,
	use_pll = False, pll_alpha=0):
	assignments = []
	if use_pll:
		assignments = simplest_ever_pll(pulses[window_start:window_end], cluster_guess, pll_alpha)
		#assignments = shifted_window_pll(pulses, window_start, window_end, x=cluster_guess, alpha=pll_alpha)
	else:
		for i in range(window_start, window_end):
			assignment = np.argmin(np.abs(pulses[i] - cluster_guess))
			assignments.append(assignment)

	# ---------------------

	MFM_code, MFM_pos = assignments_to_MFM_code(np.array(assignments))
	struct_starts = get_struct_starts(assignments, MFM_code, MFM_pos)
	unexpected_errors = 0

	for i in range(len(struct_starts)):
		struct_start = struct_starts[i][1]
		if i == len(struct_starts)-1:
			struct_end = len(MFM_code)
		else:
			# We don't know how long the struct is; assume it goes
			# all the way to the next one. (May cause trouble with
			# sector-in-sector style copy protection.)
			struct_end = struct_starts[i+1][1]

		bits, errors = decode_from_MFM(MFM_code[struct_start:struct_end])
		# We expect one clock error for each repetition of the preamble.
		# Fewer than that means we have an unexpected error (none where
		# there should be one); more also means unexpected error.
		unexpected_errors += np.abs(len(errors)-3)

	return -len(struct_starts), unexpected_errors, -np.sum(assignments)

# Do a bunch of informed tests before going to brute force.
def baseline_tests(weights, i):
	# Quick and dirty. TODO: FIX LATER.
	try:
		if i == 0:
			return kmedian.wt_optimal_k_median_from_cumul(weights, 3)
		if i >= 1 and i < 8:
			return kmedian.optimal_k_median(kmedian.wt_optimal_k_median_from_cumul(weights, 3+i), 3)
		if i >= 8 and i < 16:
			return kcenter.wt_optimal_k_center_from_cumul(weights, 3, i-8)
	except IndexError:
		return np.array([20., 30., 40.])

# optimizing on short-77-4-t69 (full length) and then running demonstrate
# on the full version gets me
#	Without dewarping:
#		Valid CRC chunks: 2692	Invalid: 408	Fraction valid: 0.868
#		Number of sectors: 17	Bad sectors: 1
# fluxengine baseline is 12.

# Note that we *cannot* use CRC validity as an optimization criterion
# for this brute force cluster check. The reason is simple: that would
# destroy the only independent signal we have of success, and thus
# overfit easily. The 16 bit CRC would be overwhelmed in short order.

# (We could possibly use CRC validity on the sectors we've already
# recovered -- or really, accurate reproduction of the entire sector.)

# ((Every time you test against the CRC, you throw a die with a certain
# chance of success. If this die roll succeeds, the CRC return success
# even though the content is wrong.))
# The chance is 1/65k, so testing against the CRC 1024 times is on average
# going to return a false positive (CRC OK when it really isn't) 1.6% of
# the time such a batch is run. (Unless my reasoning is off.)
# Doing actual directed optimization is much more likely to return a false
# result.

def brute_force_cluster_check(pulses, window, known_good=np.array([20., 30., 40.]), tries=200):
	record_error = (np.inf, 0)
	record_clusters = None
	record_alpha = 0
	alpha = 0
	weights = kmedian.wt_get_prefix_sum(pulses[:window])
	use_pll = False

	for i in range(tries):
		if i == 0:
			if known_good is not None:
				clusters = np.array(known_good)
			else:
				continue

		if i < 17:
			cluster_guess = baseline_tests(weights, i-1)
			k = None
			b = None
		else:
			k = np.random.uniform(8, 11)
			b = np.random.uniform(-5, 5)
			clusters = np.array([2*k + b, 3*k + b, 4*k + b])
			use_pll = np.random.choice([True, False])
			use_pll = False # doesn't properly work yet
			if use_pll:
				alpha = np.random.uniform(0.01, np.sqrt(0.12))**2
			else:
				alpha = 0

		errors = get_tentative_errors(pulses, 0, window, clusters, use_pll,
			alpha)
		tiebreak = (errors == record_error) and \
			np.sum(np.abs(clusters - known_good)) < \
			np.sum(np.abs(record_clusters - known_good))

		if errors < record_error or tiebreak:
			print(clusters, k, b, errors)
			record_error = errors
			record_clusters = clusters
			record_alpha = alpha

	return record_clusters, record_alpha, record_error

def do_brute_force(pulses, window_size=2000):
	assignments = []
	last_clustering = [22., 33., 44.]

	for i in range(window_size, len(pulses)-window_size, window_size):
		print("Testing %d of %d, i.e. %.3f %%" %(i, len(pulses),
			(i*100)/len(pulses)))
		cluster_guess, pll_alpha, error = brute_force_cluster_check(
			pulses[i-window_size:i+window_size], window=4000,
			known_good=last_clustering)
		last_clustering = cluster_guess
		print(cluster_guess)
		if alpha > 0:
			assignments += simplest_ever_pll(pulses[i:i+window_size],
				cluster_guess, pll_alpha)
		else:
			for j in range(window_size):
				assignment = np.argmin(np.abs(pulses[i+j] - cluster_guess))
				assignments.append(assignment)

	return assignments

# Some model-based stuff follow. This might be useful for dewarping,
# because every point informs the state space, thus avoiding the usual
# "loss of gradient" problem.

# The model is that
# Points that belong to the first band are closest to 2*k + b
# Points that belong to the second band are closest to 3*k + b
# Points that belong to the third band are closest to 4*k + b.

# (Am I making things too hard? Could all of this be replaced with a
# k-medians or k-center clustering around only these points? That would
# perhaps be an interesting variant.)

# This returns 0 if the given assignment k,b, and given points
# assignment, and otherwise returns the number of points violating
# it as the integer part, and the distance of the point furthest from
# being satisfied as the fractional part.
def test_mfm(points_by_band, k, b):
	bands_concatenated = np.concatenate(points_by_band)
	centers = np.array([2, 3, 4]) * k + b

	midpoints = (centers[1:] + centers[:-1]) / 2.0
	# Add infinity points on each side because it's never possible
	# to get low enough to not be assigned to the first band, or high
	# enough not to be assigned to the last.
	midpoints = np.concatenate([[-np.inf], midpoints, [np.inf]])

	wrong_points = 0
	max_distance = 0

	for band_idx in range(len(points_by_band)):
		for point in points_by_band[band_idx]:
			if point < midpoints[band_idx]:
				cur_distance = midpoints[band_idx] - point
				max_distance = max(max_distance, cur_distance)
				wrong_points += 1
			if point >= midpoints[band_idx+1]:
				cur_distance = point - midpoints[band_idx+1]
				max_distance = max(max_distance, cur_distance)
				wrong_points += 1

	return wrong_points + max_distance / np.max(bands_concatenated)

# Quick and dirty. There's probably something more mathematically
# elegant to be had...
def do_segmented_brute_force(pulses):
	assignments = []
	last_clustering = [22., 33., 44.]
	# Start 60 points earlier so we can get the zero train before the
	# actual preamble, too.
	promising_starts = np.array(get_potential_preambles(pulses))-60

	for i in range(len(promising_starts)):
		print("Testing %d of %d, i.e. %.3f %%" %(i, len(promising_starts),
			(i*100)/len(promising_starts)))
		start_pos = promising_starts[i]
		if i == len(promising_starts)-1:
			end_pos = len(pulses)
		else:
			end_pos = promising_starts[i+1]

		print("Interval: %d to %d (%d pulses)" % (start_pos, end_pos, end_pos-start_pos))

		# Skip very short intervals
		if end_pos - start_pos < 60:
			cluster_guess = last_clustering
		else:
			cluster_guess, pll_alpha, error = brute_force_cluster_check(
				pulses[start_pos:], window=end_pos-start_pos,
				known_good=last_clustering)
		# There should be at least one preamble struct here, so don't trust
		# the clustering if it reports none
		if error[0] != 0:
			last_clustering = cluster_guess
		print(last_clustering)
		if pll_alpha > 0:
			assignments += simplest_ever_pll(pulses[start_pos:end_pos],
				last_clustering, pll_alpha)
		else:
			for j in range(start_pos, end_pos):
				assignment = np.argmin(np.abs(pulses[j] - last_clustering))
				assignments.append(assignment)

	return assignments

# Semi-bruting from a boundary.

# Simple way of sifting clusters that can't even recognize the
# preamble. Uses intervals, and the constraints above.
def can_cluster_work(intervals, candidate_cluster):
	if candidate_cluster[0] > intervals[0][1]: return False
	if candidate_cluster[1] < 2 * intervals[0][1] - candidate_cluster[0]: return False
	if candidate_cluster[1] > 2 * intervals[1][0] - candidate_cluster[0]: return False
	if candidate_cluster[2] < 2 * intervals[1][1] - candidate_cluster[1]: return False
	if candidate_cluster[2] > 2 * intervals[2][0] - candidate_cluster[1]: return False

	return True

def semi_brute_cluster(pulses, partial_preamble_start, next_preamble_start, tries=250):
	record_error = (np.inf, 0)
	record_clusters = None
	record_alpha = 0
	alpha = 0
	use_pll = False

	points_by_band = categorize_extended_preamble(pulses, partial_preamble_start)
	known_good = np.array(cluster_on_preambles(pulses))
	if points_by_band is None:
		return None

	coverage_intervals = get_cluster_intervals(points_by_band)

	for i in range(tries):
		clusters = random_cluster(coverage_intervals)
		#use_pll = np.random.choice([True, False])
		use_pll = False
		if use_pll:
			alpha = np.random.uniform(0.01, np.sqrt(0.15))**2
		else:
			alpha = 0

		errors = get_tentative_errors(pulses,
			partial_preamble_start-additional_pulses_full_preamble(),
			next_preamble_start, clusters, use_pll,
			alpha)
		tiebreak = (errors == record_error) and \
			np.sum(np.abs(clusters - known_good)) < \
			np.sum(np.abs(record_clusters - known_good))

		if errors < record_error or tiebreak:
			print(clusters, use_pll, errors)
			record_error = errors
			record_clusters = clusters
			record_alpha = alpha

	return record_clusters, record_alpha, record_error

# I can't believe it's not a hack.
# TODO: I really should do something about that -70...
def semi_brute_per_sector(pulses, tries=500):
	gpp = get_pruned_potential_preambles(pulses)
	assignments = []
	for i in range(len(gpp)-1):
		print("Testing %d of %d, i.e. %.3f %%" %(i, len(gpp)-1,
			(i*100)/(len(gpp)-1)))
		clusters, alpha, error = semi_brute_cluster(
			pulses, gpp[i], gpp[i+1], tries=tries)
		# Cute hack where alpha=0 just does a regular assignment. Heh!
		assignments += simplest_ever_pll(
			pulses[gpp[i]-additional_pulses_full_preamble():\
				gpp[i+1]-additional_pulses_full_preamble()],
			clusters, alpha)
	return assignments

# Given some pulses and pruned preambles-clusterings, determine which
# decode without errors.
# pruned preamble-clusterings is a list of pairs of pruned preambles and
# clusterings to use for that location to the next, or None if unknown.
# The function determines for which of the Nones the current clustering
# is appropriate, and then replaces the None thus, returning the copy
# that has been such augmented.
# Something is wrong here... it's being a lot more picky than
# decode_from_assignments is.
def get_chunks_without_errors(pulses, preambles_clusterings,
	clustering, dynamic_clustering = False):

	if (clustering is not None) and len(clustering) < 3:
		raise Exception("Invalid clustering provided")

	out_preambles = np.copy(preambles_clusterings)
	correct, seen = 0, 0

	seen_clusters = []
	coverage_intervals = []

	# Generate intervals for early screening.
	for i in range(len(preambles_clusterings)):
		points_by_band = categorize_extended_preamble(pulses, preambles_clusterings[i][0])
		coverage_intervals.append(get_cluster_intervals(points_by_band))

	for i in range(len(preambles_clusterings)):
		# already classified
		if preambles_clusterings[i][1] is not None:
			continue

		seen += 1

		if not (dynamic_clustering or can_cluster_work(coverage_intervals[i], clustering)):
			continue

		start_pos = preambles_clusterings[i][0] - additional_pulses_full_preamble()
		if i == len(preambles_clusterings)-1:
			end_pos = len(pulses)
		else:
			end_pos = preambles_clusterings[i+1][0] - additional_pulses_full_preamble()

		this_clustering = clustering

		if dynamic_clustering:
			this_clustering = cluster_on_preambles(pulses,
				[preambles_clusterings[i][0]])

		structs, errors, bitsum = get_tentative_errors(pulses,
			start_pos, end_pos, this_clustering)

		# We would expect 3 errors due to the intentional preamble errors.
		if errors == 0 and structs <= -1:
			out_preambles[i][1] = clustering
			print("Preamble %d/%d (%d) classified as without error" % (i, len(preambles_clusterings), preambles_clusterings[i][0]))
			if dynamic_clustering:
				seen_clusters.append(this_clustering)
			correct += 1

	print("Classified %d/%d(%.2f%%) without error" % (correct, seen, correct*100/seen))

	# https://stackoverflow.com/questions/3724551
	if dynamic_clustering:
		seen_clusters = [list(x) for x in set(tuple(x) for x in seen_clusters)]
		print("Seen dynamic clusters", seen_clusters)

	return out_preambles, seen_clusters

def assign_on_preambles_clusterings(pulses, preambles_clusterings,
	default_clustering):

	end_pos = preambles_clusterings[0][0]# - additional_pulses_full_preamble()
	assignments = assign(pulses[:end_pos], default_clustering)
	clustering = default_clustering

	for i in range(len(preambles_clusterings)):
		start_pos = preambles_clusterings[i][0] - additional_pulses_full_preamble()
		if i == len(preambles_clusterings)-1:
			end_pos = len(pulses)
		else:
			end_pos = preambles_clusterings[i+1][0] - additional_pulses_full_preamble()

		#clustering = preambles_clusterings[i][1]
		#if clustering is None:
		#	clustering = default_clustering
		if preambles_clusterings[i][1] is not None:
			clustering = preambles_clusterings[i][1]

		assignments += assign(pulses[start_pos:end_pos], clustering)

	return assignments

def assign_clusters_to_preambles(pulses, brute_force=False):
	preambles = get_pruned_potential_preambles(pulses)
	preambles_clusterings = [[x, None] for x in preambles]
	spotted_clusters = []

	print("ACTP: Trimmed mean")
	# First try the 0.25 trimmed mean.
	dc = np.array(cluster_on_preambles(pulses, trim=0.25))
	preambles_clusterings, spotted_clusters = get_chunks_without_errors(
		pulses, preambles_clusterings, dc)

	print("ACTP: Median")
	# Then try the median
	dc = np.array(cluster_on_preambles(pulses))
	preambles_clusterings, spotted_clusters = get_chunks_without_errors(
		pulses, preambles_clusterings, dc)

	print("ACTP: k-medians")
	# Then try k-medians
	dc = kmedian.wt_optimal_k_median(pulses, 3)
	preambles_clusterings, spotted_clusters = get_chunks_without_errors(
		pulses, preambles_clusterings, dc)

	print("ACTP: Dynamic discovery")
	# Then dynamic discovery.
	preambles_clusterings, spotted_clusters = get_chunks_without_errors(
		pulses, preambles_clusterings, None, True)
	for i in spotted_clusters:
		preambles_clusterings, x = get_chunks_without_errors(
			pulses, preambles_clusterings, i)

	if not brute_force:
		return preambles_clusterings

	print("ACTP: Brute force")
	# Then brute force.
	for i in range(len(preambles_clusterings)-1):
		print("\t%d of %d (%.2f%%)" % (i, len(preambles_clusterings)-1,
			100*i/(len(preambles_clusterings)-1)))

		if preambles_clusterings[i][1] is not None:
			continue

		start_pos = preambles_clusterings[i][0] - additional_pulses_full_preamble()
		end_pos = preambles_clusterings[i+1][0] - additional_pulses_full_preamble()

		record_clusters, record_alpha, record_error = \
			semi_brute_cluster(pulses,preambles_clusterings[i][0], preambles_clusterings[i+1][0])
		print("Output", record_clusters, record_alpha, record_error)
		if record_error[0] < 0 and record_error[1] == 0:
			print("Maybe something here")
			preambles_clusterings, x = get_chunks_without_errors(
				pulses, preambles_clusterings, record_clusters)

	return preambles_clusterings


def print_stats(decoded_structure):
	invstat = get_invalid_chunk_stats(decoded_structure)
	categorized, categorized_valid = categorize_decoded(decoded_structure)

	try:
		print("\t\tValid CRC chunks: %d\tInvalid: %d\tFraction valid: %.3f" %
			(invstat[0], invstat[1], invstat[0]/(invstat[0]+invstat[1])))
		good_sector_keys = set(categorized_valid.keys())
		print("\t\tNumber of sectors: %d\tGood sectors: %d" % (len(categorized.keys()), len(good_sector_keys)))
	except ZeroDivisionError:
		print("\t\tNo chunks found, neither valid nor invalid.")

def decode_from_assignments(assignments, datalen=0):
	MFM, MFM_pos = assignments_to_MFM_code(assignments)
	struct_starts = get_struct_starts(assignments, MFM, MFM_pos)
	decoded = decode_all_from_MFM(assignments, MFM, MFM_pos,
		struct_starts, datalen=datalen)

	return decoded

def demonstrate(au_name):
	pulses = get_pulse_deltas(au_name)
	decoded = decode_all(pulses)

	print_stats(decoded)
	if len(decoded) > 0:
		print("\n\t\t%d%% of the floppy data was decoded." %
			(100*get_coverage(pulses, decoded)))

# hist_cat, hist_cat_valid = categorize_decoded(hist_decoded)

# TODO: Coverage function (what we recovered vs the entire floppy)

# hist. Good track 53 sector 8 data chunk: start 1497875
# hist. Bad track 53 sector 8 data chunk: start 1759015

# Seems to be overfitting of a stats gradient; i.e. a delay that's right
# on the border between two categories pulls the wrong category in its
# direction.

# start_indices = [x[1]["_start_idx"] for x in hist_cat[(53, 1, 13)]]
# hist_MFM, hist_MFM_pos = assignments_to_MFM_code(hist_assign)

# Error types:
# Missing flux transition		sees RNNNR, could be RNNRNR or RNRNNR
# Misclassification/close call	sees RNNR, could be RNR or RNNNR

# Extraneous flux transitions are very uncommon, but missing ones are
# more common.

#pulses_77 = get_pulse_deltas("RETRY-77-4-t77.0.hard.au")
#v77 = decode_all(pulses_77)

#demonstrate("warped.au")