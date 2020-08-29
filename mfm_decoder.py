import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict, deque
import struct, timeit, copy
import dewarp
import kmedian
import kcenter
import re

import mfm, low_level
import approx_search
from low_level import get_pulse_deltas

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
	data_length, assignment_function=low_level.assign_partial):

	global radius_test_order

	min_radius, max_radius = 1, 8
	start_pos = chunk_pos_preamble[0]

	band_intervals = low_level.get_band_intervals(low_level.classify_bands(
		pulse_deltas, start_pos,
		low_level.encode_to_assignment(chunk_pos_preamble[2])))

	candidate_clusterings = []
	with_previous_best_radius = None

	# -1 means "use k-median".
	if radius_test_order == []:
		radius_test_order = [-1] + list(range(min_radius, max_radius+1))

	# Try k-center or k-median, looking for an approach that gives zero
	# errors.
	for radius_idx in range(len(radius_test_order)):
		if radius_test_order[radius_idx] == -1:
			clustering = kmedian.wt_optimal_k_median(
				pulse_deltas[start_pos:start_pos+pulse_len],
				len(band_intervals))
		else:
			clustering = kcenter.wt_optimal_k_center(
				pulse_deltas[start_pos:start_pos+pulse_len],
				len(band_intervals), radius_test_order[radius_idx])[0]

		if not low_level.valid_clustering(band_intervals, clustering):
			continue

		try:
			next_struct = mfm.decode_chunk(pulse_deltas, chunk_pos_preamble,
				pulse_len, data_length, override_clustering=clustering,
				assignment_function=assignment_function)
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
	preamble_pos_clusters = mfm.get_preamble_positions_clusters(pulse_deltas)

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
			next_struct = mfm.decode_chunk(pulse_deltas, preamble_pos_clusters[i],
				end_upper_bound-start_pos, data_length)

			if next_struct["_num_errors"] > 0:
				given_length = next_struct["_end_pos"] - next_struct["_start_pos"]
				# The +30 to length is here because errors may make a structure
				# seem a little shorter or longer than it really is.
				#corrected_next_struct =  kcenter_semi_brute(pulse_deltas,
				#	preamble_pos_clusters[i], given_length + 30,
				#	data_length)

				#if corrected_next_struct is not None:
				#	next_struct = corrected_next_struct
				#	brute_corrected = True

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
			# Insert a placeholder structure so that later bruting
			# can fix it, if possible.
			placeholder = {}
			placeholder["header"] = "Unknown"
			placeholder["CRC_OK"] = False
			placeholder["_start_pos"] = start_pos
			placeholder["_end_pos"] = end_upper_bound
			decoded_structures.append(placeholder)
			pass
		except IndexError:
			if data_length == 0:
				print("Found DAM header without knowledge of its data length. Skipping.")
				data_length_needed = True
				pass

	return decoded_structures

# EXPERIMENTAL
# Known bugs: VIGOR (low_level_format_with_noise.flux) produces a :-X error
# which should never happen (because that signals no struct was found at all,
# but since augmented_full_brute only improves positions that are known
# to hold structs, it should always find at least a corrupted struct there.)

# USAGE:
# pulses = get_pulse_deltas(...)
# decoded = decode_all(pulses)
# augmented_decoded = augment_full_brute(pulses, decoded, 512) #assuming 512-byte data chunks.

def full_brute(pulse_deltas, chunk_pos_preamble, end_upper_bound,
	data_length=0):

	preamble_pos, clustering, preamble, short_preamble = chunk_pos_preamble

	classified_points = low_level.classify_bands(pulse_deltas, preamble_pos,
		low_level.encode_to_assignment(preamble), len(clustering))

	band_intervals = low_level.get_band_intervals(classified_points)

	almost_all_bands = get_all_bands(pulse_deltas[preamble_pos:end_upper_bound],
		band_intervals, 100)

	error_record, recordholder, recordholder_alpha = np.inf, [], 0
	record_struct = None

	for i in range(len(almost_all_bands)):
		if error_record == 0 and "CRC_OK" in record_struct and record_struct["CRC_OK"]:
			continue

		clustering = low_level.band_intervals_to_clusters(np.array(almost_all_bands[i]))
		alpha = 0
		use_alpha = False
		#alpha = 0.02 + np.sqrt(np.random.uniform(0, 0.20**2))
		#use_alpha = True

		if use_alpha:
			assn = lambda a, b, c, d: dewarp.simplest_ever_pll(a, b, c, d,
				alpha)
		else:
			alpha = 0
			assn = low_level.assign_partial

		try:
			next_struct = mfm.decode_chunk(pulse_deltas, chunk_pos_preamble,
					end_upper_bound-preamble_pos, data_length,
					override_clustering=clustering,
					assignment_function=assn)

			# Permit checking of CRC for the first 50 or so tests.
			# In general, we should never check for CRC validity
			# because each check produces a slight possibility of a
			# false positive, and we can't tell those apart.
			# But as long as we only do it for a very small number,
			# it *should* be okay...
			override = next_struct["_num_errors"] == 0 \
			 	and "CRC_OK" in next_struct \
			 	and next_struct["CRC_OK"] and i < 50

			if next_struct["_num_errors"] < error_record or override:
				error_record = next_struct["_num_errors"]
				recordholder = clustering
				recordholder_alpha = alpha
				record_struct = next_struct
				print("\t", i, ":", error_record, recordholder, recordholder_alpha)
		except IndexError:
			pass
		except KeyError:
			pass

	return record_struct

def augment_full_brute(pulse_deltas, decoded_so_far, data_length=0):
	preamble_pos_clusters = mfm.get_preamble_positions_clusters(pulse_deltas)
	augmented_decoded = []

	for i in range(len(decoded_so_far)):
		if "_num_errors" in decoded_so_far[i] and decoded_so_far[i]["_num_errors"] == 0:
			augmented_decoded.append(decoded_so_far[i])
			continue

		# This happens if the structure is unknown and we don't know
		# how many errors it has (because its length depends on its
		# type). In that case, just pretend it has infinite errors.
		if not "_num_errors" in decoded_so_far[i]:
			decoded_so_far[i]["_num_errors"] = np.inf

		start_pos = preamble_pos_clusters[i][0]
		if i == len(preamble_pos_clusters)-1:
			end_upper_bound = len(pulse_deltas)
		else:
			end_upper_bound = preamble_pos_clusters[i+1][0]

		print ("Bruting %d/%d: pos %d ... %d\t" % (i,
			len(preamble_pos_clusters), start_pos, end_upper_bound))

		possible_improved_struct = full_brute(pulse_deltas,
			preamble_pos_clusters[i], end_upper_bound+40, data_length)

		if possible_improved_struct is None:
			print(":-X could not find a single struct. This should never happen and indicates a bug, but continuing anyway...")
			augmented_decoded.append(decoded_so_far[i])
			continue

		# See above.
		if not "_num_errors" in possible_improved_struct:
			possible_improved_struct["_num_errors"] = np.inf

		if possible_improved_struct["_num_errors"] < decoded_so_far[i]["_num_errors"]:
			if possible_improved_struct["_num_errors"] == 0 and possible_improved_struct["CRC_OK"]:
				print(":-) fixed")
			else:
				print(":-/ reduced errors to", possible_improved_struct["_num_errors"])
			augmented_decoded.append(possible_improved_struct)
		else:
			print(":-( no improvement to error count")
			augmented_decoded.append(decoded_so_far[i])

	return augmented_decoded

# Additional code for the exhaustive band enumeration strategy follows.
# TODO? Make it a dynamic programming algorithm to speed up matters,
# filtering away anything that can produce an error as we go.

# Also don't judge the strategy by its speed as implemented in Python.
# I haven't optimized the algorithms: get_all_bands in particular is very
# slow as it is.

# Returns a tuple giving the bands this point might belong to.
def get_band(point_value, intervals):
	for i in range(len(intervals)):
		if point_value >= intervals[i][0] and point_value <= intervals[i][1]:
			return (i,)

	if point_value < intervals[0][0]:
		return (0,)

	for i in range(1, len(intervals)):
		if point_value > intervals[i-1][1] and point_value < intervals[i][0]:
			return(i-1, i)

	return (len(intervals)-1,)

# Outputs every possible different set of band intervals that the
# pulse data can be classified into, for semi-brute force purposes.
# If max_num is set to a finite value, it will abort early if the list
# is that full.
def get_all_bands(pulses, intervals, max_num=np.inf):
	candidate_configurations = set([tuple(map(tuple, intervals))])
	candidate_configurations_next = set()

	for i in range(len(pulses)):
		if len(candidate_configurations) > max_num:
			return list(candidate_configurations)

		delay = pulses[i]
		for interval_in in candidate_configurations:
			interval = np.array(interval_in)
			band_classification = get_band(delay, interval_in)

			for band in band_classification:
				interval = np.array(interval_in)
				interval[band][0] = min(interval[band][0], delay)
				interval[band][1] = max(interval[band][1], delay)
				candidate_configurations_next.add(tuple(map(tuple, interval)))

		candidate_configurations, candidate_configurations_next = \
			candidate_configurations_next, set()

	return list(candidate_configurations)

# END of experiment.

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

# Note that we *cannot* use CRC validity as an optimization criterion
# for unlimited brute forcing. The reason is simple: that would
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

def get_IDAM_field(decoded, field):
	for chunk in decoded:
		if chunk["header"] == "IDAM" and chunk["CRC_OK"]:
			return chunk[field]

	raise KeyError("Couldn't find a valid IDAM with field " + field)

def get_datalen(decoded): return get_IDAM_field(decoded, "datalen")
def get_track(decoded): return get_IDAM_field(decoded, "track")
def get_head(decoded): return get_IDAM_field(decoded, "head")

def demonstrate(au_name):
	pulses = get_pulse_deltas(au_name)
	decoded = decode_all(pulses)
	augmented_decoded = augment_full_brute(pulses, decoded,
		get_datalen(decoded))
	print("Reconstructing corrupted IDAMs. This might take a while...")
	reconstructed = mfm.reconstruct_idams(pulses, augmented_decoded, 
		get_track(decoded), get_head(decoded), get_datalen(decoded))

	print("\nWithout brute-forcing:")
	print_stats(decoded)
	if len(decoded) > 0:
		print("\n\t\t%d%% of the floppy data was decoded." %
			(100*get_coverage(pulses, decoded)))

	print("\nWith brute-forcing (no dewarping):")
	print_stats(augmented_decoded)
	if len(augmented_decoded) > 0:
		print("\n\t\t%d%% of the floppy data was decoded." %
			(100*get_coverage(pulses, augmented_decoded)))

	print("\nWith brute-forcing and corrupted sector detection:")
	print_stats(reconstructed)
	if len(reconstructed) > 0:
		print("\n\t\t%d%% of the floppy data was decoded." %
			(100*get_coverage(pulses, reconstructed)))

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
