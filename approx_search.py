import numpy as np
import fuzzysearch

def get_substitution_corruptions(in_list, substitutions, pos, maxpos, out_dict):
	if pos == min(maxpos, len(in_list)) or substitutions == 0:
		out_dict[tuple(in_list)] = True
		return

	original_value = in_list[pos]
		
	for i in range(-1, 3):
		if i == -1:
			get_substitution_corruptions(in_list, substitutions,
				pos+1, maxpos, out_dict)
		else:
			# Don't recurse more than needed.
			if i == original_value:
				continue

			in_list[pos] = i
			get_substitution_corruptions(in_list, substitutions-1,
				pos+1, maxpos, out_dict)

	in_list[pos] = original_value

def get_insert_corruptions(in_list, num_inserts, pos, maxpos,
	out_dict):

	if num_inserts == 0 or pos == min(maxpos, len(in_list)):
		out_dict[tuple(in_list)] = True
		return

	get_insert_corruptions(in_list, num_inserts, pos+1, maxpos, out_dict)

	for i in range(3):
            get_insert_corruptions(in_list[:pos] + [i] + in_list[pos:],
                        num_inserts-1, pos+1, maxpos, out_dict)

def get_delete_corruptions(in_list, num_deletions, pos, maxpos, out_dict):
	if num_deletions == 0 or pos == min(maxpos, len(in_list)):
		out_dict[tuple(in_list)] = True
		return

	get_delete_corruptions(in_list, num_deletions, pos+1, maxpos, out_dict)
	get_delete_corruptions(in_list[:pos] + in_list[pos+1:], num_deletions-1, pos, maxpos, out_dict)

def get_corruptions(needle, num_substs, num_inserts, num_deletes):

	out = [{}, {}, {}]

	get_substitution_corruptions(needle, num_substs, 0, len(needle),
		out[0])

	for corrupted_needle in out[0].keys():
		get_insert_corruptions(list(corrupted_needle), num_inserts, 0,
			len(corrupted_needle), out[1])

	for corrupted_needle in out[1].keys():
		get_delete_corruptions(list(corrupted_needle), num_deletes, 0,
			len(corrupted_needle), out[2])

	return out[2]

def rkarp_hash(in_list, window_length):
	# From C++11 minstd_rand LCG, as good as anything I suppose...
	# "True" Rabin-Karp lets the multiplier be a random irreducible
	# polynomial so as to thwart any attackers, but if I want precomputing
	# the hashes for corrupted MFM chunks to remain an option, I should
	# probably stick to a constant multiplier so that the user doesn't have
	# to rerun the precomputation whenever the multiplier changes.

	modulus = (1<<31) - 1
	multiplier = 48271

	hash_out = 0
	for x in in_list[:window_length]:
		hash_out = (hash_out * multiplier + x) % modulus

	return hash_out

def rkarp_check(haystack, needle_set, window_length):
	modulus = (1<<31) - 1
	multiplier = 48271

	high_multiplier = 1
	current_hash = 0

	hits = []

	for i in range(window_length-1):
		high_multiplier = (high_multiplier * multiplier) % modulus

	for i in range(len(haystack)):
		# Remove the leftmost byte
		if i >= window_length:
			current_hash = (current_hash - haystack[i-window_length] * high_multiplier) % modulus

		# Add the rightmost byte
		current_hash = (current_hash * multiplier + haystack[i]) % modulus
		start = max(i-window_length+1, 0)
		if current_hash in needle_set:
			hits.append(start)
		#print("Is %d, should be %d" % (current_hash, rkarp_hash(haystack[start:i+1], window_length)))

	return hits

def approximate_rkarp_match(haystack, needle, window_length, max_substs,
	max_inserts, max_dels, corruptions=None):

	# This can be precomputed and stored in a file for speedup purposes.
	if corruptions is None:
		corruptions = get_corruptions(needle[:window_length],
			max_substs, max_inserts, max_dels)

	min_window_length = window_length - max_dels

	chopped_corruptions = list(set([x[:min_window_length] 
		for x in corruptions.keys()]))
	rkarp_needles = set([rkarp_hash(x, min_window_length) for x in chopped_corruptions])

	return rkarp_check(haystack, rkarp_needles, min_window_length)

# This function is used for searching against a common prefix that all A1A1A1
# chunks have in common: three times the A1 (with clock error for OOB), plus
# four high bits (the F nibble that every indicator byte starts with). See
# mfm.py.

def approximate_match(haystack, needle, window_length, max_substs, max_inserts,
	max_dels, corruptions=None):

	definite_matches_by_end = {}

	possible_matches = approximate_rkarp_match(haystack, needle, 
		window_length, max_substs, max_inserts, max_dels, corruptions)

	# Required because fuzzysearch works on strings.
	needle_str = needle.astype(np.uint8).tostring()

	for rkarp_match in possible_matches:
		start_pos = rkarp_match
		end_pos = rkarp_match + len(needle) + max_inserts

		# Fuzzysearch works on strings, and we're only interested in
		# searching the window given by Rabin-Karp.
		hay_str = np.array(haystack[start_pos:end_pos+1]).astype(
			np.uint8).tostring()

		results_this = fuzzysearch.find_near_matches(needle_str,
			hay_str, max_substs, max_inserts, max_dels)

		for result in results_this:
			# Reconstruct absolute positions.
			result_absolute = (result.start + start_pos,
				result.end + end_pos, result.dist,
				result.matched)

			# Check if this result has already been registered.
			# This can happen e.g. if "ABC" matches "X[YABD]" with
			# distance 2, and then later "XY[ABD]" with distance 1.
			if (result.end + end_pos) in definite_matches_by_end:
				incumbent = definite_matches_by_end[result.end + end_pos]
				# Compare distances, best distance wins.
				if result_absolute[2] < incumbent[2]:
					definite_matches_by_end[result.end + end_pos] = result_absolute
			else:
				definite_matches_by_end[result.end + end_pos] = result_absolute
			
			#definite_matches.append((result.start + start_pos,
			#result.end + end_pos, result.dist, result.matched))

	# Get the unique results in sorted order.
	return sorted(list(definite_matches_by_end.values()))

def find_specific_chunk(chunk_needle, haystack, common_prefix, 
	common_prefix_matches, max_substs, max_inserts, max_dels):

	max_levenshtein = max_substs + max_inserts + max_dels
	results = []

	chunk_needle_str = chunk_needle.astype(np.uint8).tostring()
	haystack_str = haystack.astype(np.uint8).tostring()

	for i in range(len(common_prefix_matches)):
		# We have an approximate match against the common prefix.
		# To speed things up, now search for an approximate match
		# for the remaining part of the string only, with a max
		# Levenshtein distance that makes the distance for the whole
		# string (common prefix plus unique) within our bounds.

		common_length = len(common_prefix)
		corrupted_common_match_len = len(common_prefix_matches[i][3])
		start_pos = common_prefix_matches[i][0]
		
		needle_not_common = chunk_needle_str[common_length:]
		hay_not_common = haystack_str[start_pos+corrupted_common_match_len:
			start_pos+corrupted_common_match_len+len(chunk_needle)-common_length+max_inserts]
		fuzzy_result = fuzzysearch.find_near_matches(needle_not_common,
			hay_not_common, max_substs, max_inserts, max_dels,
			max_levenshtein - common_prefix_matches[i][2])

		# The actual fuzzy match data gives no useful information
		# because all the offsets are off. So instead just record the
		# match start position given by common_prefix_matches, and end
		# position according to the fuzzy search result.

		# Sometimes we get a false positive (match) with a nonzero
		# start index, which would count as an insertion error over the
		# whole string. So do a search of the whole string if we get a
		# match, and only *then* admit.

		if len(fuzzy_result) > 0:
			fuzzy_result = fuzzysearch.find_near_matches(
				chunk_needle_str, haystack_str[start_pos:
					start_pos+len(chunk_needle)+max_inserts],
				max_substs, max_inserts, max_dels)

		if len(fuzzy_result) > 0:
			result = fuzzy_result[0]
			results.append((result.start + start_pos, result.end + start_pos, result.dist))

	return results
