import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import struct, timeit
import kmedian

# USAGE: Create the AU file with
# fluxengine convert fluxtoau -s SOURCE/:t=TRK:s=SIDE -o somename.au
# call get_pulse_deltas on somename.au
# and call various functions (see demonstrate).

# As an example, call demonstrate("tracks/MS_Plus_warped_track.au")

# TODO? Detect an abrupt end to a sector.
# Also be able to clock backwards.

# Call the groupings (categories) "bands" to follow KF terminology.

# TODO? Do something clever with the error report when decoding MFM.

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

def scatter_plot(pulse_deltas, len, start=0, step=1):
	plt.scatter(range(0, len-start, step), pulse_deltas[start:len:step], s=10)
	plt.show()

# Returns the assignment to bands (shortest distance between reversals
# being 0). As second return value, it returns either the input pulse array
# (if correct_warping is false), or the reconstructed pulse array, corrected
# for warping, if correct_warping is true.

def get_assignments(pulse_deltas, bands=3, correct_warping=True,
	global_clusters=None):

	if global_clusters is None:
		global_clusters = kmedian.wt_optimal_k_median(pulse_deltas, bands)

	if correct_warping:
		print("Correcting warping. This may take some time...")
		pulse_clusters = kmedian.local_k_median(pulse_deltas,
			global_clusters=global_clusters, num_clusters=bands)
		assignments, residuals = kmedian.assign_to_local_clusters(
			pulse_deltas, pulse_clusters)
		pulses_reconstruction = kmedian.correct_warp(assignments,
			residuals, global_clusters)
		return assignments, pulses_reconstruction
	else:
		assignments, residuals = kmedian.assign_to_clusters(pulse_deltas,
			global_clusters)
		return assignments, pulse_deltas

# Interpret MFM.

# The encoding table is, where R is a reversal and N is none:
#	1	is represented by	NR
#	0	is represented by	RN	if we last observed a 0 bit
#	0	is represented by	NN	if we last observed a 1 bit.
# (https://info-coach.fr/atari/hardware/FD-Hard.php#Modified_Frequency_Modulation)

# A 0 coding that doesn't fit (e.g. NN with previous bit being 1)
# is an error. This either signifies out-of-band data (e.g. in the
# A1A1A1 preamble) or a literal error (demagnetized reversal, delta
# close to the boundary between one and two, or similar errors).

# TODO? Clocking backwards, when the flux record doesn't start at a
# sector boundary. The best approach is probably to try every offset
# to see which works.

def assignments_to_MFM_code(assignments):
	# A flux reversal is symbolized by a zero bit, and
	# no flux reversal by a one. The position gives the pulse number
	# that corresponds to this particular bit; it's used in later
	# lookups.

	# The function returns the MFM code as the first list, and the
	# positions as the second.

	# TODO? Probabilistic judgment of where a possible extra reversal
	# might occur?

	MFM_code = []
	position = []
	counter = 0

	for assignment in assignments:
		# Category X corresponds to 1+X no-flux-reversals, then a reversal.
		MFM_code += [0] * (assignment + 1) + [1]
		position += [counter] * (assignment + 2)
		counter += 1

	return MFM_code, position

# https://stackoverflow.com/a/12576755 modified to give index.
# Detect any instances of 'pattern' in 'mylist'. Used to find the
# A1A1A1 out-of-band preamble. TODO later if required: fuzzy search
# (Levenshtein).
def subfinder(mylist, pattern):
	matches = []
	for i in range(len(mylist)):
		if mylist[i] == pattern[0] and mylist[i:i+len(pattern)] == pattern:
			matches.append(i)
	return matches

# https://www-user.tu-chemnitz.de/~heha/basteln/PC/usbfloppy/floppy.chm/
# The minimum number of zero values is 8 bytes. Since we can't detect
# whether the first byte is an RN or an NN, we can only look for 7 bytes.

# [1,0] pattern corresponds to a 0 bit with a preceding 0. There are
# a minimum of 8 bytes of 0 before the preamble starts, and we don't know
# if the first byte has a 1-bit just ahead of it. So 7 bytes * 8 bits.

# Then we have:
# 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1,
# [ 1 ] [ 0 ] [ 1 ] [ 0 ] [ 0 ] [o:0] [ 0 ] [ 1 ]
# A1 (where "o:" indicates a clock error/out-of-band indicator)
# three times.

# (Strictly speaking, I could then add in the fact that the next nibble
# must be all ones, but I can't be bothered.)

short_A1_code = [0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1] * 3
A1_prefix_code = [1, 0] * (7 * 8) + short_A1_code

# Ditto here, only the pattern is
# 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0
# [ 1 ] [ 1 ] [ 0 ] [ 0 ] [o:0] [ 0 ] [ 1 ] [ 0 ]

short_C2_code = [0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0] * 3
C2_prefix_code = [1, 0] * (7 * 8) + short_C2_code

def get_struct_starts(assignment, MFM_code, MFM_code_pos):
	struct_starts = []

	# The format of each record is (length of matching bit string,
	# beginning MFM code index, ending MFM code index).

	# This may cause trouble when I start doing error correction, but
	# I'll cross that bridge etc...

	# Get the start of every sector structure with an A1 signature.

	for A1_idx in subfinder(MFM_code, A1_prefix_code):
		beginning = A1_idx
		end = A1_idx + len(A1_prefix_code)
		struct_starts.append((end-beginning-1, beginning, end-1))

	# Second verse, C2 signature.
	for C2_idx in subfinder(MFM_code, C2_prefix_code):
		beginning = C2_idx
		end = C2_idx + len(C2_prefix_code)
		struct_starts.append((end-beginning-1, beginning, end-1))

	# Sort the structure by beginning position.
	return sorted(struct_starts, key=lambda x: x[1])

def decode_from_MFM(MFM_code):
	outbits = []
	error = []

	# This means "I don't know what the last bit was."
	last_bit = None

	for i in range(0, len(MFM_code) - len(MFM_code)%2, 2):
		next_pair = [MFM_code[i], MFM_code[i+1]]

		if next_pair == [0, 1]:
			outbits.append(1)

		if next_pair == [1, 0]:
			outbits.append(0)
			if (last_bit != 0) and (last_bit is not None):
				error.append(i)

		if next_pair == [0, 0]:
			outbits.append(0)
			if last_bit != 1 and (last_bit is not None):
				error.append(i)

		last_bit = outbits[-1]

	return outbits, error

# https://stackoverflow.com/questions/28370991
def getbytes(bits):
	done = False
	while not done:
		byte = 0
		for _ in range(0, 8):
			try:
				bit = next(bits)
			except StopIteration:
				bit = 0
				done = True
			byte = (byte << 1) | bit
		yield byte

# https://stackoverflow.com/a/55850496
def crc16(data, offset, length):
	if data is None or offset < 0 or offset > len(data)- 1 and offset+length > len(data):
		return 0
	crc = 0xFFFF
	for i in range(0, length):
		crc ^= data[offset + i] << 8
		for j in range(0,8):
			if (crc & 0x8000) > 0:
				crc =(crc << 1) ^ 0x1021
			else:
				crc = crc << 1
	return crc & 0xFFFF

def decode_bit_train(bits):
	# Because we look for a shorter run of zeroes than the specs set
	# as standard, we may not arrive exactly at the start of such a run.
	# Thus our byte boundaries might be out of sync if we start at the zero
	# offset. Fix this by skipping directly to the first 1 bit, which
	# is the start of the A1A1A1 (or CFCFCF) run.

	# TODO: Is that really necessary now?

	first = bits.index(1)
	dec_bytes = bytearray([b for b in getbytes(iter(bits[first:]))])

	return dec_bytes

# EXPERIMENTAL
# I don't know if this is the approach I want to take... I can be much
# more generous (as long as I've already excluded the sectors that aren't
# valid), because I can rely on a CRC.
def acceptable_fuzzy_match(MFM_code):

	def check_on_pattern(long_prefix_code, short_prefix_code):
		# It must match the longer sequence (with null bytes) up to a
		# Levenshtein distance of 15

		res = fuzzysearch.find_near_matches(long_prefix_code, MFM_code,
			max_l_dist=15)

		if res == []:
			return False

		# It must match the A1A1A1/C2C2C2 part up to an L-distance of 6.
		res = fuzzysearch.find_near_matches(short_A1_code, MFM_code,
			max_l_dist=6)

		if res == []:
			return False

		# And that part must have at least two errors.
		outbits, errors = decode_from_MFM(res[0].matched)

		return len(errors) >= 2

	return check_on_pattern(A1_prefix_code, short_A1_code) or \
		check_on_pattern(C2_prefix_code, short_C2_code)

# The higher order floppy formatting is here:
# https://www-user.tu-chemnitz.de/~heha/basteln/PC/usbfloppy/floppy.chm/

# Abbreviations:
# IAM	Index Address Mark			0xFC	start of a track
# IDAM	ID Address Mark				0xFE	start of a sector
# DAM	Data Address Mark			0xFB	start of sector data
#		Deleted Data Address Mark	0xF8	start of deleted sector data

# An IDAM and DAM example:

# '0xa1L', '0xa1L', '0xa1L', '0xfeL'		IDAM marker
# '0x4L',									track
# '0x1L',									head
# '0xbL',									sector
# '0x2L'									?? (sector data length?)
# '0xd8L', '0x65L'							CRC - *of the header*
# '0x4eL' ... '0x4eL'						gap (22x)
# '0xa1L', '0xa1L', '0xa1L', '0xfbL'		DAM marker
# '0x78L', '0xbeL', '0x92L', '0x39L', ...	data (512b)

# More structure info (for finding the end of the struct, for backwards
# decoding and coverage purposes)

#				----Number of bytes----
#				Standard		Minimum		Byte pattern
# (IAM starts)
# Gap 0			80				???			4E
# Gap 2/Sync	12				8			00
# IAM			4				4			C2 C2 C2 FC
# Gap 1			50				10			4E
# (IDAM starts)
# Gap 2/Sync	12				8			00
# IDAM			4				4			A1 A1 A1 FE
# Track num		1				1			varies
# Side num		1				1			varies
# Sector num	1				1			varies
# Size code		1				1			varies
# CRC			2				2			varies
# Gap 3a		22				22			4E
# (IDAM ends)
# (DAM starts)
# Gap 3b/Sync	12				12			00
# DAM			4				4			A1 A1 A1 FB, or A1 A1 A1 F8
# Data			varies			varies		varies
# CRC			2				2			varies
# Gap 4			40				1			4E
# Gap 5			664?			16			4E
# (DAM ends)
# (up to the beginning from Gap 2)

# Apparently 4E gap bytes are irrelevant, and could be anything; i.e.
# we cannot assume that they are 4E. Sync bytes are fixed, however.

# 4.13.3: should be easily detectable with the stuff below.
# Intra-sector clock rate variation is what we deal with in dewarp.

# Since there are so many values for the gaps, I'm just going to
# consider any length possible, down to a single byte. The detector
# simply has to detect a run of, say, 4E, and then keep going for as
# long as that run keeps going. Since a CRC could conceivably be 4E4E,
# this is going to be a bit of trouble for my idea of reverse-decoding
# DAM structs (records, chunks) to get their CRC...

# Looks like I'll need two hypotheses: one where the CRC is not 4E4E, and
# one where it is. There'll probably be some standard Gap 4 + Gap 5 length
# for PC formatted floppies anyway.
# (Strictly speaking, there are three hypotheses: the CRC is 4E4E, XX4E,
# or XXXX. If the gap is of standard length, then we're done, no problem.)

# TODO: Do something when last_idam_datalen is unknown...
# TODO? look for a definite end to a struct instead of just going to
# the beginning of the next one - would handle sector-within-sector.
def decode_floppy_struct(dec_bytes, last_idam_datalen):
	signature = dec_bytes[:4]
	floppy_info = {}

	# http://dmweb.free.fr/files/Atari-Copy-Protection-V1.4.pdf p. 12
	data_length = [128, 256, 512, 1024]

	# Gracefully handle partially recovered headers.
	if signature == b'\xa1\xa1\xa1\xfe':
		floppy_info["header"] = "IDAM"
		if len(dec_bytes) > 4:	floppy_info["track"] = dec_bytes[4]
		if len(dec_bytes) > 5: floppy_info["head"] = dec_bytes[5]
		if len(dec_bytes) > 6: floppy_info["sector"] = dec_bytes[6]
		if len(dec_bytes) > 7: floppy_info["datalen"] = data_length[dec_bytes[7] & 3]
		if len(dec_bytes) > 9: floppy_info["CRC"] = struct.unpack('>H', dec_bytes[8:10])[0]
		if len(dec_bytes) > 9: floppy_info["CRC_OK"] = crc16(dec_bytes[:8], 0, 8) == floppy_info["CRC"]
		if len(dec_bytes) < 10: floppy_info["truncated"] = True

		return floppy_info

	if signature == b'\xa1\xa1\xa1\xfb':
		floppy_info["header"] = "DAM"

	if signature == b'\xa1\xa1\xa1\xf8':
		floppy_info["header"] = "DAM (deleted)"

	if (signature == b'\xa1\xa1\xa1\xfb') or (signature == b'\xa1\xa1\xa1\xf8'):
		# If this is the first sector, i.e. we don't know the last IDAM
		# data length, recurse on every possible data length and return OK
		# on the first one with CRC OK.
		if last_idam_datalen == 0:
			for i in reversed(data_length):
				speculative_data_chunk = decode_floppy_struct(dec_bytes, i)
				if "CRC_OK" in speculative_data_chunk and \
					speculative_data_chunk["CRC_OK"]:
					return speculative_data_chunk

		floppy_info["data"] = dec_bytes[4:4+last_idam_datalen]
		if len(dec_bytes) > 4+last_idam_datalen+2: floppy_info["CRC"] = struct.unpack('>H', dec_bytes[4+last_idam_datalen:4+last_idam_datalen+2])[0]
		if len(dec_bytes) > 4+last_idam_datalen+2: floppy_info["CRC_OK"] = crc16(dec_bytes[:4+last_idam_datalen], 0, 4+last_idam_datalen) == floppy_info["CRC"]
		if len(dec_bytes) <= 4+last_idam_datalen+2: floppy_info["truncated"] = True

		return floppy_info

	if signature == b'\xc2\xc2\xc2\xfc':
		floppy_info["header"] = "IAM"

		return floppy_info

	print([hex(x) for x in signature])

	raise KeyError

# It's probably better to directly search for the out-of-band A1A1A1 tokens.
# But then I would have to preprocess and classify the whole pulse stream
# at once, first. Perhaps do that later.
def decode_all(pulse_deltas):
	pulse_runs = get_runs(pulse_deltas, min_run_length=80)
	datalen = 0

	decoded_structures = []

	for i in range(len(pulse_runs)):
		if i == len(pulse_runs)-1:
			dec_bytes, error = decode_pulsetrain(pulse_deltas[pulse_runs[i][1]:], 300)
		else:
			print ("Debug (%d/%d): indices %d ... %d" % (i, len(pulse_runs), pulse_runs[i][1], pulse_runs[i+1][1]))
			dec_bytes, errors = decode_pulsetrain(pulse_deltas[pulse_runs[i][1]:pulse_runs[i+1][1]], 300)

		try:
			struct = decode_floppy_struct(dec_bytes, datalen)
			# Also mark the pulse delta array index that corresponds to
			# this particular decoded structure, for diagnosis purposes.
			struct["_start_idx"] = pulse_runs[i][1]
			# The three out-of-band errors don't count.
			# The errors do not necessarily indicate a total number of
			# errors because we decode a lot more than is required.
			# All of this should be fixed later.
			struct["_num_errors"] = len(errors) - 3
			if "datalen" in struct:
				datalen = struct["datalen"]

			decoded_structures.append(struct)
		except KeyError:
			print("Could not distinguish the floppy chunk header!")
			pass

	return decoded_structures

# For debugging. HACK. TODO: Fix later.

def decode_all_from_MFM(assignments, MFM_code, MFM_pos, struct_starts=None):
	if struct_starts is None:
		struct_starts = get_struct_starts(assignments, MFM_code, MFM_pos)
	datalen = 0

	decoded_structures = []

	for i in range(len(struct_starts)):
		struct_start = struct_starts[i][1]
		if i == len(struct_starts)-1:
			struct_end = len(MFM_code)
		else:
			# We don't know how long the struct is; assume it goes
			# all the way to the next one. (May cause trouble with
			# sector-in-sector style copy protection.)
			struct_end = struct_starts[i+1][1]

		print ("Decode (%d/%d): index %d ... %d:" % (i, len(struct_starts), struct_start, struct_end), end="")
		bits, errors = decode_from_MFM(MFM_code[struct_start:struct_end])
		# We expect one clock error for each repetition of the preamble.
		# Fewer than that means we have an unexpected error (none where
		# there should be one); more also means unexpected error.
		print ("\t%d unexpected errors." % (np.abs(len(errors)-3)))

		first = bits.index(1)
		dec_bytes = bytearray(getbytes(iter(bits[first:])))

		try:
			struct = decode_floppy_struct(dec_bytes, datalen)
			# Also mark the pulse delta array index that corresponds to
			# this particular decoded structure, for diagnosis purposes.
			struct["_start_idx"] = struct_starts[i][1]
			# The errors do not necessarily indicate a total number of
			# errors because we decode a lot more than is required.
			struct["_num_errors"] = len(errors) - 3
			if len(errors) > 3:
				struct["_first_error"] = errors[3]
			if "datalen" in struct:
				datalen = struct["datalen"]

			decoded_structures.append(struct)
		except KeyError:
			print("Could not distinguish the floppy chunk header!")
			pass

	return decoded_structures

# Group DAM data chunks by (track, head, sector) of the immediately
# preceding IDAM. More tricks may be required to fail safe if there's
# corruption (e.g. the structure as intended was IDAM1 DAM1 IDAM2 DAM2,
# but DAM1 and IDAM2 got severely corrupted, then DAM2 will be attributed
# to the sector in IDAM1). A simple check to see that the headers are
# consecutive should suffice. Do that later.

def categorize_decoded(decoded_structures):
	by_sector = defaultdict(list)
	by_sector_valid = defaultdict(list)

	for i in range(1, len(decoded_structures)):
		this, prev = decoded_structures[i], decoded_structures[i-1]
		if this["header"] == "DAM" and prev["header"] == "IDAM":
			if "CRC" in this and prev["CRC_OK"]:
				ths = (prev["track"], prev["head"], prev["sector"])
				by_sector[ths].append((prev, this))
				if this["CRC_OK"]:
					by_sector_valid[ths].append((prev, this))

	return by_sector, by_sector_valid

def get_invalid_chunk_stats(decoded):
	valid = len([ x for x in decoded if "CRC_OK" in x and x["CRC_OK"]])
	invalid = len([ x for x in decoded if "CRC_OK" in x and not x["CRC_OK"]])

	return (valid, invalid)

def demonstrate(au_name, show_plot=False, global_clusters=None):
	pulses = get_pulse_deltas(au_name)
	assign_wo_dewarp, pulses_copy = get_assignments(pulses, correct_warping=False,
		global_clusters=global_clusters)
	assign_with_dewarp, dewarped_pulses = get_assignments(pulses, correct_warping=True,
		global_clusters=global_clusters)

	if show_plot:
		scatter_plot(pulses, 30000, step=5)
		scatter_plot(dewarped_pulses, 30000, step=7)

	# Without dewarping
	without_MFM, without_MFM_pos = assignments_to_MFM_code(assign_wo_dewarp)
	without_struct_starts = get_struct_starts(assign_wo_dewarp, without_MFM,
		without_MFM_pos)
	without_decoded = decode_all_from_MFM(assign_wo_dewarp, without_MFM,
		without_MFM_pos, without_struct_starts)
	without_cat, without_cat_valid = categorize_decoded(without_decoded)

	# With dewarping
	dewarp_MFM, dewarp_MFM_pos = assignments_to_MFM_code(assign_with_dewarp)
	dewarp_struct_starts = get_struct_starts(assign_with_dewarp, dewarp_MFM,
		dewarp_MFM_pos)
	dewarp_decoded = decode_all_from_MFM(assign_with_dewarp, dewarp_MFM,
		dewarp_MFM_pos, dewarp_struct_starts)
	dewarp_cat, dewarp_cat_valid = categorize_decoded(dewarp_decoded)

	print("Stats: ")
	print("\tWithout dewarping:")
	invstat = get_invalid_chunk_stats(without_decoded)
	print("\t\tValid CRC chunks: %d\tInvalid: %d\tFraction valid: %.3f" %
		(invstat[0], invstat[1], invstat[0]/(invstat[0]+invstat[1])))
	bad_sector_keys = set(without_cat.keys()) - set(without_cat_valid.keys())
	print("\t\tNumber of sectors: %d\tBad sectors: %d" % (len(without_cat.keys()), len(bad_sector_keys)))
	print("")
	print("\tWith dewarping:")
	invstat = get_invalid_chunk_stats(dewarp_decoded)
	print("\t\tValid CRC chunks: %d\tInvalid: %d\tFraction valid: %.3f" %
		(invstat[0], invstat[1], invstat[0]/(invstat[0]+invstat[1])))
	bad_sector_keys = set(dewarp_cat.keys()) - set(dewarp_cat_valid.keys())
	print("\t\tNumber of sectors: %d\tBad sectors: %d" % (len(dewarp_cat.keys()), len(bad_sector_keys)))

def time_weighted(pulses):
	print("Timing k-median calculations... This might take a while.")
	timer = timeit.Timer(lambda: kmedian.wt_optimal_k_median(pulses, 3),
		"from __main__ import np, kmedian")
	weighted_time = timer.timeit(10) / 10
	timer = timeit.Timer(lambda: kmedian.optimal_k_median(pulses, 3),
		"from __main__ import np, kmedian")
	unweighted_time = timer.timeit(1)
	print("Calculation time: unweighted: %.5f sec, weighted: %.5f sec, speed factor: %.2f x" %
		(unweighted_time, weighted_time, unweighted_time/weighted_time))

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
