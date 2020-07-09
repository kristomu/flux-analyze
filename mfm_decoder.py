import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict, deque
import struct, timeit, copy
import kmedian
import kcenter
import re

from scipy.stats import trim_mean

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
		#global_clusters = kcenter.wt_optimal_k_center(pulse_deltas, bands, 3)[0]

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
#	0	is represented by	RN	if we last observed a 0 bit.
#	0	is represented by	NN	if we last observed a 1 bit.
# (https://info-coach.fr/atari/hardware/FD-Hard.php#Modified_Frequency_Modulation)

# A 0 coding that doesn't fit (e.g. NN with previous bit being 1)
# is an error. This either signifies out-of-band data (e.g. in the
# A1A1A1 preamble) or a literal error (demagnetized reversal, delta
# close to the boundary between one and two, or similar errors).

# TODO? Clocking backwards, when the flux record doesn't start at a
# sector boundary. The best approach is probably to try every offset
# to see which works.

# The illegal tuples are NR RN and RN NN. Complicating the fact is that
# delays do not map strictly to tuples, e.g. a delay in the third band
# may be perfectly allowed in one case and produce an error in another.
# (E.g. NR NN NR is allowed, but RN NN RN is not.)

# If we simply want to set the number of errors to zero, it should be
# no problem to happen upon the location of the first error, and then
# figure out what ranges of clustering will handle that error. Then do
# the same for each subsequent error. However, that may produce wrong
# results if there's a missing flux transition in there, after which
# *everything* after that point becomes corrupted.

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

zero_train_len = 7 * 8

short_A1_code = [0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1] * 3
A1_prefix_code = [1, 0] * zero_train_len + short_A1_code

# Ditto here, only the pattern is
# 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0
# [ 1 ] [ 1 ] [ 0 ] [ 0 ] [o:0] [ 0 ] [ 1 ] [ 0 ]

short_C2_code = [0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0] * 3
C2_prefix_code = [1, 0] * zero_train_len + short_C2_code

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
# IDAM			4				4			A1 A1 A1 Fx
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
def decode_all(pulse_deltas, datalen=0):
	pulse_runs = get_runs(pulse_deltas, min_run_length=80)

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
				# Don't trust the data length if CRC isn't OK.
				if "CRC_OK" in struct and struct["CRC_OK"]:
					datalen = struct["datalen"]

			decoded_structures.append(struct)
		except KeyError:
			print("Could not distinguish the floppy chunk header!")
			pass

	return decoded_structures

# For debugging. HACK. TODO: Fix later.

def decode_all_from_MFM(assignments, MFM_code, MFM_pos, struct_starts=None,
	datalen=0):
	if struct_starts is None:
		struct_starts = get_struct_starts(assignments, MFM_code, MFM_pos)

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
			# Don't trust the data length if CRC isn't OK.
			if "datalen" in struct:
				if "CRC_OK" in struct and struct["CRC_OK"]:
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

# warped.au
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

# "RETRY-77-4-t69.0.au" is a nice new standard.
# First you have to find the right clustering. Best so far is 22, 33, 44.
# Then dewarping helps immensely. But fluxengine gets 12 good sectors
# *without* dewarping. I'm not there yet!

# Also https://stats.stackexchange.com/questions/6920/efficient-online-linear-regression

def assign(pulses, clustering):
	assignments = []

	for i in range(len(pulses)):
		assignment = np.argmin(np.abs(pulses[i] - clustering))
		assignments.append(assignment)

	return assignments

def simplest_ever_pll(pulses, x=None, alpha=0.10):
	if x is None:
		x = kmedian.wt_optimal_k_median(pulses[:3000], 3)
	assignments = []

	for i in range(len(pulses)):
		assignment = np.argmin(np.abs(pulses[i] - x))
		assignments.append(assignment)
		residual = pulses[i] - x[assignment]

		# Relative residual
		rel_residual = residual / x[assignment]

		# Adjust the pulse delays accordingly.
		x = x + x * rel_residual * alpha

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

# detecting A1A1A1 preamble starts without knowing the clustering.
# We have a bunch of zeroes, which are encoded as RNRNRN..
# then we have the pattern NRNNNRNNRNNNRNNR thrice.
# This is equivalent to, in terms of delays, to

#  1   2   3    2   3    2
# RN RNN RNNN RNN RNNN RNN R...

# or, with two runs,

#  1  2    3    2   3    2   1   3   2    3    2
# RN RNN RNNN RNN RNNN RNN R N RNNN RNN RNNN RNN R...

# That, in turn, means that the second delay must be longer than the
# first, the third longer than the second, the fourth shorter than the
# third, and so on... and that can be searched for without knowing the
# exact delays!

# Delta wise we have, for the two runs above:
# 1->2		increase
# 2->3		increase
# 3->2		decrease
# etc..

# This will not work if there are runs of the same pulse length, because
# we can't track equality in the presence of noise. (More robust three-valued
# logic can do so, by essentially treating equality as wildcards. Not
# implemented here.)
preamble_pulse_sequence = [1, 2, 3, 2, 3, 2] + [1, 3, 2, 3, 2] * 2

# This one can not be used in get_potential_preambles because of its
# constant train of ones. However, it can be used when trying to find the
# proper band parameters, later. It's shorter than the full A1 sequence
# so that we can discover entries that have garbled zero trains in front.
extended_preamble_pulse_sequence = [1] * 5 + preamble_pulse_sequence

full_preamble_pulse_sequence = [1] * zero_train_len + preamble_pulse_sequence

def get_potential_preambles(pulses):
	pulse_differential = pulses[1:] > pulses[:-1]
	differential_preamble = np.array(preamble_pulse_sequence[1:]) > \
		np.array(preamble_pulse_sequence[:-1])

	# https://stackoverflow.com/questions/4664850
	return [m.start() for m in re.finditer(differential_preamble.tostring(),
		pulse_differential.tostring())]

# Once we've found a potential preamble, dump the numerical values into
# separate categories corresponding to the bands they must be in (if the
# match is real). This is used to determine the range of parameters where
# the assignment will classify the pulse delays into their correct bands,
# and thus will detect the preamble.
def categorize_delays(pulses, start_pos, preamble_sequence):
	points_by_band = [[] for x in range(3)]	# everything is a reference...

	# These are offset from zero, so points_by_band[0] is the first band.
	# For MFM, this first band has the value of somewhere around 2k+b,
	# for unknown values k and b (that we have to find later).
	for i in range(len(preamble_sequence)):
		points_by_band[preamble_sequence[i]-1].append(pulses[start_pos+i])

	return [np.array(x) for x in points_by_band]

# Given a position where we recognized a preamble sequence, categorize
# the points in the vicinity according to the full preamble pulse
# sequence.

def additional_pulses_extended_preamble():
	return len(extended_preamble_pulse_sequence) - len(preamble_pulse_sequence)

def additional_pulses_full_preamble():
	return len(full_preamble_pulse_sequence) - len(preamble_pulse_sequence)

def categorize_extended_preamble(pulses, partial_preamble_start):
	longer_length_full = len(extended_preamble_pulse_sequence) - len(preamble_pulse_sequence)

	return categorize_delays(pulses,
		partial_preamble_start - longer_length_full,
		extended_preamble_pulse_sequence)

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

# Get appropriate cluster intervals. It returns three tuples or None if
# it's impossible to assign clusters to tell the bands apart.

# If it doesn't return None, then row k of the matrix gives the inclusive
# interval that must be assigned to the kth band.
def get_cluster_intervals(points_by_band):
	band_intervals = np.array([[min(x), max(x)] for x in points_by_band])
	for i in range(len(band_intervals)-1):
		# If our max is greater than the next one's min, there's an
		# overlap and no classification can possibly work.
		if band_intervals[i][1] >= band_intervals[i+1][0]:
			return None
	return band_intervals

def get_pruned_potential_preambles(pulses):
	gpp = get_potential_preambles(pulses)
	accepted_preambles = []

	for preamble_pos in gpp:
		points_by_band = categorize_extended_preamble(pulses, preamble_pos)
		if get_cluster_intervals(points_by_band) is not None:
			accepted_preambles.append(preamble_pos)

	return accepted_preambles

# Surprisingly efficient function: categorize every preamble, then
# take the median of each category as the cluster center.

def cluster_on_preambles(pulses, gpp=None, trim=None):
	if gpp == None:
		gpp = get_pruned_potential_preambles(pulses)

	cats = [np.array([]) for x in range(3)]
	for i in gpp:
		cat = (categorize_extended_preamble(pulses, i))
		for j in range(len(cats)):
			cats[j] = np.concatenate([cats[j], cat[j]])

	if trim is not None:
		return ([trim_mean(x, trim) for x in cats])

	return ([np.median(x) for x in cats])

# Create a random cluster that assigns the intervals to the bands.
def random_cluster(intervals):

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
		x = np.random.uniform(0, intervals[0][1])
		y = np.random.uniform(2 * intervals[0][1] - x + eps, 2 * intervals[1][0] - x)
		z = np.random.uniform(2 * intervals[1][1] - y + eps, 2 * intervals[2][0] - y)

	return np.array([x, y, z])

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

def debug_assn_problems(pulses, preambles_clusterings, default_clustering,
	min_pos, max_pos):

	start_pos = preambles_clusterings[min_pos][0]
	end_pos = preambles_clusterings[max_pos+1][0]

	gpp_short = preambles_clusterings[min_pos:max_pos]
	pulses_short = pulses[start_pos:end_pos]

	direct_assignment = assign(pulses_short, default_clustering)
	indirect_assignment = assign_on_preambles_clusterings(pulses_short,
		gpp_short, default_clustering)

	# Indirect assignment should always be better, but isn't for some
	# reason...
	direct_stats = get_invalid_chunk_stats(decode_from_assignments(direct_assignment))
	indirect_stats = get_invalid_chunk_stats(decode_from_assignments(indirect_assignment))

	print_stats(decode_from_assignments(direct_assignment))
	print_stats(decode_from_assignments(indirect_assignment))

	print("Detected chunks, direct:", direct_stats[0]+direct_stats[1])
	print("Detected chunks, indirect:", indirect_stats[0]+indirect_stats[1])

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

def demonstrate(au_name, show_plot=False, global_clusters=None):
	pulses = get_pulse_deltas(au_name)
	assign_wo_dewarp, pulses_copy = get_assignments(pulses, correct_warping=False,
		global_clusters=global_clusters)
	#assign_with_dewarp, dewarped_pulses = get_assignments(pulses, correct_warping=True,
	#	global_clusters=global_clusters)

	assign_with_dewarp = shifted_window_pll(pulses, 100, len(pulses)-100, global_clusters)

	#if show_plot:
	#	scatter_plot(pulses, 30000, step=5)
	#	scatter_plot(dewarped_pulses, 30000, step=7)

	# Without dewarping
	without_decoded = decode_from_assignments(assign_wo_dewarp)
	without_cat, without_cat_valid = categorize_decoded(without_decoded)

	# With dewarping
	dewarp_decoded = decode_from_assignments(assign_with_dewarp)
	dewarp_cat, dewarp_cat_valid = categorize_decoded(dewarp_decoded)

	# Together
	concert_decoded = without_decoded + dewarp_decoded

	print("Stats: ")
	print("\tWithout dewarping:")
	print_stats(without_decoded)
	print("")
	print("\tWith dewarping:")
	print_stats(dewarp_decoded)
	print("")
	print("\tAuxiliary info (both combined):")
	print_stats(concert_decoded)

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
