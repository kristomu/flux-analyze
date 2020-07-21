import numpy as np
import low_level
import struct
import re

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

# -- error correction ideas --

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

# -- /error correction --

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

# Numpy version of subfinder, using the cast-to-string trick from ordinal.
def subfinder_np(values, pattern):
	return np.array([m.start() for m in re.finditer(
		np.array(pattern).astype(np.uint8).tostring(),
		np.array(values).astype(np.uint8).tostring())])

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

preambles_list = [A1_prefix_code, C2_prefix_code]
short_preambles_list = [short_A1_code, short_C2_code]

# ----------------------------------------------------------

# Finding and interpreting MFM stuff, in phases:
# (Zeroth phase: loading low-level data)
# First phase: detecting preambles,
# Second phase: determining the correct clustering,
# Third phase: decoding the relevant data,
# Fourth phase: categorizing the results.

# First and second phase stuff below, as pertains to MFM

# First:

def get_preamble_positions(pulses):
	return np.sort(np.concatenate(low_level.search_for_patterns(pulses,
		[low_level.encode_to_assignment(np.array(A1_prefix_code)),
		 low_level.encode_to_assignment(np.array(C2_prefix_code))])))

# Second, with a fixed cluster for the whole track:

def cluster_on_preambles(pulses, trim=None):
	return low_level.cluster_on_custom_preambles(pulses,
		[low_level.encode_to_assignment(np.array(A1_prefix_code)),
		 low_level.encode_to_assignment(np.array(C2_prefix_code))],
		trim)

# (Where does the brute forcing clustering stuff go? probably
# elsewhere.)

# But it's better to specify the band intervals separately. (Perhaps
# somehow choose the clustering consistent with the band intervals that
# is simultaneously the closest to the global clustering, but I don't
# know if that's required yet...)

def get_preamble_positions_clusters(pulses):
	preamble_positions_clusters = []

	for preamble, short_preamble in zip(preambles_list,
		short_preambles_list):

		encoded_preamble = np.array(low_level.encode_to_assignment(
			np.array(preamble)))

		positions = low_level.search_for_pattern(pulses, encoded_preamble)

		for position in positions:
			cluster_guess = low_level.get_cluster(pulses, position,
				encoded_preamble)

			preamble_positions_clusters.append((position, cluster_guess,
				np.array(preamble), np.array(short_preamble)))

	return sorted(preamble_positions_clusters)

# Third:

# Chunk_pos_preamble is one of the list entries returned by the previous
# function. pulse_len is an upper bound of how long (in pulse delays) this
# particular chunk is. (More sophisticated ways of doing this, later).
def decode_chunk(pulses, chunk_pos_preamble, pulse_len, last_datalen,
	override_clustering = None):

	preamble_pos, clustering, preamble, short_preamble = chunk_pos_preamble

	if override_clustering is not None:
		clustering = override_clustering

	assignments = low_level.assign_partial(
		pulses, preamble_pos, preamble_pos+pulse_len, clustering)

	# Positions map flux positions to pulse positions
	fluxes, positions = low_level.get_flux_stream_and_pos(assignments)

	# The short preamble should be the bit that starts the header for the
	# chunk itself, absent any padding, e.g. A1A1A1 for IBM and A1A1 for
	# Amiga.
	code_start = subfinder_np(fluxes, short_preamble)[0]

	# And flux_pos maps bit positions to flux positions.
	outbits, errors, flux_pos = MFM_decode_from_flux_stream(
		fluxes[code_start:])

	# When counting errors, remove the intentional errors present in the
	# preamble we're using, so that these errors don't count.
	preamble_outbits, preamble_errors, preamble_fluxpos = \
		MFM_decode_from_flux_stream(short_preamble)

	# Indicates an expected error.
	errors[np.where(preamble_errors > 0)] = -2

	# Convert to bytes and count the number of MFM errors per byte.
	outbytes = np.packbits(outbits).astype(np.uint8).tostring()
	errors_per_byte = np.resize(errors, (len(errors)//8, 8))
	num_unexpected_errors_per_byte = np.sum(errors_per_byte>0, 1)

	# Convert bytes to struct.
	chunk = decode_floppy_struct(outbytes, last_datalen)
	chunk_len = get_chunk_bytelength(outbytes, last_datalen)
	errors_in_chunk = sum(num_unexpected_errors_per_byte[:chunk_len])

	# Add some auxiliary information.
	chunk["_start_pos"] = preamble_pos
	try:
		chunk["_end_pos"] = preamble_pos + positions[code_start + \
			flux_pos[chunk_len*8]] + 1
	except IndexError:
		chunk["_end_pos"] = preamble_pos + pulse_len # truncated

	chunk["_num_errors"] = errors_in_chunk

	return chunk

# Another idea: It might be that every IDAM's CRC is unique, at least over
# the span of tracks and sectors for an 1.44M floppy. So then all we need
# to do is to capture the CRC and the rest of the IDAM can be reconstructed.

# (Enumerate all CRCs for this track.
#  Search for the ordinal patterns in the vicinity of the IDAM, *or*
#  Search for the IDAM itself, with a fuzzy Levenshtein-based search.)
# (Similar tricks are possible for looking for repeat sectors, which would
#  be a quick and dirty way to get 100% completion on
#  low_level_format_with_noise. But that victory would be very hollow.)

# To do error correction:
# errors[np.where(errors>0)]
# # Observe the 124 and 148 past the three intentional IDAM byte markers
# positions[code_start+124] # = 102
# pulses_77_long[wx[0]+102] # shows a 45
# so classify it into the previous cluster
# pulses_77_long[wx[0]+102] = 35
# redo all the above
# all errors gone, CRC OK.
# Notice that there's no way to not classify 45 to the third band because
# the third band includes 44, which is less than 45.

# There used to be error correction here, but I've removed it as it wasn't
# very good.


# The following function decodes the flux stream into an MFM bit sequence.
# The errors array is, for each bit, either -1 if the pair that generated
# the bit is legal according to MFM clock constraints, otherwise an index
# to the the first bit of the flux stream pair that violates the
# constraint.

# TODO: Abort after max_length bits have been output. Then the parsing
# routine itself can read off only enough to establish what kind of struct
# it is, then get the length of the struct, then read off the
# appropriate number of bits. May not be necessary at the moment...
def MFM_decode_from_flux_stream(flux_stream, max_length=np.inf):
	outbits = []
	error = []
	flux_pos = []

	# This means "I don't know what the last bit was."
	last_bit = None

	# TODO: PERF: don't use append, use i//2 instead, and preallocate
	# the lists as np arrays.

	for i in range(0, len(flux_stream) - len(flux_stream)%2, 2):
		if len(outbits) >= max_length:
			continue

		next_pair = [flux_stream[i], flux_stream[i+1]]

		if next_pair == [0, 1]:
			outbits.append(1)
			error.append(-i)

		if next_pair == [1, 0]:
			outbits.append(0)
			if last_bit is not None and last_bit != 0:
				error.append(i)
			else:
				error.append(-i)

		if next_pair == [0, 0]:
			outbits.append(0)
			if last_bit is not None and last_bit != 1:
				error.append(i)
			else:
				error.append(-i)

		flux_pos.append(i)

		if next_pair == [1, 1]:
			raise Exception("Hard violation of clock constraint. Aborting.")

		last_bit = outbits[-1]

	return np.array(outbits), np.array(error), np.array(flux_pos)

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

#####

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

# Return the number of bytes required to parse the structure without
# truncating. This is used to count the number of errors that matter,
# i.e. those that impinge on the chunk data itself.
def get_chunk_bytelength(dec_bytes, last_idam_datalen):
	signature = dec_bytes[:4]

	if signature == b'\xa1\xa1\xa1\xfe':
		return 10

	if (signature == b'\xa1\xa1\xa1\xfb') or (signature == b'\xa1\xa1\xa1\xf8'):
		return 4 + last_idam_datalen + 3

	if signature == b'\xc2\xc2\xc2\xfc':
		return 4

	print([hex(x) for x in signature])

	raise KeyError("get_chunk_bytelength: could not identify signature")

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
		if last_idam_datalen == 0 or last_idam_datalen is None:
			for i in reversed(data_length):
				speculative_data_chunk = decode_floppy_struct(dec_bytes, i)
				if "CRC_OK" in speculative_data_chunk and \
					speculative_data_chunk["CRC_OK"]:
					return speculative_data_chunk

		# Every attempt gave a bad result. Abort.
		if last_idam_datalen == 0 or last_idam_datalen is None:
			raise IndexError("Could not determine DAM data length!")

		floppy_info["data"] = dec_bytes[4:4+last_idam_datalen]
		if len(dec_bytes) > 4+last_idam_datalen+2: floppy_info["CRC"] = struct.unpack('>H', dec_bytes[4+last_idam_datalen:4+last_idam_datalen+2])[0]
		if len(dec_bytes) > 4+last_idam_datalen+2: floppy_info["CRC_OK"] = crc16(dec_bytes[:4+last_idam_datalen], 0, 4+last_idam_datalen) == floppy_info["CRC"]
		if len(dec_bytes) < 4+last_idam_datalen+3: floppy_info["truncated"] = True

		return floppy_info

	if signature == b'\xc2\xc2\xc2\xfc':
		floppy_info["header"] = "IAM"

		return floppy_info

	print([hex(x) for x in signature])

	raise KeyError

# EXPERIMENTAL
# I don't know if this is the approach I want to take... I can be much
# more generous (as long as I've already excluded the sectors that aren't
# valid), because I can rely on a CRC.

# track 77, head 0, sector 17, datalen 512, CRC 40454, CRC OK:
# array([1, 2, 1, 2, 1, 0, 2, 1, 2, 1, 0, 2, 1, 2, 1, 0, 0, 0, 0, 0, 0, 0,
#       1, 1, 1, 1, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0,
#       0, 0, 0, 1, 2, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0])

# test_asn = above
# sp = test_asn.astype(np.uint8).tostring()
# brute_assignment = low_level.assign(pulses_77_long, cluster_on_preambles(
#	pulses_77_long))
# hp = brute_assignment.tostring()
# fuzzysearch.find_near_matches(sp, hp, max_substitutions=5,
#	max_insertions=0, max_deletions=0)

# full_assignments = assign(pulses_77_long, cluster_on_preambles(
#	pulses_77_long))
# hp = full_assignments.astype(np.uint8).tostring()
# # Take these numbers from an intact IDAM's start pos and end pos.
# flux_sector8 = pulses_77_long[5541480:5541597+2]
# cluster_guess = low_level.get_cluster(flux_sector8, 0,
# 	low_level.encode_to_assignment(np.array(A1_prefix_code)))
# assignments = low_level.assign(flux_sector8, cluster_guess)
# sp = assignments.astype(np.uint8).tostring()
# found = fuzzysearch.find_near_matches(sp, hp, max_substitutions=5,
#	max_insertions=0, max_deletions=0)
# zero_positions = sorted([(x.start, x.end) for x in found])
# for zp in zero_positions:
#	pulses_77_cleanup[zp[0]:zp[1]] = flux_sector8

# The same thing can really be done on any CRC-OK chunk (DAM or IDAM) to
# turn corrupted ones into clear ones, with the hope that one of them may
# be adjacent to an IDAM or DAM respectively, and thus turn a sector with
# unknown data (or a data chunk with unknown designation) into one that's
# both.