Test cases for C++ flux-analyze:

- x33_Unknown_AM: This is one track from a DOS 5.0 disk. It contains some flux delays where random noise in the flux delays makes the ordinal pattern (increase decrease e.g.) just happen to be the same as a preamble (A1A1A1/C2C2C2) pattern. So C++ flux-analyze matches these, leading to false positives that cause the preceding DAMs to be truncated and thus ignored.

- 111_clock_failure_1: Produced a "found preamble but then couldn't" error with flux_analyze. Taken from a Microsoft Plus disk. This is caused by a false positive (noise interpreted as signal) that happens to look like an A1A1A1 preamble, but that then overlaps a *real* A1A1A1 preamble right after. Thus the mfm train becomes too short and the decoder got confused. Partial fix (very hacky) in place.

- 111_clock_failure_2: Produces the same error as clock_failure_1, but for a different reason. The bug seems to be related to the C2 preamble not matching something that it ought to match, but I'm still testing.

- RETRY-77-4-t69.0: This is the same file as the one in the parent directory. It contains a chunk that can't be assigned a clock due to the bands being too close to each other. My Python code deals with this by not using a clock at all, just explicitly dealing with the bands.

- x16_copying_uninitialized: This is a track from another DOS 5.0 disk. It caused an invalid read of two bytes (detected by Valgrind) due to a preamble signature being found right at the start of a track. Fixed.

- 90-DDAM: A (miscategorized) DDAM from a floppy with only 0xf6. It triggered a bug where DDAMs were set on admark.dam instead; fixed.
