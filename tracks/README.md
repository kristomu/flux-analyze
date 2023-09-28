Odd track examples:

- low_level_format_with_noise: every sector is 0xF6 (i.e. only a low level format, no file system on top). Some weird effects on the flux transitions in certain places. Useful for testing recovery ideas later, because the data pattern on the corrupted sectors is known. This track also contains IAMs (start-of-track markers). 
- MS_Plus_OK_track - OK track from a Microsoft Plus! floppy. This one can be decoded fine with non-k-median FluxEngine.
- MS_Plus_warped_track - Warped track from a Microsoft Plus! floppy. FluxEngine can't cleanly decode every sector, but mfm_decoder can.

- amiga: David Given's Amiga Workbench example track. This contains a very heavy
zero band which pulls ordinary k-median clustering towards itself, and so it falsely places two clusters on the lower band.

- RETRY-77-4-t69.0: A long track with a very heavy zero band and significant smearing near the beginning. Exercises a bug where FluxEngine's CDF sumToHere counter overflows. (To trigger, try to do a 10-clustering or thereabouts; it returns 0 as a cluster, but shouldn't. Or compare its 3-clustering to mfm_decoder's.)

- RETRY-77-4-t73.0: A quite damaged track that non-kmedian FluxEngine gets a few sectors out of, but kmedian simply gives up. It lacks the third band altogether.

- multitrack/low_level_format_with_noise.flux: The whole disk of low_level_format_with_noise. This can be used to reproduce a false positive in C++ flux-analyze in combination with brute-force dewarping: a sector is mistakenly identified as valid because the corruption of data and CRC just happen to match. This does not happen with the alternate "historical" dewarping strategy.
