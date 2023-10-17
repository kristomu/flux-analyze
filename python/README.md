# flux-analyze
Analysis and recovery of IBM MFM floppy image data

This is where I experiment with methods to more reliably decode IBM MFM floppies
based on flux transition data.

Currently, mfm_decoder contains some functions for doing so, based on the .au
format that FluxEngine exposes through fluxtoau. The main idea is to use an
optimal O(kn log n) 1D k-median algorithm, referenced in Grønlund et al. 2017,
to separate different clock tracks in the presence of noise. More info can be
found in the source files' comments.

To use, run mfm_decoder interactively in ipython, i.e. `ipython3 -i mfm_decoder.py`, then invoke `demonstrate("your-file.au")`

## Experimental local dewarping

mfm_decoder also contains functions to do a sliding window k-medians-based
correction of jitter/warping. This makes sectors easier to read in some cases,
and could potentially recover data from fragile disks where every read counts.
However, it is slow at the moment.

`demonstrate()` demonstrates the dewarping functionality, and passing
`show_plot=True` will show some before-and after plots, like this:

Before:
<img src="doc/MSPlus_warped_before.svg" width=300>
After:
<img src="doc/MSPlus_warped_after.svg" width=300>
