# flux-analyze
Analysis and recovery of IBM MFM floppy image data

This is where I experiment with methods to more reliably decode IBM MFM floppies
based on flux transition data.

There are two source subdirectories: csrc/, the C++ version, and python/, the
Python version. I'm currently working on the C++ version; the Python version is
older.

## Request for images

I may be picking up this project again, time permitting, as I've thought of new
ideas that could handle errors more robustly. If you have got any marginal or
failing floppy images that you'd be willing to share (without PII or other
private data), let me know by issue/PR or by email to my commit Author mail
address. Such images could help me tune the algorithms so that they'll have
better success at recovering data.

## Installation and usage

To compile the C++ version, use cmake. A typical invocation would be
as follows:

    $ cd ~/Sources/
    $ git clone https://github.com/kristomu/flux-analyze.git
    $ cd flux-analyze
    $ mkdir build
    $ cd build
    $ cmake -S ../csrc -B .
    $ make

Then use by invoking `flux-analyze` by specifying one or more FluxEngine .flux
filenames on the command line. If more than one flux file is specified, they'll
be treated as different images of the same floppy.

To use the Python version, run mfm_decoder interactively in ipython, i.e.
`ipython3 -i mfm_decoder.py`, then invoke `demonstrate("your-file.au")`. The
Python version uses FluxEngine .au dumps as input. The README.md file in the
Python subdirectory contains more information.

## Tracks

The tracks/ directory contains some sample tracks in FluxEngine and .au format.
The .au files are all compressed with bz2. (Uncompress before using.)

- low_level_format_with_noise: every sector is 0x00 or 0xF6. Some weird effects on the flux transitions in certain places. Useful for testing recovery ideas later, because the data pattern on the corrupted sectors is known.
- MS_Plus_OK_track - OK track from a Microsoft Plus! floppy. This one can be decoded fine with FluxEngine.
- MS_Plus_warped_track - Warped track from a Microsoft Plus! floppy. FluxEngine can't cleanly decode every sector, but both flux-analyze versions can.
