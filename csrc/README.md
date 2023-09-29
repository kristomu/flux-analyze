# flux-analyze
Analysis and recovery of IBM MFM floppy image data, C++

This tool can now be used to create disk images of 1.44M or higher format
IBM floppies, although it is a little rough around the edges.

To do so, first build (using cmake), then invoke flux-analyze followed by
one or more .flux image files (old FluxEngine format). Specifying multiple
images lets it combine flux records of the same disk from multiple floppies;
any sector that can be recovered using one of the flux images will be present
in the output.

The output files are output.img and output.mask.

output.img is the file system image. Bad sectors are represented as nulls.
output.mask is a mask that shows what bytes were recovered. If a byte at
position x is 0xff, then that byte was recovered in output.img; if it is 0,
then flux-analyze couldn't recover that sector.

Note: The tool is still very rough and takes a lot of time. It will also
output a lot of progress information. Usually, only the summary at the end
is important: it'll list what sectors it was unable to recover ("Couldn't find")
or entire tracks ("Couldn't find any sector").
