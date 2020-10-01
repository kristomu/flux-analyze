# doesn't work

This directory contains approaches that don't work, but may later be useful as
the basis for something that does.

Currently:

## mfm_dp_dewarp.cc: 

An optimal penalty dewarping approach that uses dynamic programming, tailored
for MFM. It chooses the band boundaries for each timestep (analogous to
clustering centerpoints) so as to minimize the change from one timestep to
another (as measured by Lp norm distance), under the constraint that there must
be no MFM errors. 

The dynamic programming algorithm is basically shortest path on a DAG with edge
costs equal to the one-step Lp norm difference between band boundaries. I make
refereces to a "triangle inequality trick" that you really just need to find the
valid past node that is the closest neighbor (band boundary wise) to the present
node under consideration, but I haven't implemented that yet. (Proof left as an
exercise for the reader...)

In any case, this fails because "no MFM errors" is too weak a constraint, and
so instead of the decoder directly tracking the curve of the warping, you get a
lazier (lower penalty) curve that doesn't track the warp but is just nonlinear
enough to avoid getting any MFM decoding errors.
