# Contents

This directory contains various MFM decoders. These decode flux delay streams
(timings between individual flux transitions) to MFM trains (binary sequences
that can be further decoded by using the MFM spec).

# Variable-offset decoders

Most of the variable clock decoders in this directory uses a model where
certain multiples of a slowly varying base clock represent short, medium, or
long MFM transitions; more specifically, that short transitions are two
half-clocks, medium are three, and long are four.

With an uncorrupted MFM-encoded floppy, you would expect to see a delay
distribution with peaks at 2x, 3x, and 4x some fundamental. But that's not
what I see; instead they're slightly off. One of my floppies has roughly
the following means:

	1x clock: 22.99
	1.5x clock: 35.624
	2x clock: 47.520

which suggests that either different sectors have been written using
differently calibrated drives, or that there's a constant offset (perhaps due
to hardware delays, imperfect clocks, etc). The variable-offset decoders
were my attempt to model such an "intercept and slope" state space, so that
instead of

	number of MFM nonzeroes = round(x * 2/clock) - 1

you'd have

	number of MFM nonzeroes = round(x * 2/clock - offset) - 1

with the <tt>{clock, offset}</tt> state slowly drifting with time.

The variable_offset.cc file contains one such decoder based on gradient descent.
This seems to work, but I think it's mostly by accident.

I tried to implement another decoder based on recursive least squares with
forgetting, but either I didn't get it implemented correctly, or there's a problem
with the state space model of <tt>round(x * 2/clock - offset)</tt>, because it
very much didn't work.

## The problem with (ax - b)

With the EWMA decoders, the base clock is updated by choosing <tt>clock_new</tt> so
that the current observed flux delay is bang in the middle of its band if clock_new
were the real clock. Then this clock is mixed in with the current estimate using a
forgetting factor <tt>alpha</tt>:

	proposed_mean_clock = (1-alpha) * mean_clock + alpha * new_clock;

Let's ignore the obvious problem that expanding this logic to (ax-b) requires solving
for two unknowns instead of one, because RLS can handle that in a principled way as
long as not every MFM transition is of the same type.

The more serious problem (I think) is that the optimizer can hack its reward. If the
objective function penalty is something like

	lambda^(k-i) * ||y_i - (ax_i - b)||^2

and

	y_i = round(a_i * x_i - b_i) - 1

then there's an incentive to make b large and a small, so that the inner expression
<tt>a_i * x_i - b_i</tt> becomes near-constant. Then the error term

	||y_i - (ax_i - b)||^2

goes to zero no matter what the flux delay actually is. So the x_i -> y_i calculation
instead has to be something like

	number of MFM nonzeroes = round( (x - offset) * 2/clock) - 1

i.e. the offset must be in terms of flux delay units, *not* MFM nonzeroes. Doing so
prevents the reward hack because decreasing the clock slope value means that the
offset must be blown up to get the same effect. But handling this state model
requires nonlinear optimization, which I haven't done yet.

And so I've omitted the RLS decoder for my repo for two reasons: first, I'm not sure
that I implemented it correctly, and second, I don't think it can actually work. I may
try to do nonlinear RLS, but it's not my field of expertise.