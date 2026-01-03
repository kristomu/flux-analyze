## Parity calculation

The dynamic programming algorithm works by checking every possible evolution from
one pulse transition to another (implicitly with caching), and then picking the
path with the least cost/lowest penalty.

If we want to penalize MFM errors, we need to keep track of the MFM decoding
state as well. This is tricky because whether a particular sequence of bits
(produced by the delay between two flux transitions) is legal depends on
adjacent pairs, but a flux delay may produce an odd number of bits.

Once we know the clock, each flux delay (timing between two transitions) produces
the following bits:

| Number of half-clocks	| Bit sequence |
| --------------------- | ------------ |
| 1                     | 0 1          |
| 2                     | 0 0 1        |
| 3                     | 0 0 0 1      |

(The MFM spec disallows shorter or longer intervals.)

The valid transitions for bit sequences, and what data bits they decode to, are:

| bit pair | valid if last data bit is | produces data bit |
| -------- | ------------------------- | ----------------- |
| 0 1      | (Always valid)            | 1                 |
| 0 0      | 1                         | 0                 |
| 1 0      | 0                         | 0                 |

### Verifying if an invalid transition happened

#### Parity

Our function should take the current pulse length as well as the previous pulse
length and its parity, and return whether the combination is allowed or not
(assuming that everything that can be identified by running the function on the
previous and its precedent were valid).

"Parity is odd" means that the number of bits seen (from all pulse trains so
far), up to and including the bits produced by the last ("previous") pulse,
is odd.

That is, a function call of the type f(2, 3, odd), indicates that we're dealing
with a situation like this:

        1     k-2 k-1 k   k+1 k+2 k+3 k+4
          ... 0   0   1 | 0   0   0   1

with k % 2 == 1, and the bar separating the bits produced by prior pulse delays
and the bits produced by the current ("present") pulse delay.

#### Filling in zeroes (even parity)

When decoding the validity of a sequence by the MFM rules, if the parity is
even, we're guaranteed to have at least

               ... 0 1 | present pulse goes here

This is because the previous pulse delay must at least have been one half-clock
long, and thus have at least one zero before its trailing edge.

By analogous reasoning, we know that the bit immediately after the present pulse's trailing edge must be a zero:

				past      present   future
				... 0 1 | ..... 1 | 0 ...

#### When is even parity followed by a new pulse valid

Thus, when we're dealing with even parity, we know that the last bit pair before
the present pulse is 01, which decodes to a one and is always valid by the MFM
table, no matter what bit pairs came before it. All we have to do is determine
if the bits produced by the different length pulse delays are valid, given that
the last bit pair produced a one.

Consider them exhaustively, where "present" is the length of the current pulse delay:

                ... 01 | 01         <- present = 1
                ... 01 | 00 1x      <- present = 2
                ... 01 | 00 01      <- present = 3

present = 1 is always valid because "01" bit pairs are always valid.
present = 3 is valid because 00 is OK as long as the previous bit was a one,
which it is because the previous pair was 01.

present = 2 is a little trickier. But we know that the future unknown must be a
zero, hence this is actually

				... 01 | 00 10

which is also valid because 01 unconditionally decodes to a 1, 00 decodes to a 0
on the condition that the previous decoded bit is 1, and 10 decodes to a 0 on
the condition that the previous decoded bit is 0. All of this is the case.

So even parities are always valid.

#### When is odd parity followed by a new pulse valid

Odd parity is more subtle because we now have the trailing edge from the past
pulse as the first bit:

				... 1|0 1x                  10 10           <- present = 1
				... 1|0 01                  10 01           <- present = 2
				... 1|0 00 1x               10 00 10        <- present = 3

But if our function checks the validity of every bit pair produced by the
present pulse delay, then we can inductively assume that the function, called
on the previous pulse delay, would have confirmed that the 10 that we're only
partly seeing here was indeed valid. Or that it would have counted an error that
we don't need to count here.

By that assumption, present = 1 is valid because the subsequent 10 decodes on a
0 contingent on the prior bit being a 0, which is the case because the prior bit
pair was 10, which also decodes to a zero.

Present = 2 is also valid because 01 is unconditionally valid.

Present = 3 is not. 10 decodes to 0, but the 00 that follows it needs the last
decoded bit to be a one.

### All that to a simple rule

So if all that's correct, then everything is valid *except*

	(past length, present length, parity) = (x, 3, odd),

for any x. (Of course, out-of-bounds pulse lengths for past or present are all
disallowed.)

And that's how very complex reasoning results in a one-line function.