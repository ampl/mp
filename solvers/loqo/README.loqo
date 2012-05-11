Loqo is a version of Bob Vanderbei's LP solver LOQO that can
be used either "stand-alone" or with AMPL's solve command or with
AMPL's -ob and -g options.

Stand-alone invocations of loqo have the form

	loqo stub [-AMPL] [keywd=value ...]

where stub was specified in

	ampl -obstub ...
or
	ampl -ogstub ...

Such an invocation causes loqo to read from stub.nl and, if -AMPL or
"wantsol" or "wantsol=" appears on the command line, to write stub.sol.

You can use AMPL's solve command to invoke loqo and have the solution
loqo computes made available to your AMPL session.  To do this, in the
AMPL session you may need to set the AMPL's solver option to loqo by
giving the command

	option solver loqo;

To see the current solver option, use the command

	option solver;

----------------
Controlling loqo
----------------

Loqo reads key words and values from the environment (shell) variable
loqo_options and from the command line (in which case each key word
must be followed immediately by an = sign).  Execute

	loqo -=

for a summary of the keywords loqo recognizes.

For MPS input, invoke

	loqo -m [stub] [keywd=value ...]

If stub is present, loqo expects to read files stub.spc and stub.mps .

------------------
Sample Invocations
------------------

  If you're using AMPL, just say

	option solver loqo;
	solve;

  If you've executed, say,

	ampl -objunk junk.model junk.data

then you could say

	loqo junk iterlim=30 dual=

to force loqo to solve the dual problem and to have it run for
at most 30 iterations.  Either of the invocations

	loqo_options='iterlim 30 dual' loqo junk
or
	loqo_options='iterlim 30 dual'
	export loqo_options
	loqo junk

would have the same effect; within AMPL, specifying

	option loqo_options 'iterlim 30 dual', solver loqo;
	solve;

would also have this effect.

-----------------------
solve_result_num values
=======================

Here is a table of solve_result_num values that "loqo" can return
to an AMPL session, along with the text that appears in the associated
solve_message.

	Value	Message

	0	optimal solution
	100	suboptimal solution
	201	primal infeasible
	202	dual infeasible
	203	primal and/or dual infeasible
	210	primal infeasible -- inconsistent equations
	301	primal unbounded
	302	dual unbounded
	400	iteration limit
	500	resource limit
	510	??? LOQO bug

-----------------------

Questions about this stuff?  Contact dmg@ampl.com (David M. Gay).
