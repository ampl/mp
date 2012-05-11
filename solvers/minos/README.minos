Minos is a solver based on MINOS versions 5.51 that can be used
either "stand-alone" or with AMPL's -ob and -og options to solve
linear and nonlinear problems expressed in AMPL.

To use minos with AMPL, you have several options.  You can invoke
it within an AMPL session by saying

	solve;

or, if $solver is not already minos,

	option solver minos;
	solve;

Minos understands AMPL's -ob and -og output formats; you can thus
use stand-alone invocations like

	ampl -obfoo foo.mod foo.dat
	minos foo [assignments]

With no stub argument, minos tries to read a SPECS file on stdin,
followed (unless the SPECS files dictates otherwise) by an MPS file,
so you also can invoke

	cat foo.specs foo.mps | minos [assignment ...]

Invocation, in general, is

	minos [options] [stub [-AMPL]] [assignment ...]

where stub is from `ampl -obstub` or `ampl -ogstub`.
Assignments have the form spec_phrase=value or n=filename
(where n is a 1 or 2 digit Fortran unit number, presumably one
mentioned in a spec_phrase=value assignment or the SPECS file).
No spaces may appear in either form of assignment; spec_phrases
are phrases that can appear in a SPECS file, with _ (underscore)
substituted for blank.  An assignment n=filename attaches filename
to Fortran unit n.  Assignments can also appear in the environment
variable $minos_options; minos first reads the SPECS file (if any),
then $minos_options, then the command-line assignments.  For example

	minos foo backup_basis_file=2 2=zip new_basis_file=3 3=zap

will make minos behave as though it read a SPECS file containing

	BACKUP BASIS FILE 2
	NEW BASIS FILE 3

and will connect files zip and zap with the Fortran units 2 and 3.
(It's unfortunate that you must explicitly mention Fortran unit
numbers, but this is not onerous.)  The default file name for Fortran
unit u is fort.u .

Specifying objno=n in $minos_options or the command line is the
same as specifying problem_number=n-1 (i.e., objno=2 is the
same as problem_number=1), and objno=0 means "ignore the objectve;
just seek a feasible point".

Similar comments apply to a.out's you create to solve nonlinear
problems by calling MINOS.

The argument -AMPL causes minos to emit a one-line banner; when
$solver has its default value (minos), AMPL's solve command invokes

	minos stub -AMPL

If a stub is present, minos tries to write the computed solution to
stub.sol unless the -s option was specified.  Execute

	minos '-?'

for a summary of other options.

Minos determines how much scratch storage it will need from the
stub.nl or the SPECS file and obtains it from malloc.  For
details not given here, see the "MINOS 5.1 User's Guide"
(by B. A. Murtagh and M. A. Saunders, Tech. Rep. SOL 83-20R,
http://www.stanford.edu/group/SOL/guides/minos55.pdf ).
In particular, the User's Guide gives full details about the "specs"
and MPS files mentioned below.
See also http://www.ampl.com/BOOKLETS/ampl-minos.pdf .

Default file assignments (Fortran units) are:
	5 = stdin (for "specs" and MPS files when no basename is given)
	6 = stdout (for MINOS summary file)
	7 = MINOS summary file (only if an assignment 7=filename
		appears on the command line or in $minos_options)

For invocations from AMPL's solve command or of the form

	minos stub ...

(where stub.nl is from AMPL's -ob or -og output options), you can use
outlev= to control the amount and kind of output:
	outlev=0	no chatter on stdout
	outlev=1	only report options on stdout
	outlev=2	summary file on stdout
	outlev=3	log file on stdout, no solution
	outlev=4	log file, including solution, on stdout

Source for MINOS is available from

	Stanford Business Software
	Phone:	+1 415-962-8719
	Fax:	+1 415-962-1869
	2680 Bayshore Parkway, Suite 304
	Mountain View, CA 94043


-----------------------
solve_result_num values
=======================

Here is a table of solve_result_num values that "minos" can return
to an AMPL session, along with the text that appears in the associated
solve_message.

	Value	Message

	0	optimal solution found
	100	optimal solution found?  Optimality
		tests satisfied, but reduced gradient is large
	200	infeasible problem
	200	infeasible problem (or bad starting guess)
	201	numerical error: the general constraints
		cannot be satisfied accurately
	300	unbounded (or badly scaled) problem
	400	too many iterations
	401	too many major iterations
	500	the objective has not changed for the last %ld iterations
	501	the current point cannot be improved
	510	singular basis after several factorization attempts
	520	the superbasics limit (%ld) is too small
	521	error evaluating nonlinear expressions
	522	not enough storage for the basis factors.
		Try rerunning with workspace_(total)=nnn in $minos_options 
	530	incorrect gradients from funobj
	531	incorrect gradients from funcon
	532	cannot find superbasic to replace basic variable
	533	basis factorization requested twice in a row
	534	error in basis package
	535	input basis had wrong dimensions
	536	unexpected return code (nnn)
	540	solution aborted

-----------------------

Questions about this stuff? Contact dmg@ampl.com (David M. Gay).
