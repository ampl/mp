The solver "xpressmp" uses XPRESS (see https://www.fico.com/) 
to solve integer, mixed-integer, and linear programming problems; 
it is an alternative to the ASL-based driver implemented using the 
mp library (https://github.com/ampl/mp) for communicating with AMPL 
and for reformlation of certain types of problems.
Normally xpress is invoked by AMPL's solve command, which gives the 
invocation

     xpressmp stub -AMPL

in which stub.nl is an AMPL generic output file (possibly written
by "ampl -obstub" or "ampl -ogstub").  After solving the problem,
xpressmp writes a stub.sol file for use by ampl's solve and solution
commands.  When you run ampl, this all happens automatically if you
give the AMPL commands

     option solver xpressmp;
     solve;

You can control gurobi by setting the environment variable xpress_options
appropriately (either by using ampl's option command, or by using the
shell's set and export commands before you invoke ampl).  You can put
one or more (white-space separated) phrases in $xpress_options.  To see
the possibilities, invoke

        xpressmp -=
