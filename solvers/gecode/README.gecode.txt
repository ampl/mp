The solver Gecode for AMPL uses Gecode (Generic Constraint Development
Environment) to solve constraint-based AMPL models (https://www.gecode.org/).
Normally the driver is invoked by AMPL's solve command, which gives the 
invocation

     gecode stub -AMPL

in which stub.nl is an AMPL generic output file (possibly written
by "ampl -obstub" or "ampl -ogstub").  After solving the problem,
gecode writes a stub.sol file for use by ampl's solve and solution
commands.  When you run ampl, this all happens automatically if you
give the AMPL commands

     option solver gecode;
     solve;

You can control gecode by setting the environment variable gecode_options
appropriately (either by using AMPL's option command, or by using the
shell's set and export commands before you invoke AMPL).  You can put
one or more (white-space separated) phrases in $gecode_options.  To see
the possibilities, invoke

     gecode -=
