SCIP driver for AMPL
====================

SCIP (Solving Constraint Integer Programs) is currently one of the
fastest non-commercial solvers for mixed integer programming (MIP)
and mixed integer nonlinear programming (MINLP). It is also a framework
for constraint integer programming and branch-cut-and-price.
More information can be found in https://scipopt.org/#scipoptsuite.

Normally SCIP is invoked by AMPL's solve command, which gives the
invocation

     scip stub -AMPL

in which stub.nl is an AMPL generic output file (possibly written
by "ampl -obstub" or "ampl -ogstub").  After solving the problem,
scip writes a stub.sol file for use by AMPL's solve and solution
commands.  When you run ampl, this all happens automatically if you
give the AMPL commands

     option solver scip;
     solve;

You can control scip either by setting the environment variable
scip_options appropriately (either by using ampl's option command,
or by using the shell's set and export commands before you invoke ampl),
or by passing the options on the command line:

     scip stub [-AMPL] option1=value option2=value ...

You can put one or more (white-space separated) phrases in $scip_options.
To see the possibilities, invoke

     scip -=
