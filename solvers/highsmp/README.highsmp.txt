HiGHS driver for AMPL
=====================

HiGHS (high performance software for linear optimization) is a suite of 
open source serial and parallel solvers for large-scale sparse linear 
programming (LP), mixed-integer programming (MIP), and quadratic programming 
(QP) models. More information can be found in https://highs.dev/.

Normally HiGHS is invoked by AMPL's solve command, which gives the 
invocation

     highs stub -AMPL

in which stub.nl is an AMPL generic output file (possibly written
by "ampl -obstub" or "ampl -ogstub").  After solving the problem,
highs writes a stub.sol file for use by AMPL's solve and solution
commands.  When you run ampl, this all happens automatically if you
give the AMPL commands

     option solver highs;
     solve;

You can control highs either by setting the environment variable
highs_options appropriately (either by using ampl's option command,
or by using the shell's set and export commands before you invoke ampl),
or by passing the options on the command line:

     highs stub [-AMPL] option1=value option2=value ...

You can put one or more (white-space separated) phrases in $highs_options.
To see the possibilities, invoke

     highs -=
