GCG driver for AMPL
===================

GCG is a generic decomposition solver for mixed integer linear
programs that extends the SCIP (Solving Constraint Integer Programs)
framework. It finds structures in models that can be used to apply a
Dantzig-Wolfe reformulation or Benders decomposition. Decompositions
can also be user-given, and explored and evaluated manually. For a
Dantzig-Wolfe reformulated model, a branch-price-and-cut algorithm
solves the problem, which features primal heuristics, generic and
specific pricing solvers, branching rules, dual variable stabilization,
cutting planes, etc. Like SCIP, also GCG can be used as a framework
and extended to suit one's needs.

Normally GCG is invoked by AMPL's solve command, which gives the
invocation

     gcg stub -AMPL

in which stub.nl is an AMPL generic output file (possibly written
by "ampl -obstub" or "ampl -ogstub").  After solving the problem,
gcg writes a stub.sol file for use by AMPL's solve and solution
commands.  When you run ampl, this all happens automatically if you
give the AMPL commands

     option solver gcg;
     solve;

You can control gcg either by setting the environment variable
gcg_options appropriately (either by using ampl's option command,
or by using the shell's set and export commands before you invoke ampl),
or by passing the options on the command line:

     gcg stub [-AMPL] option1=value option2=value ...

You can put one or more (white-space separated) phrases in $gcg_options.
To see the possibilities, invoke

     gcg -=
