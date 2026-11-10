MP2NL driver for AMPL
====================

MP2NL is a metadriver to run existing NL solvers through
the MP library. See the MP Modeling Guide, Tools & details > Meta-driver MP2NL:
https://mp.ampl.com/modeling-tools.html#mp2nl

Normally MP2NL is invoked by AMPL's solve command, which gives the
invocation

     mp2nl stub -AMPL

in which stub.nl is an AMPL generic output file (possibly written
by "ampl -obstub" or "ampl -ogstub").  After solving the problem,
mp2nl writes a stub.sol file for use by AMPL's solve and solution
commands.  When you run ampl, this all happens automatically if you
give the AMPL commands

     option solver mp2nl;
     solve;

You can control mp2nl either by setting the environment variable
mp2nl_options appropriately (either by using ampl's option command,
or by using the shell's set and export commands before you invoke ampl),
or by passing the options on the command line:

     mp2nl stub [-AMPL] option1=value option2=value ...

You can put one or more (white-space separated) phrases in $mp2nl_options.
To see the possibilities, invoke

     mp2nl -=
