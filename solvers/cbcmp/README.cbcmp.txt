CBC driver for AMPL
=====================

Cbc (Coin-or branch and cut) is an open-source mixed integer linear 
programming solver written in C++. More info at https://github.com/coin-or/cbc.

Normally cbc is invoked by AMPL's solve command, which gives the 
invocation

     cbc stub -AMPL

in which stub.nl is an AMPL generic output file (possibly written
by "ampl -obstub" or "ampl -ogstub").  After solving the problem,
cbc writes a stub.sol file for use by AMPL's solve and solution
commands.  When you run ampl, this all happens automatically if you
give the AMPL commands

     option solver cbc;
     solve;

You can control cbc either by setting the environment variable
cbc_options appropriately (either by using ampl's option command,
or by using the shell's set and export commands before you invoke ampl),
or by passing the options on the command line:

     cbc stub [-AMPL] option1=value option2=value ...

You can put one or more (white-space separated) phrases in $cbc_options.
To see the possibilities, invoke

     cbc -=
