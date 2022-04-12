COPT Solver
===========

The solver "copt" uses COPT (a trademark of Cardinal Operations; 
see https://www.shanshu.ai/copt) to solve integer, mixed-integer, and
linear programming problems.  Normally the solver is invoked by AMPL's
solve command, which gives the invocation

     copt stub -AMPL

in which stub.nl is an AMPL generic output file (possibly written
by "ampl -obstub" or "ampl -ogstub").  After solving the problem,
copt writes a stub.sol file for use by ampl's solve and solution
commands.  When you run ampl, this all happens automatically if you
give the AMPL commands

     option solver copt;
     solve;

You can control copt by setting the environment variable copt_options
appropriately (either by using ampl's option command, or by using the
shell's set and export commands before you invoke ampl).  You can put
one or more (white-space separated) phrases in $copt_options.  To see
the possibilities, invoke

     copt -=

Installing
----------

### Linux

On Linux systems, libcopt.so to appear in the directory where the 
binary `copt` appears, or in one of the standard places (specified by
/etc/ld.so.conf on some systems), or in a directory named in
$LD_LIBRARY_PATH.  An alternative is to add a short shell script,
such as

        #!/bin/sh
        LD_LIBRARY_PATH=/usr/local/lib
        export LD_LIBRARY_PATH
        exec /usr/local/bin/coptx "$@"

to a directory in your usual $PATH (and mark the script executable
with, e.g., "chmod +x copt").  The above script assumes that the
true `copt` binary has been moved to /usr/local/bin/coptx and that
the libcopt.so file has been moved to /usr/local/lib.

### Mac

MacOSX systems are similar to Linux systems, but with DYLD_LIBRARY_PATH
in place of LD_LIBRARY_PATH.

### Windows

On MS Windows systems, copt.exe and the relevant copt.dll must
appear somewhere in your usual search `$PATH` (or in the current
directory).

If you have questions or comments about this please contact:

     AMPL Support
     support@ampl.com
