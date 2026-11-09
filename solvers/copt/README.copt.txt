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

COPT 9 Beta and compatibility with COPT 8
----------------------------------------

The driver supports COPT 8.0.6 and COPT 9 Beta. Use the matching COPT
runtime for your driver version. The options below require COPT 9;
existing vendor defaults are retained.

  Option / COPT 9 choice                Meaning
  alg:method=7                         Hybrid simplex/barrier/PDLP selection
  alg:method=4 lp:concurrentmode=4      Deterministic concurrent CPU LP
  lp:barorder=3                        Nested dissection ND3
  lp:barorder=4..6                     Concurrent ordering strategies 1..3
  mip:startmode=3                      Accept all hierarchical multiobj starts
  mip:concurrent=N (concurrentmip=N)    Number of concurrent MIP solves

For mip:concurrent, -1 is the vendor default (automatic), 0 or 1 disables
concurrent MIP, and 2..64 requests multiple solves with different settings.
This differs from tech:miptasks, the task count within a MIP solve.
COPT does not support concurrent MIP with user callbacks or lazy constraints;
the vendor may fall back to a single solve in those situations. A request
for concurrent MIP is not a guarantee of parallel worker utilization.

Explicit legacy choices remain available: alg:method=1 (dual simplex),
alg:method=2 (barrier), alg:method=5 (simplex/barrier selection),
alg:method=4 lp:concurrentmode=0 (CPU concurrent LP), lp:barorder=-1..2,
and mip:startmode=-1..2. To disable concurrent MIP in a COPT 9 build, use
mip:concurrent=0 or 1. With COPT 8, mip:concurrent is unavailable and
the new numeric choices are rejected.
These controls preserve algorithm selection, not identical numerical
results or performance across vendor releases.

Per-worker parameter configuration for concurrent MIP is not exposed by
the driver. Use tech:optionnative and parameter-file options for ordinary
model parameters.

GPU execution requires a supported Windows/Linux CUDA environment.
General NLP is solved locally; general MINLP and native complementarity
are not supported by this interface.

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
