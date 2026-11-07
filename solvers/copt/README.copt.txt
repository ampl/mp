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

COPT 9 and compatibility with COPT 8
-----------------------------------

The driver can be built with COPT 8.0.6 or COPT 9. New COPT 9 option
choices are enabled at compilation when using COPT 9 headers. Use the
matching vendor runtime for each build; this is not runtime SDK switching.
The vendor's existing defaults are retained.

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
mip:concurrent=0 or 1. COPT 8 builds do not register mip:concurrent and
reject the new numeric choices using the vendor parameter bounds.
These controls preserve algorithm selection, not identical numerical
results or performance across vendor releases.

COPT 9 also supplies APIs for reading per-worker parameter settings. The
driver does not yet expose these APIs: the Beta manual does not specify
the file/string format, and a successful API return alone does not confirm
that settings were applied. Existing tech:optionnative and parameter-file
options continue to control ordinary model parameters.

The COPT 9 release notes describe new vendor BLAS/LAPACK implementations,
improved barrier ordering, SDP, MIP, QP and GPU PDLP performance. These
changes require no new model API. GPU execution still requires a supported
Windows/Linux CUDA environment; macOS tests exercise CPU paths.
General NLP remains a local solver; COPT 9 does not introduce a documented
general MINLP or native complementarity interface.

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
