The solver "x-gurobi" uses Gurobi (a trademark of Gurobi Optimization,
Inc.; see http://www.gurobi.com/) to solve integer, mixed-integer, and
linear programming problems; it is an alternative to the ASL-based 
driver implemented using the mp library (https://github.com/ampl/mp)
for communicating with AMPL and for reformlation of certain types of 
problems.
Normally gurobi is invoked by AMPL's solve command, which gives the 
invocation

     x-gurobi stub -AMPL

in which stub.nl is an AMPL generic output file (possibly written
by "ampl -obstub" or "ampl -ogstub").  After solving the problem,
gurobi writes a stub.sol file for use by ampl's solve and solution
commands.  When you run ampl, this all happens automatically if you
give the AMPL commands

     option solver x-gurobi;
     solve;

You can control gurobi by setting the environment variable x-gurobi_options
appropriately (either by using ampl's option command, or by using the
shell's set and export commands before you invoke ampl).  You can put
one or more (white-space separated) phrases in $x-gurobi_options.  To see
the possibilities, invoke

        x-gurobi -=

----------
INSTALLING
==========

On Linux systems, libgurobi*.so (where the value of "*" depends
on the current version of Gurobi) and the libgurobi.so.* to which
it points need to appear in the current directory when gurobi
itself appears there, or in one of the standard places (specified by
/etc/ld.so.conf on some systems), or in a directory named in
$LD_LIBRARY_PATH.  An alternative is to add a short shell script,
such as

        #!/bin/sh
        LD_LIBRARY_PATH=/usr/local/lib
        export LD_LIBRARY_PATH
        exec /usr/local/bin/gurobix "$@"

to a directory in your usual $PATH (and mark the script executable
with, e.g., "chmod +x gurobi").  The above script assumes that the
true "x-gurobi" binary has been moved to /usr/local/bin/gurobix and that
the libgurobi* files have been moved to /usr/local/lib.

MacOSX systems are similar to Linux systems, but with DYLD_LIBRARY_PATH
in place of LD_LIBRARY_PATH.  Starting 20150225, gurobi binaries
for MacOSX should find the appropriate libgurobi.so.* if it appears
in the same directory as "gurobi".

On MS Windows systems, x-gurobi.exe and the relevant gurobi*.dll must
appear somewhere in your usual search $PATH (or in the current
directory).

If you have questions about or find bugs with this stuff,
please contact:

     AMPL Support
     support@ampl.com
