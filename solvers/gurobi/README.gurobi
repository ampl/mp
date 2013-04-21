The solver "gurobi" uses Gurobi (a trademark of Gurobi Optimization,
Inc.; see http://www.gurobi.com/) to solve integer, mixed-integer, and
linear programming problems.  Normally gurobi is invoked by AMPL's
solve command, which gives the invocation

     gurobi stub -AMPL

in which stub.nl is an AMPL generic output file (possibly written
by "ampl -obstub" or "ampl -ogstub").  After solving the problem,
gurobi writes a stub.sol file for use by ampl's solve and solution
commands.  When you run ampl, this all happens automatically if you
give the AMPL commands

     option solver gurobi;
     solve;

You can control gurobi by setting the environment variable gurobi_options
appropriately (either by using ampl's option command, or by using the
shell's set and export commands before you invoke ampl).  You can put
one or more (white-space separated) phrases in $gurobi_options.  To see
the possibilities, invoke

        gurobi '-='

----------
INSTALLING
==========

On Linux systems, libgurobi*.so (where the value of "*" depends
on the current version of Gurobi) and the libgurobi.so.* to which
it points need to appear in a directory named in $LD_LIBRARY_PATH.
File "gurobi" is a very short and human-readable shell script that
sets LD_LIBRARY_PATH and invokes gurobix (the real binary file).  As
distributed, "gurobi" assumes that gurobix is in /usr/local/bin and
libgurobi*.so is in /usr/local/lib.  You can easiy change these
assumptions by editing the gurobi shell script.  This script needs
to appear somewhere in your usual search $PATH.

MacOSX systems are similar to Linux systems, but with DYLD_LIBRARY_PATH
in place of LD_LIBRARY_PATH.  The gurobi script assumes the appropriate
libgurobi.so.* is in /usr/local/lib/x86_64 -- but you can easily
change this detail.

On MS Windows systems, gurobi.exe and the relevant gurobi*.dll must
appear somewhere in your usual search $PATH (or in the current
directory).


-----------------------
solve_result_num values
=======================

Here is a table of solve_result_num values that "gurobi" can return
to an AMPL session, along with an indication of the text that appears
in the associated solve_message.

        Value   Message

          0     optimal solution
        100     suboptimal: could not satisfy optimaliter tolerances
        200     infeasible [IIS computation not attempted]
        201     infeasible [IIS returned]
        202     infeasible [IIS finder failed]
        203     infeasible; .dunbdd returned [IIS computation not attempted]
        204     infeasible; .dunbdd returned [IIS also returned]
        205     infeasible; .dunbdd returned [IIS finder failed]
        300     unbounded
        301     unbounded [unbounded or infeasible; IIS finder failed]
        302     unbounded; .unbdd returned
        303     unbounded; .unbdd returned [IIS finder failed]
        400     objective cutoff
        401     iteration limit
        402     node limit
        403     time limit
        404     solution limit
        500     Could not create the gurobi environment
        501     Gurobi call failed [message gives routine name]
        502     misc. failure [message gives details]
        503     Bad $gurobi_options
        504     Surprise VBasis[...] = ...
        505     Surprise CBasis[...] = ...
        506     Gurobi set/get parameter failed [message gives more details]
        510     cannot open logfile (specified in $gurobi_options or command line)
        511     cannot open paramfile (specified in $gurobi_options or command line)
        512     missing value in paramfile
        513     extra text in paramfile
        514     invalid parameter name in paramfile
        520     numeric error
        521     nonlinear objective
        522     nonlinear constraint
        523     quadratic objective involving division by 0
        524     indefinite quadratic objective or constraint
        525     quadratic constraint involving division by 0
        530     could not open serverlic file
        531     error in serverlic file
        532     error while tuning
        567     complementarity constraint
        600     interrupted
        601     could not talk to Gurobi compute server
        602     job rejected by Gurobi compute server
        603     no license for specified gurobi compute server
        604     surprise return while trying to use Gurobi compute server

Values 521-524 only arise in Gurobi versions >= 4.0.
Values 203-205 and 302-303 only arise in Gurobi versions >= 4.5.

*************************

If you have questions about or find bugs with this stuff,
please contact:

     David M. Gay
     dmg@ampl.com
