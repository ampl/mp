#include <climits>
#include <cfloat>

#include "mp/common.h"
#include "mp/error.h"
#include "mp/backend-std.h"

extern "C" {
  #include "xprs.h"
  #include "xslp.h"
}


namespace mp {

/// A mix-in class to add Xpress parameters.
/// Translated from '../mp/solvers/xpress/optprm.h'
/// on Thu Nov 20 15:22:33 2025
///
template <class Impl>
class CompiledOptimizerOptions {
public:
  /// Add up to 337 'Optimizer' parameters
  void AddOptimizerOptions() {

#ifdef XPRS_ALGAFTERCROSSOVER
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_algaftercrossover XPRS_ALGAFTERCROSSOVER",
      "The algorithm to be used for the final clean up step after the crossover."
      "\n\n"
      "Values (default: 1):\n"
      "\n- (1)  Automatically determined."
      "\n- (2)  Dual simplex."
      "\n- (3)  Primal simplex."
      "\n- (4)  Concurrent.",
      XPRS_ALGAFTERCROSSOVER, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_ALGAFTERCROSSOVER

#ifdef XPRS_ALGAFTERNETWORK
    MPD( AddSolverOption_MergeDuplicates("lp:xprs_algafternetwork XPRS_ALGAFTERNETWORK",
      "The algorithm to be used for the clean up step after the network simplex solver."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatically determined."
      "\n- (2)  Dual simplex."
      "\n- (3)  Primal simplex.",
      XPRS_ALGAFTERNETWORK, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_ALGAFTERNETWORK

#ifdef XPRS_ALTERNATIVEREDCOSTS
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_alternativeredcosts XPRS_ALTERNATIVEREDCOSTS",
      "Controls aggressiveness of searching for alternative reduced cost"
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  The solver decides if searching for alternative reduced cost is beneficial or not. This is the default setting."
      "\n- (0)  Searching for alternative reduced cost is disabled."
      "\n- (1)  Searching for alternative reduced cost is enabled.",
      XPRS_ALTERNATIVEREDCOSTS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_ALTERNATIVEREDCOSTS

#ifdef XPRS_AUTOCUTTING
    MPD( AddSolverOption_MergeDuplicates("cut:xprs_autocutting XPRS_AUTOCUTTING",
      "Should the Optimizer automatically decide whether to generate cutting planes at local nodes in the tree or not? If the CUTFREQ control is set, no automatic selection will be made and local cutting will be enabled. "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatic."
      "\n- (0)  Disabled."
      "\n- (1)  Enabled. ",
      XPRS_AUTOCUTTING, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_AUTOCUTTING

#ifdef XPRS_AUTOPERTURB
    MPD( AddSolverOption_MergeDuplicates("lp:xprs_autoperturb XPRS_AUTOPERTURB",
      "Simplex: This indicates whether automatic perturbation is performed. If this is set to 1, the problem will be perturbed whenever the simplex method encounters an excessive number of degenerate pivot steps, thus preventing the Optimizer being hindered by degeneracies."
      "\n\n"
      "Values (default: 1):\n"
      "\n- (0)  No perturbation performed."
      "\n- (1)  Automatic perturbation is performed.",
      XPRS_AUTOPERTURB, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_AUTOPERTURB

#ifdef XPRS_AUTOSCALING
    MPD( AddSolverOption_MergeDuplicates("num:xprs_autoscaling XPRS_AUTOSCALING",
      "Whether the Optimizer should automatically select between different scaling algorithms. If the SCALING control is set, no automatic scaling will be applied. "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatic."
      "\n- (0)  Disabled."
      "\n- (1)  Cautious strategy. Non-standard scaling will only be selected if it appears to be clearly superior."
      "\n- (2)  Moderate strategy."
      "\n- (3)  Aggressive strategy. Standard scaling will only be selected if it appears to be clearly superior.",
      XPRS_AUTOSCALING, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_AUTOSCALING

#ifdef XPRS_BACKGROUNDMAXTHREADS
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_backgroundmaxthreads XPRS_BACKGROUNDMAXTHREADS",
      "     Limit the number of threads to use in background jobs (for example in         parallel to the root cut loop).   "
      "\n\nDefault: -1, let Xpress decide.",
      XPRS_BACKGROUNDMAXTHREADS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BACKGROUNDMAXTHREADS

#ifdef XPRS_BACKGROUNDSELECT
    MPD( AddSolverOption_MergeDuplicates("heur:xprs_backgroundselect XPRS_BACKGROUNDSELECT",
      "     Bit-vector control (see Section Bit-vector controls) to select which tasks to run in background jobs (for example in parallel to         the root cut loop).   "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (0)  Feasibility jump heuristic."
      "\n- (1)  Fast branch-and-bound heuristic."
      "\n- (2)  Same as bit 1 but with some additional heuristics         enabled."
      "\n- (3)  Fix-propagate-repair heuristic.",
      XPRS_BACKGROUNDSELECT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BACKGROUNDSELECT

#ifdef XPRS_BACKTRACK
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_backtrack XPRS_BACKTRACK",
      "Branch and Bound: Specifies how to select the next node to work on when a full backtrack is performed."
      "\n\n"
      "Values (default: 3):\n"
      "\n- (-1)  Automatically determined."
      "\n- (1)  Unused."
      "\n- (2)  Select the node with the best estimated solution."
      "\n- (3)  Select the node with the best bound on the solution."
      "\n- (4)  Select the deepest node in the search tree (equivalent to depth-first search)."
      "\n- (5)  Select the highest node in the search tree (equivalent to breadth-first search)."
      "\n- (6)  Select the earliest node created."
      "\n- (7)  Select the latest node created."
      "\n- (8)  Select a node randomly."
      "\n- (9)  Select the node whose LP relaxation contains the fewest number of infeasible MIP entities."
      "\n- (10)  Combination of 2 and 9."
      "\n- (11)  Combination of 2 and 4."
      "\n- (12)  Combination of 3 and 4.",
      XPRS_BACKTRACK, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BACKTRACK

#ifdef XPRS_BACKTRACKTIE
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_backtracktie XPRS_BACKTRACKTIE",
      "Branch and Bound: Specifies how to break ties when selecting the next node to work on when a full backtrack is performed. The options are the same as for the BACKTRACK control."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Default selection."
      "\n- (1)  Unused."
      "\n- (2)  Select the node with the best estimated solution."
      "\n- (3)  Select the node with the best bound on the solution."
      "\n- (4)  Select the deepest node in the search tree (equivalent to depth-first search)."
      "\n- (5)  Select the highest node in the search tree (equivalent to breadth-first search)."
      "\n- (6)  Select the earliest node created."
      "\n- (7)  Select the latest node created."
      "\n- (8)  Select a node randomly."
      "\n- (9)  Select the node whose LP relaxation contains the fewest number of infeasible MIP entities."
      "\n- (10)  Combination of 2 and 9."
      "\n- (11)  Combination of 2 and 4."
      "\n- (12)  Combination of 3 and 4.",
      XPRS_BACKTRACKTIE, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BACKTRACKTIE

#ifdef XPRS_BARALG
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_baralg XPRS_BARALG",
      "This control determines which barrier algorithm is used to solve the problem. Notably, this is also the control to enable the primal-dual hybrid gradient algorithm. "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Determined automatically."
      "\n- (0)  Unused."
      "\n- (1)  Use the infeasible-start barrier algorithm."
      "\n- (2)  Use the homogeneous self-dual barrier algorithm."
      "\n- (3)  Start with 2 and optionally switch to 1 during the execution."
      "\n- (4)  Use the hybrid gradient algorithm.",
      XPRS_BARALG, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BARALG

#ifdef XPRS_BARCORES
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_barcores XPRS_BARCORES",
      "If set to a positive integer it determines the number of physical CPU cores assumed to be present in the system by the barrier and hybrid gradient algorithms. If the value is set to the default value (-1), Xpress will automatically detect the number of cores."
      "\n\nDefault: -1(automatically detected)",
      XPRS_BARCORES, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BARCORES

#ifdef XPRS_BARCRASH
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_barcrash XPRS_BARCRASH",
      "Newton barrier and hybrid gradient: This determines the type of crash used for the crossover. During the crash procedure, an initial basis is determined which attempts to speed up the crossover. A good choice at this stage will significantly reduce the number of iterations required to crossover to an optimal solution. The possible values increase proportionally to their time-consumption."
      "\n\n"
      "Values (default: 4):\n"
      "\n- (0)  Turns off all crash procedures."
      "\n- (1-6)  Available strategies with 1 being conservative and 6 being aggressive.",
      XPRS_BARCRASH, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BARCRASH

#ifdef XPRS_BARDUALSTOP
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_bardualstop XPRS_BARDUALSTOP",
      "Newton barrier and hybrid gradient: This is a convergence parameter, representing the tolerance for dual infeasibilities. If the difference between the constraints and their bounds in the dual problem falls below this tolerance in absolute value, optimization will stop and the current solution will be returned."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  The default value is determined automatically based on the problem size, structure and algorithm choice."
      "\n- (>=0)  The tolerance for dual infeasibilities.",
      XPRS_BARDUALSTOP, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_BARDUALSTOP

#ifdef XPRS_BARFAILITERLIMIT
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_barfailiterlimit XPRS_BARFAILITERLIMIT",
      "Newton barrier: The maximum number of consecutive iterations that fail to improve the solution in the barrier algorithm. "
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  Determined automatically"
      "\n- (>0)  Maximum number of consecutive barrier iterations allowed without progress.",
      XPRS_BARFAILITERLIMIT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BARFAILITERLIMIT

#ifdef XPRS_BARFREESCALE
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_barfreescale XPRS_BARFREESCALE",
      "Defines how the barrier algorithm scales free variables."
      "\n\nDefault: 1e-6",
      XPRS_BARFREESCALE, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_BARFREESCALE

#ifdef XPRS_BARGAPSTOP
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_bargapstop XPRS_BARGAPSTOP",
      "Newton barrier and hybrid gradient: This is a convergence parameter, representing the tolerance for the relative duality gap. When the difference between the primal and dual objective function values falls below this tolerance, the Optimizer determines that the optimal solution has been found. "
      "\n\n"
      "Values (default: 0 ):\n"
      "\n- (0)  The default value is determined automatically based on the problem size, structure and algorithm choice."
      "\n- (>=0)  The tolerance for the relative duality gap.",
      XPRS_BARGAPSTOP, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_BARGAPSTOP

#ifdef XPRS_BARGAPTARGET
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_bargaptarget XPRS_BARGAPTARGET",
      "Newton barrier: The target tolerance for the relative duality gap. The barrier algorithm will keep iterating until either BARGAPTARGET is satisfied or until no further improvements are possible. In the latter case, if BARGAPSTOP is satisfied, it will declare the problem optimal. "
      "\n\n"
      "Values (default: 0 ):\n"
      "\n- (0)  The default value is determined automatically based on the problem size, structure and algorithm choice."
      "\n- (>=0)  The target tolerance for the relative duality gap.",
      XPRS_BARGAPTARGET, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_BARGAPTARGET

#ifdef XPRS_BARHGEXTRAPOLATE
    MPD( AddSolverOption_MergeDuplicates("pdhg:xprs_barhgextrapolate XPRS_BARHGEXTRAPOLATE",
      "                         Extrapolation parameter for the hybrid gradient algorithm. Although theory suggests that a value of 1 is best, slightly smaller values perform better in general.                 "
      "\n\nDefault:                          0.15                 ",
      XPRS_BARHGEXTRAPOLATE, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_BARHGEXTRAPOLATE

#ifdef XPRS_BARHGGPU
    MPD( AddSolverOption_MergeDuplicates("pdhg:xprs_barhggpu XPRS_BARHGGPU",
      "             Whether to use a GPU for the hybrid gradient algorithm. Even though the GPU implementation of the hybrid 			gradient algorithm is identical in operation and functionality to the CPU implementation, the returned solutions 			can differ between the two versions due to the different architecture of the GPU. 			 			GPU support is not available in the deterministic concurrent LP algorithm. 		"
      "\n\n"
      "Values (default:  			0, do not use a GPU. 		):\n"
      "\n- (0)  Do not use a GPU."
      "\n- (1)  Use the GPU.",
      XPRS_BARHGGPU, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BARHGGPU

#ifdef XPRS_BARHGGPUBLOCKSIZE
    MPD( AddSolverOption_MergeDuplicates("pdhg:xprs_barhggpublocksize XPRS_BARHGGPUBLOCKSIZE",
      "             The size of CUDA blocks to use for the GPU calculations. 		"
      "\n\nDefault:  			256 		",
      XPRS_BARHGGPUBLOCKSIZE, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BARHGGPUBLOCKSIZE

#ifdef XPRS_BARHGMAXRESTARTS
    MPD( AddSolverOption_MergeDuplicates("pdhg:xprs_barhgmaxrestarts XPRS_BARHGMAXRESTARTS",
      "                 The maximum number of restarts in the hybrid gradient algorithm. Restarts play the role of iterations in the hybrid gradient algorithm.                 A log line is printed at every restart, unless BAROUTPUT is set to 0.         "
      "\n\nDefault:                  1250         ",
      XPRS_BARHGMAXRESTARTS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BARHGMAXRESTARTS

#ifdef XPRS_BARHGOPS
    MPD( AddSolverOption_MergeDuplicates("pdhg:xprs_barhgops XPRS_BARHGOPS",
      "                         Bit-vector control (see Section Bit-vector controls) options for the hybrid gradient algorithm. Bits 1, 2 and 3 control which norms of the coefficient matrix are used for solution normalization. The normalization factor is the maximum of the selected norms. By default, or if all three bits are set to 0, the infinity norm is used.                         The omega parameter referenced in bits 4, 5 and 6 is a measure of the relative magnitudes of the objective and the right-hand side.                 "
      "\n\n"
      "Values (default:                          24, the infinity norm is used for initialization and the L2 norm for measuring the solution quality.                 ):\n"
      "\n- (0)  Use an asymmetric average for the primal averaging."
      "\n- (1)  Use the 1-norm of the coefficient matrix in normalizing the initial solution."
      "\n- (2)  Use the 2-norm of the coefficient matrix in normalizing the initial solution."
      "\n- (3)  Use the infinity norm of the coefficient matrix in normalizing the initial solution."
      "\n- (4)  Use L2 norm to measure solution quality. "
      "\n- (5)  Contract omega towards 1 if the infeasibility is small enough."
      "\n- (6)  Omega is based on the infeasibility.",
      XPRS_BARHGOPS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BARHGOPS

#ifdef XPRS_BARHGPRECISION
    MPD( AddSolverOption_MergeDuplicates("pdhg:xprs_barhgprecision XPRS_BARHGPRECISION",
      " 			Whether to use single or double precision floating-point arithmetic in the hybrid gardient algorithm. The single precision 			implementation uses less memory and is, in general, faster than the double precision implementation.  			 			This control applies to both the CPU and the GPU implementation of the algorithm. The performance difference is greater 			for the GPU version. 		"
      "\n\n"
      "Values (default:  			-1, use double precision for CPU platforms and single precision on GPU platforms. 		):\n"
      "\n- (-1)  Automatically selected based on the value of BARHGGPU: 			single precision arithmetic is used if BARHGGPU is 1 (GPU execution), 			and double precision arithmetic is used	if BARHGGPU is 0 (CPU execution)."
      "\n- (0)  Use single precision arithmetic on both CPU and GPU platforms."
      "\n- (1)  Use double precision arithmetic on both CPU and GPU platforms.",
      XPRS_BARHGPRECISION, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BARHGPRECISION

#ifdef XPRS_BARHGRELTOL
    MPD( AddSolverOption_MergeDuplicates("pdhg:xprs_barhgreltol XPRS_BARHGRELTOL",
      "                         Relative feasibility tolerance for the hybrid gradient algorithm.                 "
      "\n\nDefault:                          0, determined automatically.                 ",
      XPRS_BARHGRELTOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_BARHGRELTOL

#ifdef XPRS_BARINDEFLIMIT
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_barindeflimit XPRS_BARINDEFLIMIT",
      "Newton Barrier. This limits the number of consecutive indefinite barrier iterations that will be performed. The optimizer will try to minimize (resp. maximize) a QP problem even if the Q matrix is not positive (resp. negative) semi-definite. However, the optimizer may detect that the Q matrix is indefinite and this can result in the optimizer not converging. This control specifies how many indefinite iterations may occur before the optimizer stops and reports that the problem is indefinite. It is usual to specify a value greater than one, and only stop after a series of indefinite matrices, as the problem may be found to be indefinite incorrectly on a few iterations for numerical reasons."
      "\n\nDefault: 15",
      XPRS_BARINDEFLIMIT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BARINDEFLIMIT

#ifdef XPRS_BARITERATIVE
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_bariterative XPRS_BARITERATIVE",
      " The maximum number of barrier iterations in which an iterative solver is used instead of the Cholesky decomposition. "
      "\n\n"
      "Values (default: -2):\n"
      "\n- (-2)  Automatically determined."
      "\n- (-1)  Turn iterative solver off."
      "\n- (0)  Use iterative solver for the starting point computation."
      "\n- (n>0)  Try to apply iterative solver for the first n barrier iterations.",
      XPRS_BARITERATIVE, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BARITERATIVE

#ifdef XPRS_BARITERLIMIT
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_bariterlimit XPRS_BARITERLIMIT",
      "Newton barrier: The maximum number of iterations. While the simplex method usually performs a number of iterations which is proportional to the number of constraints (rows) in a problem, the barrier method standardly finds the optimal solution to a given accuracy after a number of iterations which is independent of the problem size. The penalty is rather that the time for each iteration increases with the size of the problem. BARITERLIMIT specifies the maximum number of iterations which will be carried out by the barrier."
      "\n\nDefault: 500",
      XPRS_BARITERLIMIT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BARITERLIMIT

#ifdef XPRS_BARKERNEL
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_barkernel XPRS_BARKERNEL",
      "Newton barrier: Defines how centrality is weighted in the barrier algorithm."
      "\n\n"
      "Values (default: 0.0):\n"
      "\n- (>=+1.0)  Increases the emphasis on centrality when larger value is set."
      "\n- (<=-1.0)  Selects a value adaptively in every iteration from [+1, -BARKERNEL].",
      XPRS_BARKERNEL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_BARKERNEL

#ifdef XPRS_BARLARGEBOUND
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_barlargebound XPRS_BARLARGEBOUND",
      "Threshold for the barrier to handle large bounds."
      "\n\nDefault: 0",
      XPRS_BARLARGEBOUND, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_BARLARGEBOUND

#ifdef XPRS_BAROBJPERTURB
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_barobjperturb XPRS_BAROBJPERTURB",
      "Defines how the barrier perturbs the objective."
      "\n\n"
      "Values (default: 1e-6):\n"
      "\n- (>0)  Let the optimizer decide if the objective is perturbed or not and use the parameter value as the scale of the perturbation."
      "\n- (0)  Turn off objective perturbation."
      "\n- (<0)  Always perturb the objective by the absolute value of the parameter.",
      XPRS_BAROBJPERTURB, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_BAROBJPERTURB

#ifdef XPRS_BAROBJSCALE
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_barobjscale XPRS_BAROBJSCALE",
      "Defines how the barrier scales the objective."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Let the optimizer decide."
      "\n- (0)  Scale by geometric mean."
      "\n- (>=0)  Scale such that the largest objective coefficient's largest element does not exceed this number. In quadratic problems, the quadratic diagonal is used as reference valuses instead of the linear objective.",
      XPRS_BAROBJSCALE, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_BAROBJSCALE

#ifdef XPRS_BARORDER
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_barorder XPRS_BARORDER",
      "Newton barrier: This controls the Cholesky factorization in the Newton-Barrier. "
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  Choose automatically."
      "\n- (1)  Minimum degree method. This selects diagonal elements with the smallest number of nonzeros in their rows or columns."
      "\n- (2)  Minimum local fill method. This considers the adjacency graph of nonzeros in the matrix and seeks to eliminate nodes that minimize the creation of new edges."
      "\n- (3)  Nested dissection method. This considers the adjacency graph and recursively seeks to separate it into non-adjacent pieces.",
      XPRS_BARORDER, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BARORDER

#ifdef XPRS_BARORDERTHREADS
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_barorderthreads XPRS_BARORDERTHREADS",
      "If set to a positive integer it determines the number of concurrent threads for the sparse matrix ordering algorithm in the Newton-barrier method. "
      "\n\n"
      "Values (default: 0 ):\n"
      "\n- (0)  The default value is determined automatically based on the problem size, structure and algorithm choice."
      "\n- (>=0)  The number of concurrent threads for the sparse matrix ordering algorithm in the Newton-barrier method.",
      XPRS_BARORDERTHREADS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BARORDERTHREADS

#ifdef XPRS_BAROUTPUT
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_baroutput XPRS_BAROUTPUT",
      "Newton barrier and hybrid gradient: This specifies the level of solution output provided. Output is provided either after each iteration of the algorithm, or else can be turned off completely by this parameter. "
      "\n\n"
      "Values (default: 1):\n"
      "\n- (0)  No output."
      "\n- (1)  At each iteration.",
      XPRS_BAROUTPUT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BAROUTPUT

#ifdef XPRS_BARPERTURB
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_barperturb XPRS_BARPERTURB",
      "Newton barrier: In numerically challenging cases it is often advantageous to apply perturbations on the KKT system to improve its numerical properties. BARPERTURB controlls how much perturbation is allowed during the barrier iterations. By default no perturbation is allowed. Set this parameter with care as larger perturbations may lead to less efficient iterates and the best settings are problem-dependent. "
      "\n\nDefault: 0",
      XPRS_BARPERTURB, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_BARPERTURB

#ifdef XPRS_BARPRESOLVEOPS
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_barpresolveops XPRS_BARPRESOLVEOPS",
      "Newton barrier: This bit-vector (see Section Bit-vector controls) controls the Newton-Barrier specific presolve operations. "
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  Use standard presolve."
      "\n- (1)  Extra effort is spent in barrier specific presolve."
      "\n- (2)  Do full matrix eliminations (reduce matrix size).",
      XPRS_BARPRESOLVEOPS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BARPRESOLVEOPS

#ifdef XPRS_BARPRIMALSTOP
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_barprimalstop XPRS_BARPRIMALSTOP",
      "Newton barrier and hybrid gradient: This is a convergence parameter, indicating the tolerance for primal infeasibilities. If the difference between the constraints and their bounds in the primal problem falls below this tolerance in absolute value, the Optimizer will terminate and return the current solution."
      "\n\n"
      "Values (default: 0 ):\n"
      "\n- (0)  The default value is determined automatically based on the problem size, structure and algorithm choice."
      "\n- (>=0)  The tolerance for primal infeasibilities.",
      XPRS_BARPRIMALSTOP, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_BARPRIMALSTOP

#ifdef XPRS_BARREFITER
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_barrefiter XPRS_BARREFITER",
      "Newton barrier: After terminating the barrier algorithm, further refinement steps can be performed. Such refinement steps are especially helpful if the solution is near to the optimum and can improve primal feasibility and decrease the complementarity gap. It is also often advantageous for the crossover algorithm. BARREFITER specifies the maximum number of such refinement iterations. "
      "\n\nDefault: 0",
      XPRS_BARREFITER, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BARREFITER

#ifdef XPRS_BARREGULARIZE
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_barregularize XPRS_BARREGULARIZE",
      "This bit-vector control (see Section Bit-vector controls) determines how the barrier algorithm applies regularization on the KKT system."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (0)  Standard regularization is turned on/off."
      "\n- (1)  Reduced regularization is turned on/off. This option reduces the perturbation effect of the standard regularization."
      "\n- (2)  Forces to keep dependent rows in the KKT system."
      "\n- (3)  Forces to preserve degenerate rows in the KKT system."
      "\n- (4)  Restrict regularization to infeasible iterates."
      "\n- (5)  Disable iterative regularizations."
      "\n- (6)  Apply iterative regularization more often.",
      XPRS_BARREGULARIZE, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BARREGULARIZE

#ifdef XPRS_BARRHSSCALE
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_barrhsscale XPRS_BARRHSSCALE",
      "Defines how the barrier scales the right hand side."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Let the optimizer decide."
      "\n- (0)  Scale by geometric mean."
      "\n- (>=0)  Scale such that the largest right hand side coefficient's largest element does not exceed this number.",
      XPRS_BARRHSSCALE, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_BARRHSSCALE

#ifdef XPRS_BARSOLUTION
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_barsolution XPRS_BARSOLUTION",
      "This determines whether the barrier has to decide which is the best solution found or return the solution computed by the last barrier iteration. "
      "\n\n"
      "Values (default: 0):\n"
      "\n- (-1)  (callback only: do not save current soulution as the best one)."
      "\n- (0)  return the best solution found (in callback: let the barrier decide the current solution is the best or not)."
      "\n- (1)  return the last barrier iteration (in callback: save current solution as the best solution so far).",
      XPRS_BARSOLUTION, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BARSOLUTION

#ifdef XPRS_BARSTART
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_barstart XPRS_BARSTART",
      "Controls the computation of the starting point and warm-starting for the Newton barrier and the hybrid gradient algorithms."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (-1)  Uses the existing solution for warm-start if one is available."
      "\n- (0)  Warm-start is disabled; the starting point is determined automatically from the next three options."
      "\n- (1)  Uses simple heuristics to compute the starting point based on the magnitudes of the matrix entries."
      "\n- (2)  Uses the pseudoinverse of the constraint matrix to determine primal and dual initial solutions.                     Less sensitive to scaling and numerically more robust, but in several case less efficient than 1."
      "\n- (3)  Uses the unit starting point for the homogeneous self-dual barrier algorithm.",
      XPRS_BARSTART, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BARSTART

#ifdef XPRS_BARSTARTWEIGHT
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_barstartweight XPRS_BARSTARTWEIGHT",
      "Newton barrier: This sets a weight for the warm-start point when warm-start is set for the barrier algorithm. Using larger weight gives more emphasis for the supplied starting point. "
      "\n\nDefault: 0.85",
      XPRS_BARSTARTWEIGHT, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_BARSTARTWEIGHT

#ifdef XPRS_BARSTEPSTOP
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_barstepstop XPRS_BARSTEPSTOP",
      "Newton barrier: A convergence parameter, representing the minimal step size. On each iteration of the barrier algorithm, a step is taken along a computed search direction. If that step size is smaller than BARSTEPSTOP, the Optimizer will terminate and return the current solution."
      "\n\nDefault: 1.0E-16",
      XPRS_BARSTEPSTOP, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_BARSTEPSTOP

#ifdef XPRS_BARTHREADS
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_barthreads XPRS_BARTHREADS",
      "         If set to a positive integer it determines the number of threads implemented to run the Newton-barrier and hybrid gradient algorithms.         If the value is set to the default value (-1), the THREADS control will determine the number of threads used."
      "\n\nDefault: -1(determined by the THREADS control)",
      XPRS_BARTHREADS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BARTHREADS

#ifdef XPRS_BIGM
    MPD( AddSolverOption_MergeDuplicates("lp:xprs_bigm XPRS_BIGM",
      "The infeasibility penalty used if the 'Big M' method is implemented. "
      "\n\nDefault: Dependent on the matrix characteristics.",
      XPRS_BIGM, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_BIGM

#ifdef XPRS_BIGMMETHOD
    MPD( AddSolverOption_MergeDuplicates("lp:xprs_bigmmethod XPRS_BIGMMETHOD",
      "Simplex: This specifies whether to use the 'Big M' method, or the standard phase I (achieving feasibility) and phase II (achieving optimality). In the 'Big M' method, the objective coefficients of the variables are considered during the feasibility phase, possibly leading to an initial feasible basis which is closer to optimal. The side-effects involve possible round-off errors due to the presence of the 'Big M' factor in the problem. "
      "\n\n"
      "Values (default: 1):\n"
      "\n- (0)  For phase I / phase II."
      "\n- (1)  If 'Big M' method to be used.",
      XPRS_BIGMMETHOD, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BIGMMETHOD

#ifdef XPRS_BRANCHCHOICE
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_branchchoice XPRS_BRANCHCHOICE",
      "Once a MIP entity has been selected for branching, this control determines which of the branches is solved first."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  Minimum estimate branch first."
      "\n- (1)  Maximum estimate branch first."
      "\n- (2)  If an incumbent solution exists, solve the branch satisfied by that solution first. Otherwise solve the minimum estimate branch first (option 0)."
      "\n- (3)  Solve first the branch that forces the value of the branching variable to move farther away from the value it had at the root node. If the branching entity is not a simple variable, solve the minimum estimate branch first (option 0).",
      XPRS_BRANCHCHOICE, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BRANCHCHOICE

#ifdef XPRS_BRANCHDISJ
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_branchdisj XPRS_BRANCHDISJ",
      "Branch and Bound: Determines whether the optimizer should attempt to branch on general split disjunctions during the branch and bound search. "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatic selection of the strategy."
      "\n- (0)  Disabled."
      "\n- (1)  Cautious strategy. Disjunctive branches will be created only for general integers with a wide range."
      "\n- (2)  Moderate strategy."
      "\n- (3)  Aggressive strategy. Disjunctive branches will be created for both binaries and integers.",
      XPRS_BRANCHDISJ, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BRANCHDISJ

#ifdef XPRS_BRANCHSTRUCTURAL
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_branchstructural XPRS_BRANCHSTRUCTURAL",
      "Branch and Bound: Determines whether the optimizer should search for special structure in the problem to branch on during the branch and bound search. "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatically determined."
      "\n- (0)  Disabled."
      "\n- (1)  Enabled.",
      XPRS_BRANCHSTRUCTURAL, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BRANCHSTRUCTURAL

#ifdef XPRS_BREADTHFIRST
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_breadthfirst XPRS_BREADTHFIRST",
      "The number of nodes to include in the best-first search before switching to the local first search (NODESELECTION=4)."
      "\n\nDefault: 11",
      XPRS_BREADTHFIRST, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_BREADTHFIRST

#ifdef XPRS_CACHESIZE
    MPD( AddSolverOption_MergeDuplicates("sys:xprs_cachesize XPRS_CACHESIZE",
      "This parameter is deprecated and will be removed in a future release. Newton Barrier: L2 or L3 (see notes) cache size in kB (kilobytes) of the CPU. On Intel (or compatible) platforms a value of -1 may be used to determine the cache size automatically. If the CPU model is new then the cache size may not be correctly detected by an older release of the software. "
      "\n\nDefault: -1",
      XPRS_CACHESIZE, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_CACHESIZE

#ifdef XPRS_CALLBACKCHECKTIMEDELAY
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_callbackchecktimedelay XPRS_CALLBACKCHECKTIMEDELAY",
      "Minimum delay in milliseconds between two consecutive executions of the CHECKTIME callback in the same solution process"
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  Callback delay is disabled - the callback is executed every time;"
      "\n- (n>0)  Callback invocation is suppressed if less than n milliseconds have passed since the last invocation.",
      XPRS_CALLBACKCHECKTIMEDELAY, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_CALLBACKCHECKTIMEDELAY

#ifdef XPRS_CALLBACKCHECKTIMEWORKDELAY
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_callbackchecktimeworkdelay XPRS_CALLBACKCHECKTIMEWORKDELAY",
      "Minimum delay in work units between two consecutive executions of the CHECKTIME callback in the same solution process"
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  Callback delay is disabled - the callback is executed every time;"
      "\n- (n>0)  Callback invocation if less than n work units, which may be a fraction of a work unit, have passed since the last invocation.",
      XPRS_CALLBACKCHECKTIMEWORKDELAY, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_CALLBACKCHECKTIMEWORKDELAY

#ifdef XPRS_CALLBACKFROMMAINTHREAD
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_callbackfrommainthread XPRS_CALLBACKFROMMAINTHREAD",
      "Branch and Bound: specifies whether the MIP callbacks should only be called on the main thread."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  Invoke callbacks on worker threads during parallel MIP;"
      "\n- (1)  Only ever invoke a callback on the thread that called XPRSmipoptimize.",
      XPRS_CALLBACKFROMMAINTHREAD, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_CALLBACKFROMMAINTHREAD

#ifdef XPRS_CHECKINPUTDATA
    MPD( AddSolverOption_MergeDuplicates("prob:xprs_checkinputdata XPRS_CHECKINPUTDATA",
      "Check input arrays for bad data."
      "\n\n"
      "Values (default: 1 (on)):\n"
      "\n- (0)  Do not check."
      "\n- (1)  Check input arrays.",
      XPRS_CHECKINPUTDATA, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_CHECKINPUTDATA

#ifdef XPRS_CHOLESKYALG
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_choleskyalg XPRS_CHOLESKYALG",
      "Newton barrier: type of Cholesky factorization used; bit-vector-control (see Section Bit-vector controls)."
      "\n\n"
      "Values (default: -1 (automatic)):\n"
      "\n- (0)  matrix blocking: 0: automatic setting; 1: manual setting."
      "\n- (1)  if manual selection of matrix blocking: 0: multi-pass; 1: single-pass."
      "\n- (2)  nonseparable QP relaxation: 0: off; 1: on."
      "\n- (3)  corrector weight: 0: automatic setting; 1: manual setting."
      "\n- (4)  if manual selection of corrector weight: 0: off; 1: on."
      "\n- (5)  refinement: 0: automatic setting; 1: manual setting."
      "\n- (6)  preconditioned conjugate gradient method (PCGM): 0: PCGM off; 1: PCGM on."
      "\n- (7)  Preconditioned quasi minimal residual (QMR) to refine solution: 0: QMR off; 1: QMR on."
      "\n- (8)  Perform refinement on the augmented system 0: off; 1: on."
      "\n- (9)  Force highest accuracy in refinement 0: off; 1: on.",
      XPRS_CHOLESKYALG, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_CHOLESKYALG

#ifdef XPRS_CHOLESKYTOL
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_choleskytol XPRS_CHOLESKYTOL",
      "Newton barrier: The tolerance for pivot elements in the Cholesky decomposition of the normal equations coefficient matrix, computed at each iteration of the barrier algorithm. If the absolute value of the pivot element is less than or equal to CHOLESKYTOL, it merits special treatment in the Cholesky decomposition process."
      "\n\nDefault: 1.0E-15",
      XPRS_CHOLESKYTOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_CHOLESKYTOL

#ifdef XPRS_CLAMPING
    MPD( AddSolverOption_MergeDuplicates("sol:xprs_clamping XPRS_CLAMPING",
      "This bit-vector control (see Section Bit-vector controls) allows for the adjustment of returned solution values such that they are always within bounds."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (-1)  Determined automatically."
      "\n- (0)  Adjust primal solution to always be within primal bounds. Slacks if provided will be adjusted accordingly."
      "\n- (1)  Adjust primal slack values to always be within constraint bounds."
      "\n- (2)  Adjust dual solution to always be within the dual bounds implied by the slacks. Reduced costs, if provided, will be adjusted accordingly."
      "\n- (3)  Adjust reduced costs to always be within dual bounds implied by the primal solution.",
      XPRS_CLAMPING, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_CLAMPING

#ifdef XPRS_COMPUTE
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_compute XPRS_COMPUTE",
      "Controls whether the next solve is performed directly or on an Insight Compute Interface. "
      "\n\n"
      "Values (default: Depends on environment):\n"
      "\n- (0)  Solve locally."
      "\n- (1)  Solve using an Insight Compute Interface.",
      XPRS_COMPUTE, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_COMPUTE

#ifdef XPRS_COMPUTEEXECSERVICE
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_computeexecservice XPRS_COMPUTEEXECSERVICE",
      "Selects the Insight execution service that will be used for solving remote optimizations."
      "\n\nDefault: Empty string",
      XPRS_COMPUTEEXECSERVICE) );
#endif  // ifdef XPRS_COMPUTEEXECSERVICE

#ifdef XPRS_COMPUTEJOBPRIORITY
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_computejobpriority XPRS_COMPUTEJOBPRIORITY",
      "Selects the priority that will be used for remote optimization jobs."
      "\n\nDefault: 0",
      XPRS_COMPUTEJOBPRIORITY, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_COMPUTEJOBPRIORITY

#ifdef XPRS_COMPUTELOG
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_computelog XPRS_COMPUTELOG",
      "Controls how the run log is fetched when a solve is performed on an Insight Compute Interface. "
      "\n\n"
      "Values (default: 1):\n"
      "\n- (0)  Run log will not be fetched"
      "\n- (1)  Run log will be fetched in real-time"
      "\n- (2)  Run log will be fetched at the end of the solve"
      "\n- (3)  Run log will be fetched at the end of the solve if the solve fails with an error",
      XPRS_COMPUTELOG, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_COMPUTELOG

#ifdef XPRS_CONCURRENTTHREADS
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_concurrentthreads XPRS_CONCURRENTTHREADS",
      " Determines the number of threads used by the concurrent solver. "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Determined automatically"
      "\n- (>0)  Number of threads to use.",
      XPRS_CONCURRENTTHREADS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_CONCURRENTTHREADS

#ifdef XPRS_CONFLICTCUTS
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_conflictcuts XPRS_CONFLICTCUTS",
      "Branch and Bound: Specifies how cautious or aggressive the optimizer should be when searching for and applying conflict cuts. Conflict cuts are in-tree cuts derived from nodes found to be infeasible or cut off, which can be used to cut off other branches of the search tree."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatic."
      "\n- (0)  Disable conflict cuts."
      "\n- (1)  Cautious application of conflict cuts."
      "\n- (2)  Medium application of conflict cuts."
      "\n- (3)  Aggressive application of conflict cuts.",
      XPRS_CONFLICTCUTS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_CONFLICTCUTS

#ifdef XPRS_CORESPERCPU
    MPD( AddSolverOption_MergeDuplicates("sys:xprs_corespercpu XPRS_CORESPERCPU",
      "Used to override the detected value of the number of cores on a CPU. The cache size (either detected or specified via the CACHESIZE control) used in Barrier methods will be divided by this amount, and this scaled-down value will be the amount of cache allocated to each Barrier thread"
      "\n\nDefault: -1",
      XPRS_CORESPERCPU, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_CORESPERCPU

#ifdef XPRS_COVERCUTS
    MPD( AddSolverOption_MergeDuplicates("cut:xprs_covercuts XPRS_COVERCUTS",
      "Branch and Bound: The number of rounds of lifted cover inequalities at the root node. A lifted cover inequality is an additional constraint that can be particularly effective at reducing the size of the feasible region without removing potential integral solutions. The process of generating these can be carried out a number of times, further reducing the feasible region, albeit incurring a time penalty. There is usually a good payoff from generating these at the root node, since these inequalities then apply to every subsequent node in the tree search."
      "\n\nDefault: -1 — determined automatically.",
      XPRS_COVERCUTS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_COVERCUTS

#ifdef XPRS_CPIALPHA
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_cpialpha XPRS_CPIALPHA",
      "decay term for confined primal integral computation."
      "\n\nDefault: 0",
      XPRS_CPIALPHA, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_CPIALPHA

#ifdef XPRS_CPUPLATFORM
    MPD( AddSolverOption_MergeDuplicates("sys:xprs_cpuplatform XPRS_CPUPLATFORM",
      "Newton Barrier: Selects the AMD, Intel x86 or ARM vectorization instruction set that Barrier should run optimized code for. On AMD / Intel x86 platforms the SSE2, AVX and AVX2 instruction sets are supported while on ARM platforms the NEON architecture extension can be activated. "
      "\n\n"
      "Values (default: -2, using AVX2 instructions if supported by the CPU):\n"
      "\n- (-2)  Highest supported [Generic, SSE2, AVX or AVX2]."
      "\n- (-1)  Highest supported solve path consistent code [Generic, SSE2 or AVX]."
      "\n- (0)  Use generic code compatible with all CPUs."
      "\n- (1)  Use SSE2 / NEON optimized code."
      "\n- (2)  Use AVX optimized code."
      "\n- (3)  Use AVX2 optimized code.",
      XPRS_CPUPLATFORM, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_CPUPLATFORM

#ifdef XPRS_CPUTIME
    MPD( AddSolverOption_MergeDuplicates("sys:xprs_cputime XPRS_CPUTIME",
      "How time should be measured when timings are reported in the log and when checking against time limits"
      "\n\n"
      "Values (default: 0):\n"
      "\n- (-1)  Disable the timer."
      "\n- (0)  Use elapsed time."
      "\n- (1)  Use process time.",
      XPRS_CPUTIME, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_CPUTIME

#ifdef XPRS_CRASH
    MPD( AddSolverOption_MergeDuplicates("lp:xprs_crash XPRS_CRASH",
      "Simplex: This determines the type of crash used when the algorithm begins. During the crash procedure, an initial basis is determined which is as close to feasibility and triangularity as possible. A good choice at this stage will significantly reduce the number of iterations required to find an optimal solution. The possible values increase proportionally to their time-consumption."
      "\n\n"
      "Values (default: 2):\n"
      "\n- (0)  Turns off all crash procedures."
      "\n- (1)  For singletons only (one pass)."
      "\n- (2)  For singletons only (multi pass)."
      "\n- (3)  Multiple passes through the matrix considering slacks."
      "\n- (4)  Multiple (≤10) passes through the matrix but only doing slacks at the very end."
      "\n- (n>10)  As for value 4 but performing at most n - 10 passes."
      "\n- (0)  Perform standard crash."
      "\n- (1)  Perform additional numerical checks during crash."
      "\n- (2)  Extend the set of column candidates for crash. "
      "\n- (3)  Extend the set of row candidates for crash. "
      "\n- (4)  Force crash, i.e., consider all suitable columns/rows as candidates for crash.",
      XPRS_CRASH, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_CRASH

#ifdef XPRS_CROSSOVER
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_crossover XPRS_CROSSOVER",
      "Newton barrier and hybrid gradient: This control determines whether the barrier method will cross over to the simplex method when at optimal solution has been found, to provide an end basis (see XPRSgetbasis, XPRSwritebasis) and advanced sensitivity analysis information (see XPRSobjsa, XPRSrhssa, XPRSbndsa)."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Determined automatically."
      "\n- (0)  No crossover."
      "\n- (1)  Primal crossover first."
      "\n- (2)  Dual crossover first.",
      XPRS_CROSSOVER, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_CROSSOVER

#ifdef XPRS_CROSSOVERACCURACYTOL
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_crossoveraccuracytol XPRS_CROSSOVERACCURACYTOL",
      "Newton barrier: This control determines how crossover adjusts the default relative pivot tolerance. When re-inversion is necessary, crossover will compare the recalculated working basic solution with the assumed ones just before re-inversion took place. If the error is above this threshold, crossover will adjust the relative pivot tolerance to address the build-up of numerical inaccuracies."
      "\n\nDefault: 1e-6",
      XPRS_CROSSOVERACCURACYTOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_CROSSOVERACCURACYTOL

#ifdef XPRS_CROSSOVERITERLIMIT
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_crossoveriterlimit XPRS_CROSSOVERITERLIMIT",
      "Newton barrier and hybrid gradient: The maximum number of iterations that will be performed in the crossover procedure before the optimization process terminates. "
      "\n\nDefault: 2147483647",
      XPRS_CROSSOVERITERLIMIT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_CROSSOVERITERLIMIT

#ifdef XPRS_CROSSOVEROPS
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_crossoverops XPRS_CROSSOVEROPS",
      "Newton barrier and hybrid gradient: a bit-vector (see Section Bit-vector controls) for adjusting the behavior of the crossover procedure."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  Returned solution when the crossover terminates prematurely: 0: Return the last basis from the crossover; 1: Return the barrier solution."
      "\n- (1)  Select the crossover stages to be performed: 0: Perform both crossover stages; 1: Skip second crossover stage."
      "\n- (2)  Set crossover behaviour: 0: Force to perform all pivots; 1: Skip pivots that are numerically less reliable."
      "\n- (3)  Set crossover behaviour: 0: Perform standard crossover; 1: Perform a slower, but numerically more careful crossover.",
      XPRS_CROSSOVEROPS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_CROSSOVEROPS

#ifdef XPRS_CROSSOVERTHREADS
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_crossoverthreads XPRS_CROSSOVERTHREADS",
      "Determines the maximum number of threads that parallel crossover is allowed to use. If CROSSOVERTHREADS is set to the default value (-1), the BARTHREADS control will determine the number of threads used."
      "\n\nDefault: -1 (determined by the BARTHREADS control)",
      XPRS_CROSSOVERTHREADS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_CROSSOVERTHREADS

#ifdef XPRS_CUTDEPTH
    MPD( AddSolverOption_MergeDuplicates("cut:xprs_cutdepth XPRS_CUTDEPTH",
      "Branch and Bound: Sets the maximum depth in the tree search at which cuts will be generated. Generating cuts can take a lot of time, and is often less important at deeper levels of the tree since tighter bounds on the variables have already reduced the feasible region. A value of 0 signifies that no cuts will be generated."
      "\n\nDefault: -1 — determined automatically.",
      XPRS_CUTDEPTH, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_CUTDEPTH

#ifdef XPRS_CUTFACTOR
    MPD( AddSolverOption_MergeDuplicates("cut:xprs_cutfactor XPRS_CUTFACTOR",
      "Limit on the number of cuts and cut coefficients the optimizer is allowed to add to the matrix during tree search. The cuts and cut coefficients are limited by CUTFACTOR times the number of rows and coefficients in the initial matrix."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Let the optimizer decide on the maximum amount of cuts based on CUTSTRATEGY."
      "\n- (>=0)  Multiple of number of rows and coefficients to use.",
      XPRS_CUTFACTOR, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_CUTFACTOR

#ifdef XPRS_CUTFREQ
    MPD( AddSolverOption_MergeDuplicates("cut:xprs_cutfreq XPRS_CUTFREQ",
      "Branch and Bound: This specifies the frequency at which cuts are generated in the tree search. If the depth of the node modulo CUTFREQ is zero, then cuts will be generated. "
      "\n\nDefault: -1 — determined automatically.",
      XPRS_CUTFREQ, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_CUTFREQ

#ifdef XPRS_CUTSELECT
    MPD( AddSolverOption_MergeDuplicates("cut:xprs_cutselect XPRS_CUTSELECT",
      "A bit-vector (see Section Bit-vector controls) providing detailed control of the cuts created for the root node of a MIP solve. Use TREECUTSELECT to control cuts during the tree search."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (5)  Clique cuts."
      "\n- (6)  Mixed Integer Rounding (MIR) cuts."
      "\n- (7)  Lifted cover cuts."
      "\n- (8)  Turn on row aggregation for MIR cuts."
      "\n- (11)  Flow path cuts."
      "\n- (12)  Implication cuts."
      "\n- (13)  Turn on automatic Lift-and-Project cutting strategy."
      "\n- (14)  Disable cutting from cut rows."
      "\n- (15)  Lifted GUB cover cuts."
      "\n- (16)  Zero-half cuts."
      "\n- (17)  Indicator constraint cuts."
      "\n- (18)  Strong Chvatal-Gomory cuts."
      "\n- (20)  Farkas cuts.",
      XPRS_CUTSELECT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_CUTSELECT

#ifdef XPRS_CUTSTRATEGY
    MPD( AddSolverOption_MergeDuplicates("cut:xprs_cutstrategy XPRS_CUTSTRATEGY",
      "Branch and Bound: This specifies the cut strategy. A more aggressive cut strategy, generating a greater number of cuts, will result in fewer nodes to be explored, but with an associated time cost in generating the cuts. The fewer cuts generated, the less time taken, but the greater subsequent number of nodes to be explored."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatic selection of the cut strategy."
      "\n- (0)  No cuts."
      "\n- (1)  Conservative cut strategy."
      "\n- (2)  Moderate cut strategy."
      "\n- (3)  Aggressive cut strategy.",
      XPRS_CUTSTRATEGY, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_CUTSTRATEGY

#ifdef XPRS_DEFAULTALG
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_defaultalg XPRS_DEFAULTALG",
      "This selects the algorithm that will be used to solve LPs, standalone or during MIP optimization."
      "\n\n"
      "Values (default: 1):\n"
      "\n- (1)  Automatically determined."
      "\n- (2)  Dual simplex."
      "\n- (3)  Primal simplex."
      "\n- (4)  Newton barrier (or hybrid gradient, if BARALG=4 is set).",
      XPRS_DEFAULTALG, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_DEFAULTALG

#ifdef XPRS_DENSECOLLIMIT
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_densecollimit XPRS_DENSECOLLIMIT",
      "Newton barrier: Columns with more than DENSECOLLIMIT elements are considered to be dense. Such columns will be handled specially in the Cholesky factorization of this matrix."
      "\n\nDefault: 0 — determined automatically.",
      XPRS_DENSECOLLIMIT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_DENSECOLLIMIT

#ifdef XPRS_DETERMINISTIC
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_deterministic XPRS_DETERMINISTIC",
      "Selects whether to use a deterministic or opportunistic mode when solving a problem using multiple threads."
      "\n\n"
      "Values (default: 1):\n"
      "\n- (0)  Use opportunistic mode."
      "\n- (1)  Use deterministic mode."
      "\n- (2)  Use deterministic mode, except allow the initial concurrent continuous solve of a MIP to be opportunistic.",
      XPRS_DETERMINISTIC, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_DETERMINISTIC

#ifdef XPRS_DETERMINISTICLOG
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_deterministiclog XPRS_DETERMINISTICLOG",
      "Suppress non-deterministic log information in the standard MIP log. In particular, wall clock time stamps are replaced by (deterministic) work units."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  Report wall clock time stamps and other non-deterministic log information (the default)"
      "\n- (1)  Suppress non-deterministic log information. In particular, report deterministic work units instead of time.",
      XPRS_DETERMINISTICLOG, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_DETERMINISTICLOG

#ifdef XPRS_DUALGRADIENT
    MPD( AddSolverOption_MergeDuplicates("lp:xprs_dualgradient XPRS_DUALGRADIENT",
      "Simplex: This specifies the dual simplex pricing method."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Determined automatically."
      "\n- (0)  Devex."
      "\n- (1)  Steepest edge."
      "\n- (2)  Direct steepest edge."
      "\n- (3)  Sparse Devex.",
      XPRS_DUALGRADIENT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_DUALGRADIENT

#ifdef XPRS_DUALIZE
    MPD( AddSolverOption_MergeDuplicates("lp:xprs_dualize XPRS_DUALIZE",
      "For a linear problem or the initial linear relaxation of a MIP, determines whether to form and solve the dual problem."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Determine automatically which version would be faster."
      "\n- (0)  Solve the original problem."
      "\n- (1)  Solve the dualized problem.",
      XPRS_DUALIZE, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_DUALIZE

#ifdef XPRS_DUALIZEOPS
    MPD( AddSolverOption_MergeDuplicates("lp:xprs_dualizeops XPRS_DUALIZEOPS",
      "Bit-vector control (see Section Bit-vector controls) for adjusting the behavior when a problem is dualized."
      "\n\n"
      "Values (default: 1 (bit 0 is set)):\n"
      "\n- (0)  Swap the simplex algorithm to run. If dual simplex is selected for the original problem then primal simplex will be run on the dualized problem, and simiarly if primal simplex is selected.",
      XPRS_DUALIZEOPS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_DUALIZEOPS

#ifdef XPRS_DUALPERTURB
    MPD( AddSolverOption_MergeDuplicates("lp:xprs_dualperturb XPRS_DUALPERTURB",
      "The factor by which the problem will be perturbed prior to optimization by dual simplex. A value of 0.0 results in no perturbation prior to optimization.  Note the interconnection to the AUTOPERTURB control. If AUTOPERTURB is set to 1, the decision whether to perturb or not is left to the Optimizer. When the problem is automatically perturbed in dual simplex, however, the value of DUALPERTURB will be used for perturbation. "
      "\n\nDefault: -1 — determined automatically. ",
      XPRS_DUALPERTURB, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_DUALPERTURB

#ifdef XPRS_DUALSTRATEGY
    MPD( AddSolverOption_MergeDuplicates("lp:xprs_dualstrategy XPRS_DUALSTRATEGY",
      "This bit-vector control (see Section Bit-vector controls) specifies the dual simplex strategy."
      "\n\n"
      "Values (default: 1):\n"
      "\n- (0)  Switch to primal when re-optimization goes dual infeasible and numerically unstable."
      "\n- (1)  When dual intend to switch to primal, stop the solve instead of switching to primal."
      "\n- (2)  Use more aggressive cut-off in MIP search."
      "\n- (3)  Use dual simplex to remove cost perturbations."
      "\n- (4)  Enable more aggressive dual pivoting strategy."
      "\n- (5)  Keep using dual simplex even when it's numerically unstable.",
      XPRS_DUALSTRATEGY, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_DUALSTRATEGY

#ifdef XPRS_DUALTHREADS
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_dualthreads XPRS_DUALTHREADS",
      "Determines the maximum number of threads that dual simplex is allowed to use. If DUALTHREADS is set to the default value (-1), the THREADS control will determine the number of threads used. "
      "\n\nDefault: -1 (determined by the THREADS control)",
      XPRS_DUALTHREADS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_DUALTHREADS

#ifdef XPRS_EIGENVALUETOL
    MPD( AddSolverOption_MergeDuplicates("qp:xprs_eigenvaluetol XPRS_EIGENVALUETOL",
      "A quadratic matrix is considered not to be positive semi-definite, if its smallest eigenvalue is smaller than the negative of this value."
      "\n\nDefault: 1E-6",
      XPRS_EIGENVALUETOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_EIGENVALUETOL

#ifdef XPRS_ELIMFILLIN
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_elimfillin XPRS_ELIMFILLIN",
      "Amount of fill-in allowed when performing an elimination in presolve ."
      "\n\nDefault: 7",
      XPRS_ELIMFILLIN, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_ELIMFILLIN

#ifdef XPRS_ELIMTOL
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_elimtol XPRS_ELIMTOL",
      "The Markowitz tolerance for the elimination phase of the presolve."
      "\n\nDefault: 0.001",
      XPRS_ELIMTOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_ELIMTOL

#ifdef XPRS_ESCAPENAMES
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_escapenames XPRS_ESCAPENAMES",
      "If characters illegal to an mps or lp file should be escaped to guarantee readability, and whether escaped characters should be transformed back when reading such a file."
      "\n\n"
      "Values (default: 1):\n"
      "\n- (0)  Illegal characters are not escaped."
      "\n- (1)  Illegal characters are escaped.",
      XPRS_ESCAPENAMES, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_ESCAPENAMES

#ifdef XPRS_ETATOL
    MPD( AddSolverOption_MergeDuplicates("sim:xprs_etatol XPRS_ETATOL",
      "Tolerance on eta elements. During each iteration, the basis inverse is premultiplied by an elementary matrix, which is the identity except for one column - the eta vector. Elements of eta vectors whose absolute value is smaller than ETATOL are taken to be zero in this step."
      "\n\nDefault: 1.0E-13",
      XPRS_ETATOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_ETATOL

#ifdef XPRS_EXTRACOLS
    MPD( AddSolverOption_MergeDuplicates("prob:xprs_extracols XPRS_EXTRACOLS",
      "The initial number of extra columns to allow for in the matrix. If columns are to be added to the matrix, then, for maximum efficiency, space should be reserved for the columns before the matrix is input by setting the EXTRACOLS control. If this is not done, resizing will occur automatically, but more space may be allocated than the user actually requires."
      "\n\nDefault: 0",
      XPRS_EXTRACOLS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_EXTRACOLS

#ifdef XPRS_EXTRAELEMS
    MPD( AddSolverOption_MergeDuplicates("prob:xprs_extraelems XPRS_EXTRAELEMS",
      "The initial number of extra matrix elements to allow for in the matrix, including coefficients for cuts. If rows or columns are to be added to the matrix, then, for maximum efficiency, space should be reserved for the extra matrix elements before the matrix is input by setting the EXTRAELEMS control. If this is not done, resizing will occur automatically, but more space may be allocated than the user actually requires."
      "\n\nDefault: Hardware/platform dependent.",
      XPRS_EXTRAELEMS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_EXTRAELEMS

#ifdef XPRS_EXTRAMIPENTS
    MPD( AddSolverOption_MergeDuplicates("prob:xprs_extramipents XPRS_EXTRAMIPENTS",
      "The initial number of extra MIP entities to allow for."
      "\n\nDefault: 0",
      XPRS_EXTRAMIPENTS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_EXTRAMIPENTS

#ifdef XPRS_EXTRAROWS
    MPD( AddSolverOption_MergeDuplicates("prob:xprs_extrarows XPRS_EXTRAROWS",
      "The initial number of extra rows to allow for in the matrix, including cuts. If rows are to be added to the matrix, then, for maximum efficiency, space should be reserved for the rows before the matrix is input by setting the EXTRAROWS control. If this is not done, resizing will occur automatically, but more space may be allocated than the user actually requires."
      "\n\nDefault: 0",
      XPRS_EXTRAROWS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_EXTRAROWS

#ifdef XPRS_EXTRASETELEMS
    MPD( AddSolverOption_MergeDuplicates("prob:xprs_extrasetelems XPRS_EXTRASETELEMS",
      "The initial number of extra elements in sets to allow for in the matrix. If sets are to be added to the matrix, then, for maximum efficiency, space should be reserved for the set elements before the matrix is input by setting the EXTRASETELEMS control. If this is not done, resizing will occur automatically, but more space may be allocated than the user actually requires."
      "\n\nDefault: 0",
      XPRS_EXTRASETELEMS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_EXTRASETELEMS

#ifdef XPRS_EXTRASETS
    MPD( AddSolverOption_MergeDuplicates("prob:xprs_extrasets XPRS_EXTRASETS",
      "The initial number of extra sets to allow for in the matrix. If sets are to be added to the matrix, then, for maximum efficiency, space should be reserved for the sets before the matrix is input by setting the EXTRASETS control. If this is not done, resizing will occur automatically, but more space may be allocated than the user actually requires."
      "\n\nDefault: 0",
      XPRS_EXTRASETS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_EXTRASETS

#ifdef XPRS_FEASIBILITYJUMP
    MPD( AddSolverOption_MergeDuplicates("heur:xprs_feasibilityjump XPRS_FEASIBILITYJUMP",
      "MIP: Decides if the Feasibility Jump heuristic should be run. The value for this control is either -1 (let Xpress decide), 0 (off) or a value that indicates for which type of models the heuristic should be run. "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Use automatic settings."
      "\n- (0)  Turned off."
      "\n- (1)  Run the heuristic on models with all integer variables."
      "\n- (2)  Run the heuristic on models in which all non-integer variables have bounds [0,1]."
      "\n- (3)  Run the heuristic on models in which all non-integer variables have integer bounds.",
      XPRS_FEASIBILITYJUMP, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_FEASIBILITYJUMP

#ifdef XPRS_FEASIBILITYPUMP
    MPD( AddSolverOption_MergeDuplicates("heur:xprs_feasibilitypump XPRS_FEASIBILITYPUMP",
      "Branch and Bound: Decides if the Feasibility Pump heuristic should be run at the root node."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatic."
      "\n- (0)  Turned off."
      "\n- (1)  Always try the Feasibility Pump."
      "\n- (2)  Try the Feasibility Pump only if other heuristics have failed to find an integer solution.",
      XPRS_FEASIBILITYPUMP, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_FEASIBILITYPUMP

#ifdef XPRS_FEASTOL
    MPD( AddSolverOption_MergeDuplicates("tol:xprs_feastol XPRS_FEASTOL",
      "   This tolerance determines when a solution is treated as feasible.   If the amount by which a constraint's activity violates its right-hand side or ranged bound   is less in absolute magnitude than FEASTOL, then the constraint is treated as satisfied.   Similarly, if the amount by which a column violates its bounds is less in absolute magnitude   than FEASTOL, those bounds are also treated as satisfied. "
      "\n\nDefault: 1.0E-06",
      XPRS_FEASTOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_FEASTOL

#ifdef XPRS_FEASTOLPERTURB
    MPD( AddSolverOption_MergeDuplicates("sim:xprs_feastolperturb XPRS_FEASTOLPERTURB",
      "   This tolerance determines how much a feasible primal basic solution is   allowed to be perturbed when performing basis changes. The tolerance FEASTOL is always   considered as an upper limit for the perturbations, but in some cases smaller value can be more   desirable. "
      "\n\nDefault: 1.0E-06",
      XPRS_FEASTOLPERTURB, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_FEASTOLPERTURB

#ifdef XPRS_FEASTOLTARGET
    MPD( AddSolverOption_MergeDuplicates("sol:xprs_feastoltarget XPRS_FEASTOLTARGET",
      "This specifies the target feasibility tolerance for the solution refiner."
      "\n\nDefault: 0 — use the value specified by FEASTOL.",
      XPRS_FEASTOLTARGET, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_FEASTOLTARGET

#ifdef XPRS_FORCEOUTPUT
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_forceoutput XPRS_FORCEOUTPUT",
      "Certain names in the problem object may be incompatible with different file formats (such as names containing spaces for LP files). If the optimizer might be unable to read back a problem because of non-standard names, it will first attempt to write it out using an extended naming convention. If the names would not be possible to extend so that they would be reproducible and recognizable, it will give an error message and won't create the file. If the optimizer might be unable to read back a problem because of non-standard names, it will give an error message and won't create the file. This option may be used to force output anyway."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  Check format compatibility, and in case of failure try to extend names so that they are reproducible and recognizable."
      "\n- (1)  Force output using problem names as is."
      "\n- (2)  Always use 'x(' original name ')' in LP files to create a representation that can be read by Xpress. Default for problem having spaces in names"
      "\n- (3)  Substitute spaces by the '_' character in LP files",
      XPRS_FORCEOUTPUT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_FORCEOUTPUT

#ifdef XPRS_FORCEPARALLELDUAL
    MPD( AddSolverOption_MergeDuplicates("lp:xprs_forceparalleldual XPRS_FORCEPARALLELDUAL",
      "Dual simplex: specifies whether the dual simplex solver should always use the parallel simplex algorithm. By default, when using a single thread, the dual simplex solver will execute a dedicated sequential simplex algorithm. "
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  Disabled."
      "\n- (1)  Enabled. Force the dual simplex solver to use the parallel algorithm.",
      XPRS_FORCEPARALLELDUAL, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_FORCEPARALLELDUAL

#ifdef XPRS_GENCONSABSTRANSFORMATION
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_genconsabstransformation XPRS_GENCONSABSTRANSFORMATION",
      "This control specifies the reformulation method for absolute value general constraints at the beginning of the search. "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatic."
      "\n- (0)  Use a formulation based on indicator constraints."
      "\n- (1)  Use a formulation based on SOS1-constraints.",
      XPRS_GENCONSABSTRANSFORMATION, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_GENCONSABSTRANSFORMATION

#ifdef XPRS_GENCONSDUALREDUCTIONS
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_genconsdualreductions XPRS_GENCONSDUALREDUCTIONS",
      "This parameter specifies whether dual reductions should be applied to reduce the number of columns and rows added when transforming general constraints to MIP structs. "
      "\n\n"
      "Values (default: 1):\n"
      "\n- (0)  Disabled. No dual reductions, add columns and rows."
      "\n- (1)  Enabled. Only add neccessary columns and rows, drop those implied by the objective sense.",
      XPRS_GENCONSDUALREDUCTIONS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_GENCONSDUALREDUCTIONS

#ifdef XPRS_GLOBALBOUNDINGBOX
    MPD( AddSolverOption_MergeDuplicates("global:xprs_globalboundingbox XPRS_GLOBALBOUNDINGBOX",
      "If a nonlinear problem cannot be solved due to appearing unbounded, it can automatically be regularized by the application of a bounding box on the variables. If this control is set to a negative value, in a second solving attempt all original variables will be bounded by the absolute value of this control. If set to a positive value, there will be a third solving attempt afterwards, if necessary, in which also all auxiliary variables are bounded by this value."
      "\n\n"
      "Values (default: 1.0E+06):\n"
      "\n- (0)  Disabled. Problem will return unbounded."
      "\n- (n<0)  Enabled. Apply lower and upper bounds of this magnitude to all original variables if initial LP is unbounded."
      "\n- (n>0)  Enabled. Apply lower and upper bounds of this magnitude to all original and auxiliary variables if initial LP and first regularization are unbounded.",
      XPRS_GLOBALBOUNDINGBOX, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_GLOBALBOUNDINGBOX

#ifdef XPRS_GLOBALLSHEURSTRATEGY
    MPD( AddSolverOption_MergeDuplicates("global:xprs_globallsheurstrategy XPRS_GLOBALLSHEURSTRATEGY",
      "                         When integer-feasible (for MINLP, any solution for NLP) but nonlinear-infeasible solutions                         are encountered within a global solve, the integer variables can be fixed and a local solver (as defined                         by the LOCALSOLVER control) can be called on the remaining continuous problem. This                         control defines the frequency and effort of such local solves.                 "
      "\n\n"
      "Values (default:                          -1                 ):\n"
      "\n- (-1)  Automatic selection of the strategy."
      "\n- (0)  Never run a local solver on a partially fixed solution."
      "\n- (1)  Conservative strategy."
      "\n- (2)  Moderate strategy."
      "\n- (3)  Aggressive strategy.",
      XPRS_GLOBALLSHEURSTRATEGY, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_GLOBALLSHEURSTRATEGY

#ifdef XPRS_GLOBALNLPCUTS
    MPD( AddSolverOption_MergeDuplicates("global:xprs_globalnlpcuts XPRS_GLOBALNLPCUTS",
      "         Limit on the number of rounds of outer approximation and convexification cuts generated for the root node, when solving an (MI)NLP to global optimality. "
      "\n\nDefault: -1 — determined automatically.",
      XPRS_GLOBALNLPCUTS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_GLOBALNLPCUTS

#ifdef XPRS_GLOBALNUMINITNLPCUTS
    MPD( AddSolverOption_MergeDuplicates("global:xprs_globalnuminitnlpcuts XPRS_GLOBALNUMINITNLPCUTS",
      "         Specifies the maximum number of tangent cuts when setting up the initial relaxation during a global solve.     By default, the algorithm chooses the number of cuts automatically. Adding more cuts tightens the problem,         resulting in a smaller branch-and-bound tree, at the cost of slowing down each LP solve. "
      "\n\nDefault: -1 — determined automatically.",
      XPRS_GLOBALNUMINITNLPCUTS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_GLOBALNUMINITNLPCUTS

#ifdef XPRS_GLOBALPRESOLVEOBBT
    MPD( AddSolverOption_MergeDuplicates("global:xprs_globalpresolveobbt XPRS_GLOBALPRESOLVEOBBT",
      " Controls the amount of work performed by Optimization-Based Bound Tightening (OBBT). "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatic. The solver decides how much effort goes into OBBT."
      "\n- (0)  Disabled. No OBBT will be performed."
      "\n- (1)  OBBT is run for a small subset of the variables, approximately equal to the square root of the number of total variables used in the solve, i.e., original and auxiliary."
      "\n- (2)  OBBT is run on a larger portion of the variables while keeping the computational effort limited w.r.t. the whole global solve."
      "\n- (3)  OBBT is run on all variables. This is most computationally taxing as a large number of LPs will be solved.",
      XPRS_GLOBALPRESOLVEOBBT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_GLOBALPRESOLVEOBBT

#ifdef XPRS_GLOBALSPATIALBRANCHCUTTINGEFFORT
    MPD( AddSolverOption_MergeDuplicates("global:xprs_globalspatialbranchcuttingeffort XPRS_GLOBALSPATIALBRANCHCUTTINGEFFORT",
      "Limits the effort that is spent on creating cuts during spatial branching."
      "\n\n"
      "Values (default: -1.0):\n"
      "\n- (-1)  The algorithm chooses the effort limit automatically (default)."
      "\n- (0)  Disables cuts on branching entities."
      "\n- (0<n<1)  Relative effort to spend on cutting on branching entities. Higher values lead to more cuts.",
      XPRS_GLOBALSPATIALBRANCHCUTTINGEFFORT, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_GLOBALSPATIALBRANCHCUTTINGEFFORT

#ifdef XPRS_GLOBALSPATIALBRANCHIFPREFERORIG
    MPD( AddSolverOption_MergeDuplicates("global:xprs_globalspatialbranchifpreferorig XPRS_GLOBALSPATIALBRANCHIFPREFERORIG",
      "Whether spatial branchings on original variables should be preferred over branching on auxiliary variables that were introduced by the reformulation of the global solver."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatically determined."
      "\n- (0)  Always consider both original and auxiliary variables."
      "\n- (1)  Always prefer branching on original variables over auxiliaries."
      "\n- (2)  Always prefer branching on auxiliary variables over originals.",
      XPRS_GLOBALSPATIALBRANCHIFPREFERORIG, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_GLOBALSPATIALBRANCHIFPREFERORIG

#ifdef XPRS_GLOBALSPATIALBRANCHPROPAGATIONEFFORT
    MPD( AddSolverOption_MergeDuplicates("global:xprs_globalspatialbranchpropagationeffort XPRS_GLOBALSPATIALBRANCHPROPAGATIONEFFORT",
      "Limits the effort that is spent on propagation during spatial branching. "
      "\n\n"
      "Values (default: -1.0):\n"
      "\n- (-1)  The algorithm chooses the effort limit automatically (default)."
      "\n- (0)  Disables propagation on branching entities."
      "\n- (n>0)  Relative effort to spend on propagating on branching entities. Higher values lead to more propagation.",
      XPRS_GLOBALSPATIALBRANCHPROPAGATIONEFFORT, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_GLOBALSPATIALBRANCHPROPAGATIONEFFORT

#ifdef XPRS_GLOBALTREENLPCUTS
    MPD( AddSolverOption_MergeDuplicates("global:xprs_globaltreenlpcuts XPRS_GLOBALTREENLPCUTS",
      "         Limit on the number of rounds of outer approximation and convexification cuts generated for each node in the tree, when solving an (MI)NLP to global optimality. "
      "\n\nDefault: -1 — determined automatically.",
      XPRS_GLOBALTREENLPCUTS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_GLOBALTREENLPCUTS

#ifdef XPRS_GOMCUTS
    MPD( AddSolverOption_MergeDuplicates("cut:xprs_gomcuts XPRS_GOMCUTS",
      "Branch and Bound: The number of rounds of Gomory or lift-and-project cuts at the root node."
      "\n\nDefault: -1 — determined automatically.",
      XPRS_GOMCUTS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_GOMCUTS

#ifdef XPRS_GPUPLATFORM
    MPD( AddSolverOption_MergeDuplicates("sys:xprs_gpuplatform XPRS_GPUPLATFORM",
      " 			Controls what kind of GPU support is enabled overall in Xpress. Individual solver components may disable GPU support. 		"
      "\n\n"
      "Values (default:  			1, use GPU support if available, unless disabled by another control 		):\n"
      "\n- (0)  Do not use any GPU support."
      "\n- (1)  Use GPU support based on CUDA (if available).",
      XPRS_GPUPLATFORM, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_GPUPLATFORM

#ifdef XPRS_HEURBEFORELP
    MPD( AddSolverOption_MergeDuplicates("heur:xprs_heurbeforelp XPRS_HEURBEFORELP",
      "Branch and Bound: Determines whether primal heuristics should be run before the initial LP relaxation has been solved."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatic - let the optimizer decide if heuristics should be run."
      "\n- (0)  Disabled."
      "\n- (1)  Enabled.",
      XPRS_HEURBEFORELP, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_HEURBEFORELP

#ifdef XPRS_HEURDEPTH
    MPD( AddSolverOption_MergeDuplicates("heur:xprs_heurdepth XPRS_HEURDEPTH",
      "Branch and Bound: Sets the maximum depth in the tree search at which heuristics will be used to find MIP solutions. It may be worth stopping the heuristic search for solutions after a certain depth in the tree search. A value of 0 signifies that heuristics will not be used.  This control no longer has any effect and will be removed from future releases."
      "\n\nDefault: -1",
      XPRS_HEURDEPTH, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_HEURDEPTH

#ifdef XPRS_HEURDIVEITERLIMIT
    MPD( AddSolverOption_MergeDuplicates("heur:xprs_heurdiveiterlimit XPRS_HEURDIVEITERLIMIT",
      "Branch and Bound: Simplex iteration limit for reoptimizing during the diving heuristic."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (>=1)  Fixed iteration limit."
      "\n- (0)  No iteration limit."
      "\n- (<0)  Automatic selection of the iteration limit based on the problem size. The absolute value is used as a multiplier on the automatic selection.",
      XPRS_HEURDIVEITERLIMIT, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_HEURDIVEITERLIMIT

#ifdef XPRS_HEURDIVERANDOMIZE
    MPD( AddSolverOption_MergeDuplicates("heur:xprs_heurdiverandomize XPRS_HEURDIVERANDOMIZE",
      "The level of randomization to apply in the diving heuristic. The diving heuristic uses priority weights on rows and columns to determine the order in which to e.g. round fractional columns, or the direction in which to round them. This control determines by how large a random factor these weights should be changed."
      "\n\n"
      "Values (default: 0.0):\n"
      "\n- (0.0-1.0)  Amount of randomization (0.0=none, 1.0=full)",
      XPRS_HEURDIVERANDOMIZE, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_HEURDIVERANDOMIZE

#ifdef XPRS_HEURDIVESOFTROUNDING
    MPD( AddSolverOption_MergeDuplicates("heur:xprs_heurdivesoftrounding XPRS_HEURDIVESOFTROUNDING",
      "Branch and Bound: Enables a more cautious strategy for the diving heuristic, where it tries to push binaries and integer variables to their bounds using the objective, instead of directly fixing them. This can be useful when the default diving heuristics fail to find any feasible solutions."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatic selection."
      "\n- (0)  Do not use soft rounding."
      "\n- (1)  Cautious use of the soft rounding strategy."
      "\n- (2)  More aggressive use of the soft rounding strategy.",
      XPRS_HEURDIVESOFTROUNDING, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_HEURDIVESOFTROUNDING

#ifdef XPRS_HEURDIVESPEEDUP
    MPD( AddSolverOption_MergeDuplicates("heur:xprs_heurdivespeedup XPRS_HEURDIVESPEEDUP",
      "Branch and Bound: Changes the emphasis of the diving heuristic from solution quality to diving speed."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-2)  Automatic selection biased towards quality"
      "\n- (-1)  Automatic selection biased towards speed."
      "\n- (0-4)  manual emphasis bias from emphasis on quality (0) to emphasis on speed (4).",
      XPRS_HEURDIVESPEEDUP, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_HEURDIVESPEEDUP

#ifdef XPRS_HEURDIVESTRATEGY
    MPD( AddSolverOption_MergeDuplicates("heur:xprs_heurdivestrategy XPRS_HEURDIVESTRATEGY",
      "Branch and Bound: Chooses the strategy for the diving heuristic."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatic selection of strategy."
      "\n- (0)  Disables the diving heuristic."
      "\n- (1-18)  Available pre-set strategies for rounding infeasible MIP entities and reoptimizing during the heuristic dive.",
      XPRS_HEURDIVESTRATEGY, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_HEURDIVESTRATEGY

#ifdef XPRS_HEUREMPHASIS
    MPD( AddSolverOption_MergeDuplicates("heur:xprs_heuremphasis XPRS_HEUREMPHASIS",
      "Branch and Bound: This control specifies an emphasis for the search w.r.t. primal heuristics and   other procedures that affect the speed of convergence of the primal-dual gap.   For problems where the goal is to achieve a small gap but not neccessarily solving them to optimality,   it is recommended to set HEUREMPHASIS to 1.   This setting triggers many additional heuristic calls, aiming for reducing the gap at the beginning   of the search, typically at the expense of an increased time for proving optimality. "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Optimizer default strategy."
      "\n- (0)  Disables all heuristics."
      "\n- (1)  Focus on reducing the primal-dual gap in the early part of the search."
      "\n- (2)  Extremely aggressive search heuristics.",
      XPRS_HEUREMPHASIS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_HEUREMPHASIS

#ifdef XPRS_HEURFORCESPECIALOBJ
    MPD( AddSolverOption_MergeDuplicates("heur:xprs_heurforcespecialobj XPRS_HEURFORCESPECIALOBJ",
      "Branch and Bound: This specifies whether local search heuristics without objective or with an auxiliary objective should always be used, despite the automatic selection of the Optimiezr. Deactivated by default. "
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  Disabled."
      "\n- (1)  Enabled. Run special objective heuristics on large problems and even if incumbent exists.",
      XPRS_HEURFORCESPECIALOBJ, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_HEURFORCESPECIALOBJ

#ifdef XPRS_HEURFREQ
    MPD( AddSolverOption_MergeDuplicates("heur:xprs_heurfreq XPRS_HEURFREQ",
      "Branch and Bound: This specifies the frequency at which heuristics are used in the tree search. Heuristics will only be used at a node if the depth of the node is a multiple of HEURFREQ."
      "\n\nDefault: -1",
      XPRS_HEURFREQ, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_HEURFREQ

#ifdef XPRS_HEURMAXSOL
    MPD( AddSolverOption_MergeDuplicates("heur:xprs_heurmaxsol XPRS_HEURMAXSOL",
      "Branch and Bound: This specifies the maximum number of heuristic solutions that will be found in the tree search.  This control no longer has any effect and will be removed from future releases."
      "\n\nDefault: -1",
      XPRS_HEURMAXSOL, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_HEURMAXSOL

#ifdef XPRS_HEURNODES
    MPD( AddSolverOption_MergeDuplicates("heur:xprs_heurnodes XPRS_HEURNODES",
      "Branch and Bound: This specifies the maximum number of nodes at which heuristics are used in the tree search.  This control no longer has any effect and will be removed from future releases."
      "\n\nDefault: -1",
      XPRS_HEURNODES, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_HEURNODES

#ifdef XPRS_HEURSEARCHBACKGROUNDSELECT
    MPD( AddSolverOption_MergeDuplicates("heur:xprs_heursearchbackgroundselect XPRS_HEURSEARCHBACKGROUNDSELECT",
      "     Bit-vector control (see Section Bit-vector controls) to select which large neighborhood searches to run in the background         (for example in parallel to the root cut loop).   "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (1)  Enable L heuristic.",
      XPRS_HEURSEARCHBACKGROUNDSELECT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_HEURSEARCHBACKGROUNDSELECT

#ifdef XPRS_HEURSEARCHCOPYCONTROLS
    MPD( AddSolverOption_MergeDuplicates("heur:xprs_heursearchcopycontrols XPRS_HEURSEARCHCOPYCONTROLS",
      "     Select how user-set controls should affect local search heuristics.   "
      "\n\n"
      "Values (default: 1):\n"
      "\n- (0)  Do not copy any user-set controls"
      "\n- (1)  Automatic - Let the Optimizer decide which user-set controls to copy"
      "\n- (2)  Copy all user-set controls",
      XPRS_HEURSEARCHCOPYCONTROLS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_HEURSEARCHCOPYCONTROLS

#ifdef XPRS_HEURSEARCHEFFORT
    MPD( AddSolverOption_MergeDuplicates("heur:xprs_heursearcheffort XPRS_HEURSEARCHEFFORT",
      "Adjusts the overall level of the local search heuristics."
      "\n\nDefault: 1.0",
      XPRS_HEURSEARCHEFFORT, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_HEURSEARCHEFFORT

#ifdef XPRS_HEURSEARCHFREQ
    MPD( AddSolverOption_MergeDuplicates("heur:xprs_heursearchfreq XPRS_HEURSEARCHFREQ",
      "Branch and Bound: This specifies how often the local search heuristic should be run in the tree."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatic."
      "\n- (0)  Disabled in the tree."
      "\n- (n>0)  Number of nodes between each run.",
      XPRS_HEURSEARCHFREQ, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_HEURSEARCHFREQ

#ifdef XPRS_HEURSEARCHROOTCUTFREQ
    MPD( AddSolverOption_MergeDuplicates("heur:xprs_heursearchrootcutfreq XPRS_HEURSEARCHROOTCUTFREQ",
      " How frequently to run the local search heuristic during root cutting. This is given as how many cut rounds to perform between runs of the heuristic.  Set to zero to avoid applying the heuristic during root cutting.   Branch and Bound: This specifies how often the local search heuristic should be run in the tree."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatic."
      "\n- (0)  Disabled heuristic during cutting."
      "\n- (n>0)  Number of cutting rounds between each run.",
      XPRS_HEURSEARCHROOTCUTFREQ, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_HEURSEARCHROOTCUTFREQ

#ifdef XPRS_HEURSEARCHROOTSELECT
    MPD( AddSolverOption_MergeDuplicates("heur:xprs_heursearchrootselect XPRS_HEURSEARCHROOTSELECT",
      "A bit-vector control (see Section Bit-vector controls) for selecting which local search heuristics to apply on the root node of a MIP solve. Use HEURSEARCHTREESELECT to control local search heuristics during the tree search."
      "\n\n"
      "Values (default: 117):\n"
      "\n- (0)  Local search with a large neighborhood. Potentially slow but is good for finding solutions that  differs significantly from the incumbent."
      "\n- (1)  Local search with a small neighborhood centered around a node LP solution."
      "\n- (2)  Local search with a small neighborhood centered around an integer solution. This heuristic will  often provide smaller, incremental improvements to an incumbent solution."
      "\n- (3)  Local search with a neighborhood set up through the combination of multiple integer solutions."
      "\n- (4)  Unused"
      "\n- (5)  Local search without an objective function. Called seldom and only when no feasible solution is available."
      "\n- (6)  Local search with an auxiliary objective function. Called seldom and only when no feasible solution is available.",
      XPRS_HEURSEARCHROOTSELECT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_HEURSEARCHROOTSELECT

#ifdef XPRS_HEURSEARCHTREESELECT
    MPD( AddSolverOption_MergeDuplicates("heur:xprs_heursearchtreeselect XPRS_HEURSEARCHTREESELECT",
      "A bit-vector control (see Section Bit-vector controls) for selecting which local search heuristics to apply during the tree search of a MIP solve. Use HEURSEARCHROOTSELECT to control local search heuristics on the root node."
      "\n\n"
      "Values (default: 17):\n"
      "\n- (0)  Local search with a large neighborhood. Potentially slow but is good for finding solutions that  differs significantly from the incumbent."
      "\n- (1)  Local search with a small neighborhood centered around a node LP solution."
      "\n- (2)  Local search with a small neighborhood centered around an integer solution. This heuristic will  often provide smaller, incremental improvements to an incumbent solution."
      "\n- (3)  Local search with a neighborhood set up through the combination of multiple integer solutions."
      "\n- (4)  Unused"
      "\n- (5)  Local search without an objective function. Called seldom and only when no feasible solution is available."
      "\n- (6)  Local search with an auxiliary objective function. Called seldom and only when no feasible solution is available.",
      XPRS_HEURSEARCHTREESELECT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_HEURSEARCHTREESELECT

#ifdef XPRS_HEURSHIFTPROP
    MPD( AddSolverOption_MergeDuplicates("heur:xprs_heurshiftprop XPRS_HEURSHIFTPROP",
      "Determines whether the Shift-and-propagate primal heuristic should be executed. If enabled, Shift-and-propagate is an LP-free primal heuristic that is executed immediately after presolve."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  The solver decides if Shift-and-propagate should be run. This is the default setting."
      "\n- (0)  Shift-and-propagate is disabled."
      "\n- (1)  Shift-and-propagate is enabled.",
      XPRS_HEURSHIFTPROP, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_HEURSHIFTPROP

#ifdef XPRS_HEURTHREADS
    MPD( AddSolverOption_MergeDuplicates("heur:xprs_heurthreads XPRS_HEURTHREADS",
      "Branch and Bound: The number of threads to dedicate to running heuristics during the root solve."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (-1)  Automatically determined from the THREADS control."
      "\n- (0)  Disabled."
      "\n- (>=1)  Number of additional threads to dedicate to parallel heuristics.",
      XPRS_HEURTHREADS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_HEURTHREADS

#ifdef XPRS_HISTORYCOSTS
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_historycosts XPRS_HISTORYCOSTS",
      "Branch and Bound: How to update the pseudo cost for a MIP entity when a strong branch or a regular branch is applied."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatically determined."
      "\n- (0)  No update."
      "\n- (1)  Update using only regular branches from the root to the current node."
      "\n- (2)  Same as 1, but update with strong branching results as well."
      "\n- (3)  Update using any regular branching or strong branching information from all nodes solves before the current node.",
      XPRS_HISTORYCOSTS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_HISTORYCOSTS

#ifdef XPRS_IFCHECKCONVEXITY
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_ifcheckconvexity XPRS_IFCHECKCONVEXITY",
      "Determines if the convexity of the problem is checked before optimization. Applies to quadratic, mixed integer quadratic and quadratically constrained problems. Checking convexity takes some time, thus for problems that are known to be convex it might be reasonable to switch the checking off."
      "\n\n"
      "Values (default: 1):\n"
      "\n- (0)  Turn off convexity checking."
      "\n- (1)  Turn on convexity checking.",
      XPRS_IFCHECKCONVEXITY, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_IFCHECKCONVEXITY

#ifdef XPRS_IISLOG
    MPD( AddSolverOption_MergeDuplicates("inf:xprs_iislog XPRS_IISLOG",
      "                         Selects how much information should be printed during the IIS procedure. Please refer to Appendix The IIS Log for a more detailed description of the IIS logging format.                 "
      "\n\n"
      "Values (default:                  1, a progress log is printed         ):\n"
      "\n- (0)  The IIS procedure does not produce any output."
      "\n- (1)  Prints general information and a progress log of the deletion filter, including bounds on the size of the IIS and timing information."
      "\n- (2)  Complete logging, including the full progress log of all the subproblem solves in the deletion filter. This setting is recommended only for debugging as it may produce a lot of output.",
      XPRS_IISLOG, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_IISLOG

#ifdef XPRS_IISOPS
    MPD( AddSolverOption_MergeDuplicates("inf:xprs_iisops XPRS_IISOPS",
      "                 Selects which part of the restrictions (bounds, constraints, entities) should always be kept in an IIS. This is useful if certain types of restrictions cannot be violated, thus they are known not to be the cause of infeasibility.                 The IIS obtained this way is irreducible only for the non-protected restrictions.                  This bit-vector control (see Section Bit-vector controls) has an effect only on the deletion filter of the IIS procedure.         "
      "\n\n"
      "Values (default:                  0, all restrictions are valid candidates for removal         ):\n"
      "\n- (0)  Keep binary integralities."
      "\n- (1)  Keep the 0 lower bounds of variables."
      "\n- (2)  Keep fixed variables."
      "\n- (3)  Keep all variable bounds."
      "\n- (4)  Keep all general integer entities, except binaries."
      "\n- (5)  Keep all equality constraints."
      "\n- (6)  Keep all general constraints."
      "\n- (7)  Keep all piecewise linear constraints."
      "\n- (8)  Keep all specially ordered set (SOS) constraints."
      "\n- (9)  Keep all indicator constraints."
      "\n- (10)  Keep all delayed rows."
      "\n- (11)  Keep all constraints.",
      XPRS_IISOPS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_IISOPS

#ifdef XPRS_INDLINBIGM
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_indlinbigm XPRS_INDLINBIGM",
      "During presolve, indicator constraints will be linearized using a BigM coefficient whenever that BigM coefficient is small enough. This control defines the largest BigM for which such a linearized version will be added to the problem in addition to the original constraint. If the BigM is even smaller than INDPRELINBIGM, then the original indicator constraint will additionally be dropped from the problem. "
      "\n\nDefault: 1.0E+05",
      XPRS_INDLINBIGM, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_INDLINBIGM

#ifdef XPRS_INDPRELINBIGM
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_indprelinbigm XPRS_INDPRELINBIGM",
      "During presolve, indicator constraints will be linearized using a BigM coefficient whenever that BigM coefficient is small enough. This control defines the largest BigM for which the original constraint will be replaced by the linearized version. If the BigM is larger than INDPRELINBIGM but smaller than INDLINBIGM, the linearized row will be added but the original indicator constraint is kept as a numerically stable way to check feasibility. "
      "\n\nDefault: 100.0",
      XPRS_INDPRELINBIGM, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_INDPRELINBIGM

#ifdef XPRS_INPUTTOL
    MPD( AddSolverOption_MergeDuplicates("prob:xprs_inputtol XPRS_INPUTTOL",
      "The tolerance on input values elements. If any value is less than or equal to INPUTTOL in absolute value, it is treated as zero. For the internal zero tolerance see MATRIXTOL. "
      "\n\nDefault: 0.0",
      XPRS_INPUTTOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_INPUTTOL

#ifdef XPRS_INVERTFREQ
    MPD( AddSolverOption_MergeDuplicates("sim:xprs_invertfreq XPRS_INVERTFREQ",
      "Simplex: The frequency with which the basis will be inverted. The basis is maintained in a factorized form and on most simplex iterations it is incrementally updated to reflect the step just taken. This is considerably faster than computing the full inverted matrix at each iteration, although after a number of iterations the basis becomes less well-conditioned and it becomes necessary to compute the full inverted matrix. The value of INVERTFREQ specifies the maximum number of iterations between full inversions."
      "\n\nDefault: -1 — the frequency is determined automatically.",
      XPRS_INVERTFREQ, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_INVERTFREQ

#ifdef XPRS_INVERTMIN
    MPD( AddSolverOption_MergeDuplicates("sim:xprs_invertmin XPRS_INVERTMIN",
      "Simplex: The minimum number of iterations between full inversions of the basis matrix. See the description of INVERTFREQ for details."
      "\n\nDefault: 3",
      XPRS_INVERTMIN, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_INVERTMIN

#ifdef XPRS_IOTIMEOUT
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_iotimeout XPRS_IOTIMEOUT",
      "The maximum number of seconds to wait for an I/O operation before it is cancelled."
      "\n\nDefault: 30",
      XPRS_IOTIMEOUT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_IOTIMEOUT

#ifdef XPRS_KEEPBASIS
    MPD( AddSolverOption_MergeDuplicates("lp:xprs_keepbasis XPRS_KEEPBASIS",
      "Simplex: This determines whether the basis should be kept when reoptimizing a problem.  The choice is between using a crash basis created at the beginning of simplex or using a basis from a previous solve (if such exists). By default, this control gets (re)set automatically in various situations. By default, it will be automatically set to 1 after a solve that produced a valid basis. This will automatically warmstart a subsequent solve. Explicitly loading a starting basis will also set this control to 1. If the control is explicitly set to 0, any existing basis will be ignored for a new solve, and the Optimizer will start from an ad-hoc crash basis. "
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  Problem optimization starts from scratch, i.e., any previous basis is ignored."
      "\n- (1)  The previous basis should be used as a starting basis."
      "\n- (2)  Use the previous basis only if it is valid for the current problem (the number of basic variables must match the number of rows)."
      "\n- (3)  Use the previous basis only if it is valid and numerically stable in the current problem.",
      XPRS_KEEPBASIS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_KEEPBASIS

#ifdef XPRS_KEEPNROWS
    MPD( AddSolverOption_MergeDuplicates("prob:xprs_keepnrows XPRS_KEEPNROWS",
      "How nonbinding rows should be handled by the MPS reader."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Delete N type rows from the matrix."
      "\n- (0)  Delete elements from N type rows leaving empty N type rows in the matrix."
      "\n- (1)  Keep N type rows.",
      XPRS_KEEPNROWS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_KEEPNROWS

#ifdef XPRS_L1CACHE
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_l1cache XPRS_L1CACHE",
      "This parameter is deprecated and will be removed in a future release. Newton barrier: L1 cache size in kB (kilo bytes) of the CPU. On Intel (or compatible) platforms a value of -1 may be used to determine the cache size automatically. "
      "\n\nDefault: Hardware/platform dependent.",
      XPRS_L1CACHE, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_L1CACHE

#ifdef XPRS_LNPBEST
    MPD( AddSolverOption_MergeDuplicates("cut:xprs_lnpbest XPRS_LNPBEST",
      "Number of infeasible MIP entities to create lift-and-project cuts for during each round of Gomory cuts at the root node (see GOMCUTS)."
      "\n\nDefault: 50",
      XPRS_LNPBEST, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_LNPBEST

#ifdef XPRS_LNPITERLIMIT
    MPD( AddSolverOption_MergeDuplicates("cut:xprs_lnpiterlimit XPRS_LNPITERLIMIT",
      "Number of iterations to perform in improving each lift-and-project cut."
      "\n\nDefault: -1 — determined automatically.",
      XPRS_LNPITERLIMIT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_LNPITERLIMIT

#ifdef XPRS_LOCALCHOICE
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_localchoice XPRS_LOCALCHOICE",
      "Controls when to perform a local backtrack between the two child nodes during a dive in the branch and bound tree."
      "\n\n"
      "Values (default: 3):\n"
      "\n- (1)  Never backtrack from the first child, unless it is dropped (infeasible or cut off)."
      "\n- (2)  Always solve both child nodes before deciding which child to continue with."
      "\n- (3)  Automatically determined. ",
      XPRS_LOCALCHOICE, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_LOCALCHOICE

#ifdef XPRS_LPFLAGS
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_lpflags XPRS_LPFLAGS",
      "A bit-vector control (see Section Bit-vector controls) which defines the algorithm for solving an LP problem or the initial LP relaxation of a MIP problem."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  Use the dual simplex method."
      "\n- (1)  Use the primal simplex method."
      "\n- (2)  Use the barrier method (or hybrid gradient method if BARALG=4 is set)."
      "\n- (3)  Use the network simplex method.",
      XPRS_LPFLAGS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_LPFLAGS

#ifdef XPRS_LPFOLDING
    MPD( AddSolverOption_MergeDuplicates("lp:xprs_lpfolding XPRS_LPFOLDING",
      "Simplex and barrier: whether to fold an LP problem before solving it. "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatic."
      "\n- (0)  Disable LP folding."
      "\n- (1)  Enable LP folding. Attempt to fold all LP problems and MIP initial relaxations. ",
      XPRS_LPFOLDING, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_LPFOLDING

#ifdef XPRS_LPITERLIMIT
    MPD( AddSolverOption_MergeDuplicates("lp:xprs_lpiterlimit XPRS_LPITERLIMIT",
      "The maximum number of iterations that will be performed by primal simplex or dual simplex before the optimization process terminates. For MIP problems, this is the maximum total number of iterations over all nodes explored by the Branch and Bound method."
      "\n\nDefault: 2147483647",
      XPRS_LPITERLIMIT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_LPITERLIMIT

#ifdef XPRS_LPLOG
    MPD( AddSolverOption_MergeDuplicates("lp:xprs_lplog XPRS_LPLOG",
      "Simplex: The frequency at which the simplex log is printed. "
      "\n\n"
      "Values (default: 100):\n"
      "\n- (n<0)  Detailed output every -n iterations."
      "\n- (0)  Log displayed at the end of the optimization only."
      "\n- (n>0)  Summary output every n iterations.",
      XPRS_LPLOG, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_LPLOG

#ifdef XPRS_LPLOGDELAY
    MPD( AddSolverOption_MergeDuplicates("lp:xprs_lplogdelay XPRS_LPLOGDELAY",
      "Time interval between two LP log lines. "
      "\n\nDefault: 1.0",
      XPRS_LPLOGDELAY, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_LPLOGDELAY

#ifdef XPRS_LPLOGSTYLE
    MPD( AddSolverOption_MergeDuplicates("lp:xprs_lplogstyle XPRS_LPLOGSTYLE",
      "Simplex: The style of the simplex log."
      "\n\n"
      "Values (default: 1):\n"
      "\n- (0)  Simplex log is printed based on simplex iteration count, at a fixed frequency as specified by the LPLOG control."
      "\n- (1)  Simplex log is printed based on an estimation of elapsed time, determined by an internal deterministic timer.",
      XPRS_LPLOGSTYLE, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_LPLOGSTYLE

#ifdef XPRS_LPREFINEITERLIMIT
    MPD( AddSolverOption_MergeDuplicates("sol:xprs_lprefineiterlimit XPRS_LPREFINEITERLIMIT",
      "This specifies the simplex iteration limit the solution refiner can spend in attempting to increase the accuracy of an LP solution."
      "\n\nDefault: -1 — determined automatically.",
      XPRS_LPREFINEITERLIMIT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_LPREFINEITERLIMIT

#ifdef XPRS_MARKOWITZTOL
    MPD( AddSolverOption_MergeDuplicates("lp:xprs_markowitztol XPRS_MARKOWITZTOL",
      "The Markowitz tolerance used for the factorization of the basis matrix."
      "\n\nDefault: 0.01",
      XPRS_MARKOWITZTOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_MARKOWITZTOL

#ifdef XPRS_MATRIXTOL
    MPD( AddSolverOption_MergeDuplicates("tol:xprs_matrixtol XPRS_MATRIXTOL",
      "The zero tolerance on matrix elements. If the value of a matrix element is less than or equal to MATRIXTOL in absolute value, it is treated as zero. The control applies when solving a problem, for an input tolerance see INPUTTOL. "
      "\n\nDefault: 1.0E-09",
      XPRS_MATRIXTOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_MATRIXTOL

#ifdef XPRS_MAXCHECKSONMAXCUTTIME
    MPD( AddSolverOption_MergeDuplicates("lim:xprs_maxchecksonmaxcuttime XPRS_MAXCHECKSONMAXCUTTIME",
      "This control is intended for use where optimization runs that are terminated using the MAXCUTTIME control are  required to be reproduced exactly. This control is necessary because of the inherent difficulty in terminating algorithmic  software in a consistent way using temporal criteria. The control value relates to the number of times the optimizer checks  the MAXCUTTIME criterion up to and including the check when the termination of cutting was activated. To use the control the user first must  obtain the value of the CHECKSONMAXCUTTIME attribute after the run returns. This attribute value is the number of times the  optimizer checked the MAXCUTTIME criterion during the last call to the optimization routine  XPRSmipoptimize. Note that this attribute value will be negative if the  optimizer terminated cutting on the MAXCUTTIME criterion. To ensure accurate reproduction of a run the user should first  ensure that MAXCUTTIME is set to its default value or to a large value so the run does not terminate again on MAXCUTTIME  and then simply set the control MAXCHECKSONMAXCUTTIME to the absolute value of the CHECKSONMAXCUTTIME value. "
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  Not active."
      "\n- (n>0)  The number of times the optimizer should check the MAXCUTTIME criterion before triggering a termination.",
      XPRS_MAXCHECKSONMAXCUTTIME, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MAXCHECKSONMAXCUTTIME

#ifdef XPRS_MAXCHECKSONMAXTIME
    MPD( AddSolverOption_MergeDuplicates("lim:xprs_maxchecksonmaxtime XPRS_MAXCHECKSONMAXTIME",
      "This control is intended for use where optimization runs that are terminated using the TIMELIMIT (or the deprecated MAXTIME) control are required to be reproduced exactly. This control is necessary because of the inherent difficulty in terminating algorithmic software in a consistent way using temporal criteria. The control value relates to the number of times the optimizer checks the TIMELIMIT criterion up to and including the check when the termination was activated. To use the control the user first must obtain the value of the CHECKSONMAXTIME attribute after the run returns. This attribute value is the number of times the optimizer checked the TIMELIMIT criterion during the last call to the optimization routine XPRSmipoptimize. Note that this attribute value will be negative if the optimizer terminated on the TIMELIMIT criterion. To ensure that a reproduction of a run terminates in the same way the user should first ensure that TIMELIMIT is set to its default value or to a large value so the run does not terminate again on TIMELIMIT and then simply set the control MAXCHECKSONMAXTIME to the absolute value of the CHECKSONMAXTIME value. "
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  Not active."
      "\n- (n>0)  The number of times the optimizer should check the TIMELIMIT (or MAXTIME) criterion before triggering a termination.",
      XPRS_MAXCHECKSONMAXTIME, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MAXCHECKSONMAXTIME

#ifdef XPRS_MAXCUTTIME
    MPD( AddSolverOption_MergeDuplicates("cut:xprs_maxcuttime XPRS_MAXCUTTIME",
      "The maximum amount of time allowed for generation of cutting planes and reoptimization. The limit is checked during generation and no further cuts are added once this limit has been exceeded."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  No time limit."
      "\n- (>0)  Stop cut generation after the given number of seconds.",
      XPRS_MAXCUTTIME, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_MAXCUTTIME

#ifdef XPRS_MAXIIS
    MPD( AddSolverOption_MergeDuplicates("inf:xprs_maxiis XPRS_MAXIIS",
      "This function controls the number of Irreducible Infeasible Sets to be found using the XPRSiisall (IIS-a)."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Search for all IIS."
      "\n- (0)  Do not search for IIS."
      "\n- (n>0)  Search for the first n IIS.",
      XPRS_MAXIIS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MAXIIS

#ifdef XPRS_MAXIMPLIEDBOUND
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_maximpliedbound XPRS_MAXIMPLIEDBOUND",
      "Presolve: When tighter bounds are calculated during MIP preprocessing, only bounds whose absolute value are smaller than MAXIMPLIEDBOUND will be applied to the problem."
      "\n\nDefault: 1.0E+08",
      XPRS_MAXIMPLIEDBOUND, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_MAXIMPLIEDBOUND

#ifdef XPRS_MAXLOCALBACKTRACK
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_maxlocalbacktrack XPRS_MAXLOCALBACKTRACK",
      "Branch-and-Bound: How far back up the current dive path the optimizer is allowed to look for a local backtrack candidate node."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatic."
      "\n- (n>0)  Local backtrack limit.",
      XPRS_MAXLOCALBACKTRACK, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MAXLOCALBACKTRACK

#ifdef XPRS_MAXMCOEFFBUFFERELEMS
    MPD( AddSolverOption_MergeDuplicates("prob:xprs_maxmcoeffbufferelems XPRS_MAXMCOEFFBUFFERELEMS",
      "The maximum number of matrix coefficients to buffer before flushing into the internal representation of the problem.  Buffering coefficients can offer a significant performance gain when you are building a matrix using XPRSchgcoef or XPRSchgmcoef, but can lead to a significant memory overhead for large matrices, which this control allows you to influence. "
      "\n\nDefault: 2147483647",
      XPRS_MAXMCOEFFBUFFERELEMS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MAXMCOEFFBUFFERELEMS

#ifdef XPRS_MAXMEMORYHARD
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_maxmemoryhard XPRS_MAXMEMORYHARD",
      "This control sets the maximum amount of memory in megabytes the optimizer should allocate. If this limit is exceeded, the solve will terminate. This control is designed to make the optimizer stop in a controlled manner, so that the problem object is valid once termination occurs. The solve state will be set to incomplete. This is different to an out of memory condition in which case the optimizer returns an error. The optimizer may still allocate memory once the limit is exceeded to be able to finsish the operations and stop in a controlled manner. When RESOURCESTRATEGY is enabled, the control also has the same effect as MAXMEMORYSOFT and will cause the optimizer to try preserving memory when possible."
      "\n\nDefault: 0 (no limit)",
      XPRS_MAXMEMORYHARD, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MAXMEMORYHARD

#ifdef XPRS_MAXMEMORYSOFT
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_maxmemorysoft XPRS_MAXMEMORYSOFT",
      "When RESOURCESTRATEGY is enabled, this control sets the maximum amount of memory in megabytes the optimizer targets to allocate. This may change the solving path, but will not cause the solve to terminate early. To set a hard version of the same, please set MAXMEMORYHARD."
      "\n\nDefault: 0 (no limit)",
      XPRS_MAXMEMORYSOFT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MAXMEMORYSOFT

#ifdef XPRS_MAXMIPSOL
    MPD( AddSolverOption_MergeDuplicates("lim:xprs_maxmipsol XPRS_MAXMIPSOL",
      "Branch and Bound: This specifies a limit on the number of integer solutions to be found by the Optimizer. It is possible that during optimization the Optimizer will find the same objective solution from different nodes. However, MAXMIPSOL refers to the total number of integer solutions found, and not necessarily the number of distinct solutions."
      "\n\nDefault: 0",
      XPRS_MAXMIPSOL, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MAXMIPSOL

#ifdef XPRS_MAXMIPTASKS
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_maxmiptasks XPRS_MAXMIPTASKS",
      "Branch-and-Bound: The maximum number of tasks to run in parallel during a MIP solve."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Task limit determined automatically from MIPTHREADS."
      "\n- (>0)  Fixed task limit.",
      XPRS_MAXMIPTASKS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MAXMIPTASKS

#ifdef XPRS_MAXNODE
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_maxnode XPRS_MAXNODE",
      "Branch and Bound: The maximum number of nodes that will be explored."
      "\n\nDefault: 2147483647",
      XPRS_MAXNODE, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MAXNODE

#ifdef XPRS_MAXPAGELINES
    MPD( AddSolverOption_MergeDuplicates("log:xprs_maxpagelines XPRS_MAXPAGELINES",
      "Number of lines between page breaks in printable output."
      "\n\nDefault: 23",
      XPRS_MAXPAGELINES, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MAXPAGELINES

#ifdef XPRS_MAXSCALEFACTOR
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_maxscalefactor XPRS_MAXSCALEFACTOR",
      "This determines the maximum scaling factor that can be applied during scaling. The maximum is provided as an exponent of a power of 2."
      "\n\n"
      "Values (default: 64):\n"
      "\n- (0-64)  The maximum is provided an exponent of a power of 2.",
      XPRS_MAXSCALEFACTOR, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MAXSCALEFACTOR

#ifdef XPRS_MAXSTALLTIME
    MPD( AddSolverOption_MergeDuplicates("lim:xprs_maxstalltime XPRS_MAXSTALLTIME",
      "The maximum time in seconds that the Optimizer will continue to search for improving solution after finding a new incumbent."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  No stall time limit."
      "\n- (>0)  If an integer solution has been found, stop MIP search after the given number of seconds without a new incumbent. No effect as long as no solution was found.",
      XPRS_MAXSTALLTIME, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_MAXSTALLTIME

#ifdef XPRS_MAXTIME
    MPD( AddSolverOption_MergeDuplicates("lim:xprs_maxtime XPRS_MAXTIME",
      "This parameter is deprecated and will be removed in a future release. The maximum time in seconds that the Optimizer will run before it terminates, including the problem setup time and solution time. For MIP problems, this is the total time taken to solve all nodes. "
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  No time limit."
      "\n- (n>0)  If an integer solution has been found, stop MIP search after n seconds, otherwise continue until an integer solution is finally found."
      "\n- (n<0)  Stop in LP or MIP search after n seconds.",
      XPRS_MAXTIME, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MAXTIME

#ifdef XPRS_MAXTREEFILESIZE
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_maxtreefilesize XPRS_MAXTREEFILESIZE",
      "The maximum size, in megabytes, to which the tree file may grow, or 0 for no limit.  When the tree file reaches this limit, a second tree file will be created.  Useful if you are using a filesystem that puts a maximum limit on the size of a file. "
      "\n\nDefault: 0",
      XPRS_MAXTREEFILESIZE, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MAXTREEFILESIZE

#ifdef XPRS_MCFCUTSTRATEGY
    MPD( AddSolverOption_MergeDuplicates("cut:xprs_mcfcutstrategy XPRS_MCFCUTSTRATEGY",
      " Level of Multi-Commodity Flow (MCF) cutting planes separation: This specifies how aggressively MCF cuts should be separated. If the separation of MCF cuts is enabled, Xpress will try to detect a MCF network structure in the problem and, if such a structure is identified, it will separate specific cutting planes exploiting the identified network. "
      "\n\n"
      "Values (default:  -1 ):\n"
      "\n- (-1)  Automatic - let the Optimizer decide."
      "\n- (0)  Separation of MCF cuts disabled."
      "\n- (1)  Moderate separation of MCF cuts."
      "\n- (2)  Aggressive separation of MCF cuts.",
      XPRS_MCFCUTSTRATEGY, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MCFCUTSTRATEGY

#ifdef XPRS_MIPABSCUTOFF
    MPD( AddSolverOption_MergeDuplicates("lim:xprs_mipabscutoff XPRS_MIPABSCUTOFF",
      "Branch and Bound: If the user knows that they are interested only in values of the objective function which are better than some value, this can be assigned to MIPABSCUTOFF. This allows the Optimizer to ignore solving any nodes which may yield worse objective values, saving solution time. When a MIP solution is found a new cut off value is calculated and the value can be obtained from the CURRMIPCUTOFF attribute.  The value of CURRMIPCUTOFF is calculated using the MIPRELCUTOFF and MIPADDCUTOFF controls."
      "\n\nDefault: 1.0E+40 (for minimization problems); -1.0E+40 (for maximization problems).",
      XPRS_MIPABSCUTOFF, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_MIPABSCUTOFF

#ifdef XPRS_MIPABSGAPNOTIFY
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_mipabsgapnotify XPRS_MIPABSGAPNOTIFY",
      "Branch and bound: if the gapnotify callback has been set using XPRSaddcbgapnotify, then this callback will be triggered during the tree search when the absolute gap reaches or passes the value you set of the MIPRELGAPNOTIFY control."
      "\n\nDefault: -1.0",
      XPRS_MIPABSGAPNOTIFY, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_MIPABSGAPNOTIFY

#ifdef XPRS_MIPABSGAPNOTIFYBOUND
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_mipabsgapnotifybound XPRS_MIPABSGAPNOTIFYBOUND",
      "Branch and bound: if the gapnotify callback has been set using XPRSaddcbgapnotify, then this callback will be triggered during the tree search when the best bound reaches or passes the value you set of the MIPRELGAPNOTIFYBOUND control."
      "\n\nDefault: 1.0E+20 (for minimization problems); -1.0E+20 (for maximization problems)",
      XPRS_MIPABSGAPNOTIFYBOUND, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_MIPABSGAPNOTIFYBOUND

#ifdef XPRS_MIPABSGAPNOTIFYOBJ
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_mipabsgapnotifyobj XPRS_MIPABSGAPNOTIFYOBJ",
      "Branch and bound: if the gapnotify callback has been set using XPRSaddcbgapnotify, then this callback will be triggered during the tree search when the best solution value reaches or passes the value you set of the MIPRELGAPNOTIFYOBJ control."
      "\n\nDefault: -1.0E+20 (for minimization problems); 1.0E+20 (for maximization problems)",
      XPRS_MIPABSGAPNOTIFYOBJ, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_MIPABSGAPNOTIFYOBJ

#ifdef XPRS_MIPABSSTOP
    MPD( AddSolverOption_MergeDuplicates("lim:xprs_mipabsstop XPRS_MIPABSSTOP",
      "Branch and Bound: The absolute tolerance determining whether the tree search will continue or not. It will terminate if    |MIPOBJVAL - BESTBOUND| ≤ MIPABSSTOP where MIPOBJVAL is the value of the best solution's objective function, and BESTBOUND is the current best solution bound. For example, to stop the tree search when a MIP solution has been found and the Optimizer can guarantee it is within 100 of the optimal solution, set MIPABSSTOP to 100."
      "\n\nDefault: 0.0",
      XPRS_MIPABSSTOP, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_MIPABSSTOP

#ifdef XPRS_MIPADDCUTOFF
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_mipaddcutoff XPRS_MIPADDCUTOFF",
      "Branch and Bound: The amount to add to the objective function of the best integer solution found to give the new CURRMIPCUTOFF. Once an integer solution has been found whose objective function is equal to or better than CURRMIPCUTOFF, improvements on this value may not be interesting unless they are better by at least a certain amount. If MIPADDCUTOFF is nonzero, it will be added to CURRMIPCUTOFF each time an integer solution is found which is better than this new value. This cuts off sections of the tree whose solutions would not represent substantial improvements in the objective function, saving processor time. The control MIPABSSTOP provides a similar function but works in a different way."
      "\n\nDefault: -1.0E-05",
      XPRS_MIPADDCUTOFF, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_MIPADDCUTOFF

#ifdef XPRS_MIPCOMPONENTS
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_mipcomponents XPRS_MIPCOMPONENTS",
      "       Determines whether disconnected components in a MIP should be solved as separate MIPs.       There can be significant performence benefits from solving disconnected components individual instead of being part of the main branch-and-bound search.     "
      "\n\n"
      "Values (default:        -1     ):\n"
      "\n- (-1)  Automatic - let the solver decide."
      "\n- (0)  Disable solving disconnected components separately."
      "\n- (1)  Solve disconnected components separately.",
      XPRS_MIPCOMPONENTS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MIPCOMPONENTS

#ifdef XPRS_MIPCONCURRENTNODES
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_mipconcurrentnodes XPRS_MIPCONCURRENTNODES",
      "     Sets the node limit for when a winning solve is selected when concurrent MIP solves are enabled. When multiple MIP solves are started,     they each run up to the MIPCONCURRENTNODES node limit and only one winning solve is selected for contuinuing the search with.   "
      "\n\n"
      "Values (default:      -1   ):\n"
      "\n- (-1)  Automatic - let the solver decide on a node limit."
      "\n- (>0)  Number of nodes each concurrent solve should complete before a winner is selected.",
      XPRS_MIPCONCURRENTNODES, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MIPCONCURRENTNODES

#ifdef XPRS_MIPCONCURRENTSOLVES
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_mipconcurrentsolves XPRS_MIPCONCURRENTSOLVES",
      "     Selects the number of concurrent solves to start for a MIP. Each solve will use a unique random seed for its random number generator, but will otherwise apply the same user controls.     The first concurrent solve to complete will have solved the MIP and all the concurrent solves will be terminated at this point.     Using concurrent solves can be advantageous when a MIP displays a high level of performance variability.   "
      "\n\n"
      "Values (default:      0   ):\n"
      "\n- (-1)  Enabled. The number of concurrent solves depends on MIPTHREADS."
      "\n- (0, 1)  Disabled"
      "\n- (n>1)  Enabled. The number of concurrent solves to start is given by n.",
      XPRS_MIPCONCURRENTSOLVES, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MIPCONCURRENTSOLVES

#ifdef XPRS_MIPDUALREDUCTIONS
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_mipdualreductions XPRS_MIPDUALREDUCTIONS",
      "Branch and Bound: Limits operations that can reduce the MIP solution space."
      "\n\n"
      "Values (default: 1):\n"
      "\n- (2)  Allow dual reductions on continuous variables only."
      "\n- (1)  Allow all dual reductions."
      "\n- (0)  Prevent all dual reductions.",
      XPRS_MIPDUALREDUCTIONS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MIPDUALREDUCTIONS

#ifdef XPRS_MIPFRACREDUCE
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_mipfracreduce XPRS_MIPFRACREDUCE",
      "Branch and Bound: Specifies how often the optimizer should run a heuristic to reduce the number of fractional integer variables in the node LP solutions."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatic."
      "\n- (0)  Disabled."
      "\n- (1)  Run before and after cutting on the root node."
      "\n- (2)  Run also during root cutting."
      "\n- (3)  Run also during the tree search.",
      XPRS_MIPFRACREDUCE, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MIPFRACREDUCE

#ifdef XPRS_MIPKAPPAFREQ
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_mipkappafreq XPRS_MIPKAPPAFREQ",
      "Branch and Bound: Specifies how frequently the basis condition number (also known as kappa) should be calculated during the branch-and-bound search."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  Do not calculate condition numbers."
      "\n- (1)  Calculate conditions numbers on every node, including after each round of root cutting."
      "\n- (n>1)  Calculate a condition number once per node of every n'th level of the branch-and-bound tree.",
      XPRS_MIPKAPPAFREQ, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MIPKAPPAFREQ

#ifdef XPRS_MIPLOG
    MPD( AddSolverOption_MergeDuplicates("log:xprs_miplog XPRS_MIPLOG",
      "MIP log print control."
      "\n\n"
      "Values (default: -100):\n"
      "\n- (-n)  Print out summary log at each nth node."
      "\n- (0)  No printout during MIP tree search."
      "\n- (1)  Only print out summary statement at the end."
      "\n- (2)  Print out detailed log at all solutions found."
      "\n- (3)  Print out detailed log at each node.",
      XPRS_MIPLOG, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MIPLOG

#ifdef XPRS_MIPPRESOLVE
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_mippresolve XPRS_MIPPRESOLVE",
      "Branch and Bound: Type of integer processing to be performed. If set to 0, no processing will be performed. This is a bit-vector control (see Section Bit-vector controls)."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (0)  Reduced cost fixing will be performed at each node. This can simplify the node before  it is solved, by deducing that certain variables' values can be fixed based on additional bounds imposed on  other variables at this node."
      "\n- (1)  Primal reductions will be performed at each node. Uses constraints of the node to  tighten the range of variables, often resulting in fixing their values. This greatly simplifies the  problem and may even determine optimality or infeasibility of the node before the simplex method commences."
      "\n- (2)  [Unused] This bit is no longer used to control probing. Refer to the integer control PREPROBING for setting probing level during presolve."
      "\n- (3)  If node preprocessing is allowed to change bounds on continuous columns."
      "\n- (4)  Dual reductions will be performed at each node."
      "\n- (5)  Allow global (non-bound) tightening of the problem during the tree search."
      "\n- (6)  The objective function will be used to find reductions at each node."
      "\n- (7)  [Unused] This bit is no longer used to control restarts. Refer to the integer control MIPRESTART for disabling tree restarts."
      "\n- (8)  Allow that symmetry is used to presolve the node problem.",
      XPRS_MIPPRESOLVE, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MIPPRESOLVE

#ifdef XPRS_MIPRAMPUP
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_miprampup XPRS_MIPRAMPUP",
      "Controls the strategy used by the parallel MIP solver during the ramp-up phase of a branch-and-bound tree search."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatically determined."
      "\n- (0)  No special treatment during the ramp-up phase. Always run with the maximal number of tasks."
      "\n- (1)  Limit the number of tasks until the initial dives have completed.",
      XPRS_MIPRAMPUP, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MIPRAMPUP

#ifdef XPRS_MIPREFINEITERLIMIT
    MPD( AddSolverOption_MergeDuplicates("sol:xprs_miprefineiterlimit XPRS_MIPREFINEITERLIMIT",
      "This defines an effort limit expressed as simplex iterations for the MIP solution refiner. The limit is per reoptimizations in the MIP refiner."
      "\n\nDefault: -1 — determined automatically.",
      XPRS_MIPREFINEITERLIMIT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MIPREFINEITERLIMIT

#ifdef XPRS_MIPRELCUTOFF
    MPD( AddSolverOption_MergeDuplicates("lim:xprs_miprelcutoff XPRS_MIPRELCUTOFF",
      "Branch and Bound: Percentage of the incumbent value to be added to the value of the objective function when an integer solution is found, to give the new value of CURRMIPCUTOFF. The effect is to cut off the search in parts of the tree whose best possible objective function would not be substantially better than the current solution. The control MIPRELSTOP provides a similar functionality but works in a different way."
      "\n\nDefault: 1.0E-04",
      XPRS_MIPRELCUTOFF, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_MIPRELCUTOFF

#ifdef XPRS_MIPRELGAPNOTIFY
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_miprelgapnotify XPRS_MIPRELGAPNOTIFY",
      "Branch and bound: if the gapnotify callback has been set using XPRSaddcbgapnotify, then this callback will be triggered during the branch and bound tree search when the relative gap reaches or passes the value you set of the MIPRELGAPNOTIFY control."
      "\n\nDefault: -1.0",
      XPRS_MIPRELGAPNOTIFY, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_MIPRELGAPNOTIFY

#ifdef XPRS_MIPRELSTOP
    MPD( AddSolverOption_MergeDuplicates("lim:xprs_miprelstop XPRS_MIPRELSTOP",
      "Branch and Bound: This determines when the branch and bound tree search will terminate. Branch and bound tree search will stop if:    |MIPOBJVAL - BESTBOUND| ≤ MIPRELSTOP x max(|BESTBOUND|,|MIPOBJVAL|) where MIPOBJVAL is the value of the best solution's objective function and BESTBOUND is the current best solution bound. For example, to stop the tree search when a MIP solution has been found and the Optimizer can guarantee it is within 5% of the optimal solution, set MIPRELSTOP to 0.05."
      "\n\nDefault: 0.0001",
      XPRS_MIPRELSTOP, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_MIPRELSTOP

#ifdef XPRS_MIPRESTART
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_miprestart XPRS_MIPRESTART",
      "Branch and Bound: controls strategy for in-tree restarts."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Determined automatically (XPRS_MIPRESTART_DEFAULT)."
      "\n- (0)  Disable in-tree restarts (XPRS_MIPRESTART_OFF)."
      "\n- (1)  Allow in-tree restarts at normal aggressiveness (XPRS_MIPRESTART_MODERATE)."
      "\n- (2)  Allow in-tree restarts at higher aggressiveness (more likely to trigger a restart) (XPRS_MIPRESTART_AGGRESSIVE).",
      XPRS_MIPRESTART, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MIPRESTART

#ifdef XPRS_MIPRESTARTFACTOR
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_miprestartfactor XPRS_MIPRESTARTFACTOR",
      "Branch and Bound: Fine tune initial conditions to trigger an in-tree restart. Use a value > 1 to increase the aggressiveness with which the Optimizer restarts. Use a value < 1 to relax the aggressiveness with which the Optimizer restarts. Note that this control does not affect the initial condition on the gap, which must be set separately. "
      "\n\nDefault: 1.0",
      XPRS_MIPRESTARTFACTOR, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_MIPRESTARTFACTOR

#ifdef XPRS_MIPRESTARTGAPTHRESHOLD
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_miprestartgapthreshold XPRS_MIPRESTARTGAPTHRESHOLD",
      "Branch and Bound: Initial gap threshold to delay in-tree restart. The restart is delayed initially if the gap, given as a fraction between 0 and 1, is below this threshold. The optimizer adjusts the threshold every time a restart is delayed. Note that there are other criteria that can delay or prevent a restart."
      "\n\nDefault: 0.02",
      XPRS_MIPRESTARTGAPTHRESHOLD, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_MIPRESTARTGAPTHRESHOLD

#ifdef XPRS_MIPTERMINATIONMETHOD
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_mipterminationmethod XPRS_MIPTERMINATIONMETHOD",
      " Branch and Bound: How a MIP solve should be stopped on early termination when there are still active tasks in the system. This can happen when, for example, a time or node limit is reached. "
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  Terminate tasks at the earliest opportunity. This can result in some unfinished node solves being discarded, although never integer solutions."
      "\n- (1)  Allow tasks to complete their current work but prevent new tasks from being started.",
      XPRS_MIPTERMINATIONMETHOD, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MIPTERMINATIONMETHOD

#ifdef XPRS_MIPTHREADS
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_mipthreads XPRS_MIPTHREADS",
      "If set to a positive integer it determines the number of threads implemented to run the parallel MIP code. If MIPTHREADS is set to the default value (-1), the THREADS control will determine the number of threads used. "
      "\n\nDefault: -1 (determined by the THREADS control)",
      XPRS_MIPTHREADS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MIPTHREADS

#ifdef XPRS_MIPTOL
    MPD( AddSolverOption_MergeDuplicates("tol:xprs_miptol XPRS_MIPTOL",
      "Branch and Bound: This is the tolerance within which a decision variable's value is considered to be integral."
      "\n\nDefault: 5.0E-06",
      XPRS_MIPTOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_MIPTOL

#ifdef XPRS_MIPTOLTARGET
    MPD( AddSolverOption_MergeDuplicates("sol:xprs_miptoltarget XPRS_MIPTOLTARGET",
      "Target MIPTOL value used by the automatic MIP solution refiner as defined by REFINEOPS. Negative and zero values are ignored."
      "\n\nDefault: 0.0",
      XPRS_MIPTOLTARGET, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_MIPTOLTARGET

#ifdef XPRS_MIQCPALG
    MPD( AddSolverOption_MergeDuplicates("qp:xprs_miqcpalg XPRS_MIQCPALG",
      "This control determines which algorithm is to be used to solve mixed integer quadratic constrained and mixed integer second order cone problems."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Determined automatically."
      "\n- (0)  Use the barrier algorithm in the branch and bound algorithm."
      "\n- (1)  Use outer approximations in the branch and bound algorithm.",
      XPRS_MIQCPALG, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MIQCPALG

#ifdef XPRS_MPS18COMPATIBLE
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_mps18compatible XPRS_MPS18COMPATIBLE",
      "Provides compatibility of MPS file output for older MPS readers."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (Bit 0)  Do not write objective sense (OBJSENSE section)."
      "\n- (Bit 1)  Fixed binaries are written as fixed only (unless used as a base variable for an indicator constraint).",
      XPRS_MPS18COMPATIBLE, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MPS18COMPATIBLE

#ifdef XPRS_MPSBOUNDNAME
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_mpsboundname XPRS_MPSBOUNDNAME",
      "When reading an MPS file, this control determines which entries from the BOUNDS section will be read. As with all string controls, this is of length 64 characters plus a null terminator, \0."
      "\n\nDefault: 64 blanks",
      XPRS_MPSBOUNDNAME) );
#endif  // ifdef XPRS_MPSBOUNDNAME

#ifdef XPRS_MPSECHO
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_mpsecho XPRS_MPSECHO",
      "Determines whether comments in MPS matrix files are to be printed out during matrix input."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  MPS comments are not to be echoed."
      "\n- (1)  MPS comments are to be echoed.",
      XPRS_MPSECHO, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MPSECHO

#ifdef XPRS_MPSFORMAT
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_mpsformat XPRS_MPSFORMAT",
      "Specifies the format of MPS files."
      "\n\n"
      "Values (default: 1):\n"
      "\n- (-1)  To determine the file type automatically."
      "\n- (0)  For fixed format."
      "\n- (1)  If MPS files are assumed to be in free format by input.",
      XPRS_MPSFORMAT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MPSFORMAT

#ifdef XPRS_MPSOBJNAME
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_mpsobjname XPRS_MPSOBJNAME",
      "When reading an MPS file, this control determines which neutral row will be read as the objective function. If this control is set when reading a multi-objective MPS file, only the named objective will be read; all other objectives will be ignored. As with all string controls, this is of length 64 characters plus a null terminator, \0. "
      "\n\nDefault: 64 blanks",
      XPRS_MPSOBJNAME) );
#endif  // ifdef XPRS_MPSOBJNAME

#ifdef XPRS_MPSRANGENAME
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_mpsrangename XPRS_MPSRANGENAME",
      "When reading an MPS file, this control determines which entries from the RANGES section will be read. As with all string controls, this is of length 64 characters plus a null terminator, \0."
      "\n\nDefault: 64 blanks",
      XPRS_MPSRANGENAME) );
#endif  // ifdef XPRS_MPSRANGENAME

#ifdef XPRS_MPSRHSNAME
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_mpsrhsname XPRS_MPSRHSNAME",
      "When reading an MPS file, this control determines which entries from the RHS section will be read. As with all string controls, this is of length 64 characters plus a null terminator, \0."
      "\n\nDefault: 64 blanks",
      XPRS_MPSRHSNAME) );
#endif  // ifdef XPRS_MPSRHSNAME

#ifdef XPRS_MULTIOBJLOG
    MPD( AddSolverOption_MergeDuplicates("obj:multi:xprs_multiobjlog XPRS_MULTIOBJLOG",
      "Log level for multi-objective optimization."
      "\n\n"
      "Values (default: 2):\n"
      "\n- (0)  No logging."
      "\n- (1)  Print a summary of each problem that is solved as part of the multi-objective optimization."
      "\n- (2)  In addition to summaries, print messages produced by each solve at the level determined by OUTPUTLOG.",
      XPRS_MULTIOBJLOG, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MULTIOBJLOG

#ifdef XPRS_MULTIOBJOPS
    MPD( AddSolverOption_MergeDuplicates("obj:multi:xprs_multiobjops XPRS_MULTIOBJOPS",
      "Modifies the behaviour of the optimizer when solving multi-objective problems."
      "\n\n"
      "Values (default: 7 (all bits are set, see Section Bit-vector controls for bit-vector controls)):\n"
      "\n- (0)  XPRS_MULTIOBJOPS_ENABLEDMulti-objective enabled. If this bit is not set, multi-objective problems will treated as single-objective problems, and only objective 0 will be optimized."
      "\n- (1)  XPRS_MULTIOBJOPS_PRESOLVEApply multi-objective modifications during presolve. If this bit is not set, the original problem will be modified when solving each subsequent objective, and these modifications will remain in the problem after the solve has completed."
      "\n- (2)  XPRS_MULTIOBJOPS_RCFIXINGReduced cost fixing. If this bit is set, optimality of earlier objectives will be preserved by fixing all non-basic variables with non-zero reduced costs to their bounds. If not set, optimality of earlier objectives will be preserved by adding constraints to the problem.",
      XPRS_MULTIOBJOPS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MULTIOBJOPS

#ifdef XPRS_MUTEXCALLBACKS
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_mutexcallbacks XPRS_MUTEXCALLBACKS",
      "Branch and Bound: This determines whether the callback routines are mutexed from within the optimizer."
      "\n\n"
      "Values (default: 1):\n"
      "\n- (0)  Callbacks are not mutexed."
      "\n- (1)  Callbacks are mutexed.",
      XPRS_MUTEXCALLBACKS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_MUTEXCALLBACKS

#ifdef XPRS_NETSTALLLIMIT
    MPD( AddSolverOption_MergeDuplicates("sim:xprs_netstalllimit XPRS_NETSTALLLIMIT",
      "Limit the number of degenerate pivots of the network simplex algorithm, before switching to either primal or dual simplex, depending on ALGAFTERNETWORK."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatically determined limit"
      "\n- (0)  No limit."
      "\n- (n>0)  Limit to n network simplex iterations.",
      XPRS_NETSTALLLIMIT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_NETSTALLLIMIT

#ifdef XPRS_NODEPROBINGEFFORT
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_nodeprobingeffort XPRS_NODEPROBINGEFFORT",
      "Adjusts the overall level of node probing."
      "\n\nDefault: 1.0",
      XPRS_NODEPROBINGEFFORT, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_NODEPROBINGEFFORT

#ifdef XPRS_NODESELECTION
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_nodeselection XPRS_NODESELECTION",
      "Branch and Bound: This determines which nodes will be considered for solution once the current node has been solved."
      "\n\n"
      "Values (default: Dependent on the matrix characteristics.):\n"
      "\n- (1)  Local first: Choose between descendant and sibling nodes if available; choose from all outstanding nodes otherwise."
      "\n- (2)  Best first: Choose from all outstanding nodes."
      "\n- (3)  Local depth first: Choose between descendant and sibling nodes if available; choose from the deepest nodes otherwise."
      "\n- (4)  Best first, then local first: Best first is used for the first BREADTHFIRST nodes, after which local first is used."
      "\n- (5)  Pure depth first: Choose from the deepest outstanding nodes.",
      XPRS_NODESELECTION, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_NODESELECTION

#ifdef XPRS_NUMERICALEMPHASIS
    MPD( AddSolverOption_MergeDuplicates("num:xprs_numericalemphasis XPRS_NUMERICALEMPHASIS",
      "How much emphasis to place on numerical stability instead of solve speed. "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatic. The emphasis might be influenced by the setting of other controls."
      "\n- (0)  Emphasize speed."
      "\n- (1)  Mild emphasis on numerical stability."
      "\n- (2)  Medium emphasis on numerical stability."
      "\n- (3)  Strong emphasis on numerical stability.",
      XPRS_NUMERICALEMPHASIS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_NUMERICALEMPHASIS

#ifdef XPRS_OBJSCALEFACTOR
    MPD( AddSolverOption_MergeDuplicates("num:xprs_objscalefactor XPRS_OBJSCALEFACTOR",
      "Custom objective scaling factor, expressed as a power of 2. When set, it overwrites the automatic objective scaling factor. A value of 0 means no objective scaling. This control is applied for the full solve, and is independent of any extra scaling that may occur specifically for the barrier or simplex solvers. As it is a power of 2, to scale by 16, set the value of the control to 4."
      "\n\nDefault: 0",
      XPRS_OBJSCALEFACTOR, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_OBJSCALEFACTOR

#ifdef XPRS_OPTIMALITYTOL
    MPD( AddSolverOption_MergeDuplicates("tol:xprs_optimalitytol XPRS_OPTIMALITYTOL",
      "Simplex: This is the zero tolerance for reduced costs. On each iteration, the simplex method searches for a variable to enter the basis which has a negative reduced cost. The candidates are only those variables which have reduced costs less than the negative value of OPTIMALITYTOL."
      "\n\nDefault: 1.0E-06",
      XPRS_OPTIMALITYTOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_OPTIMALITYTOL

#ifdef XPRS_OPTIMALITYTOLTARGET
    MPD( AddSolverOption_MergeDuplicates("sol:xprs_optimalitytoltarget XPRS_OPTIMALITYTOLTARGET",
      "This specifies the target optimality tolerance for the solution refiner."
      "\n\nDefault: 0 — use the value specified by OPTIMALITYTOL.",
      XPRS_OPTIMALITYTOLTARGET, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_OPTIMALITYTOLTARGET

#ifdef XPRS_OUTPUTCONTROLS
    MPD( AddSolverOption_MergeDuplicates("log:xprs_outputcontrols XPRS_OUTPUTCONTROLS",
      "This control toggles the printing of all control settings at the beginning of the search. This includes the printing of controls that have been explicitly assigned to their default value. All unset controls are omitted as they keep their default value."
      "\n\n"
      "Values (default: 1):\n"
      "\n- (0)  Turn off printing of user-specified control settings."
      "\n- (1)  Print controls.",
      XPRS_OUTPUTCONTROLS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_OUTPUTCONTROLS

#ifdef XPRS_OUTPUTLOG
    MPD( AddSolverOption_MergeDuplicates("log:xprs_outputlog XPRS_OUTPUTLOG",
      "This controls the level of output produced by the Optimizer during optimization. In the Console Optimizer, OUTPUTLOG controls which messages are sent to the screen (stdout). When using the Optimizer library, no output is sent to the screen. If the user wishes output to be displayed, they must define a callback function and print messages to the screen themselves. In this case, OUTPUTLOG controls which messages are sent to the user output callback."
      "\n\n"
      "Values (default: 1):\n"
      "\n- (0)  Turn all output off. Use XPRS_OUTPUTLOG_NO_OUTPUT from xprs.h."
      "\n- (1)  Print all messages. Use XPRS_OUTPUTLOG_FULL_OUTPUT from xprs.h."
      "\n- (3)  Print error and warning messages. Use XPRS_OUTPUTLOG_ERRORS_AND_WARNINGS from xprs.h."
      "\n- (4)  Print error messages only. Use XPRS_OUTPUTLOG_ERRORS from xprs.h.",
      XPRS_OUTPUTLOG, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_OUTPUTLOG

#ifdef XPRS_OUTPUTMASK
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_outputmask XPRS_OUTPUTMASK",
      "Mask to restrict the row and column names written to file. As with all string controls, this is of length 64 characters plus a null terminator, \0."
      "\n\nDefault: 64 '?'s",
      XPRS_OUTPUTMASK) );
#endif  // ifdef XPRS_OUTPUTMASK

#ifdef XPRS_OUTPUTTOL
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_outputtol XPRS_OUTPUTTOL",
      "Zero tolerance on print values."
      "\n\nDefault: 1.0E-05",
      XPRS_OUTPUTTOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_OUTPUTTOL

#ifdef XPRS_PENALTY
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_penalty XPRS_PENALTY",
      "Minimum absolute penalty variable coefficient. BIGM and PENALTY are set by the input routine (XPRSreadprob (READPROB)) but may be reset by the user prior to  XPRSlpoptimize (LPOPTIMIZE). "
      "\n\nDefault: Dependent on the matrix characteristics.",
      XPRS_PENALTY, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_PENALTY

#ifdef XPRS_PIVOTTOL
    MPD( AddSolverOption_MergeDuplicates("sim:xprs_pivottol XPRS_PIVOTTOL",
      "Simplex: The zero tolerance for matrix elements. On each iteration, the simplex method seeks a nonzero matrix element to pivot on. Any element with absolute value less than PIVOTTOL is treated as zero for this purpose."
      "\n\nDefault: 1.0E-09",
      XPRS_PIVOTTOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_PIVOTTOL

#ifdef XPRS_PPFACTOR
    MPD( AddSolverOption_MergeDuplicates("sim:xprs_ppfactor XPRS_PPFACTOR",
      "The partial pricing candidate list sizing parameter."
      "\n\nDefault: 1.0",
      XPRS_PPFACTOR, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_PPFACTOR

#ifdef XPRS_PREANALYTICCENTER
    MPD( AddSolverOption_MergeDuplicates("bar:xprs_preanalyticcenter XPRS_PREANALYTICCENTER",
      " Determines if analytic centers should be computed and used for variable fixing and the generation of alternative reduced costs (-1: Auto 0: Off, 1: Fixing, 2: Redcost, 3: Both) "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatic."
      "\n- (0)  Disable analytic center presolving."
      "\n- (1)  Use analytic center for variable fixing only."
      "\n- (2)  Use analytic center for reduced cost computation only."
      "\n- (3)  Use analytic centers for both, variable fixing and reduced cost computation.",
      XPRS_PREANALYTICCENTER, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PREANALYTICCENTER

#ifdef XPRS_PREBASISRED
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_prebasisred XPRS_PREBASISRED",
      "Determines if a lattice basis reduction algorithm should be attempted as part of presolve"
      "\n\n"
      "Values (default: 0):\n"
      "\n- (-1)  Automatic."
      "\n- (0)  Disable basis reduction."
      "\n- (1)  Enable basis reduction.",
      XPRS_PREBASISRED, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PREBASISRED

#ifdef XPRS_PREBNDREDCONE
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_prebndredcone XPRS_PREBNDREDCONE",
      "Determines if second order cone constraints should be used for inferring bound reductions on variables when solving a MIP."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatic."
      "\n- (0)  Disable bound reductions from second order cone constraints."
      "\n- (1)  Enable bound reductions from second order cone constraints.",
      XPRS_PREBNDREDCONE, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PREBNDREDCONE

#ifdef XPRS_PREBNDREDQUAD
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_prebndredquad XPRS_PREBNDREDQUAD",
      "Determines if convex quadratic constraints should be used for inferring bound reductions on variables when solving a MIP."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatic."
      "\n- (0)  Disable bound reductions from quadratic constraints."
      "\n- (1)  Enable bound reductions from quadratic constraints.",
      XPRS_PREBNDREDQUAD, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PREBNDREDQUAD

#ifdef XPRS_PRECLIQUESTRATEGY
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_precliquestrategy XPRS_PRECLIQUESTRATEGY",
      "Determines how much effort to spend on clique covers in presolve."
      "\n\nDefault: -1",
      XPRS_PRECLIQUESTRATEGY, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PRECLIQUESTRATEGY

#ifdef XPRS_PRECOEFELIM
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_precoefelim XPRS_PRECOEFELIM",
      "Presolve: Specifies whether the optimizer should attempt to recombine constraints in order to reduce the number of non zero coefficients when presolving a mixed integer problem. "
      "\n\n"
      "Values (default: 2):\n"
      "\n- (0)  Disabled."
      "\n- (1)  Remove as many coefficients as possible."
      "\n- (2)  Cautious eliminations. Will not perform a reduction if it might destroy problem structure useful to e.g. heuristics or cutting.",
      XPRS_PRECOEFELIM, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PRECOEFELIM

#ifdef XPRS_PRECOMPONENTS
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_precomponents XPRS_PRECOMPONENTS",
      "Presolve: determines whether small independent components should be detected and solved as individual subproblems during root node processing. "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatically determined."
      "\n- (0)  Disable detection of independent components."
      "\n- (1)  Enable detection of independent components.",
      XPRS_PRECOMPONENTS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PRECOMPONENTS

#ifdef XPRS_PRECOMPONENTSEFFORT
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_precomponentseffort XPRS_PRECOMPONENTSEFFORT",
      "Presolve: adjusts the overall effort for the independent component presolver. This control affects working limits for the subproblem solving as well as thresholds when it is called. Increase to put more emphasis on component presolving. "
      "\n\nDefault: 1.0",
      XPRS_PRECOMPONENTSEFFORT, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_PRECOMPONENTSEFFORT

#ifdef XPRS_PRECONEDECOMP
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_preconedecomp XPRS_PRECONEDECOMP",
      "Presolve: decompose regular and rotated cones with more than two elements and apply Outer Approximation on the resulting components. "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatically determined."
      "\n- (0)  Disable cone decomposition."
      "\n- (1)  Enable cone decomposition by replacing large  cones with small ones in the presolved problem."
      "\n- (2)  Similar to 1, plus decomposition is enabled even  if the cone variable is fixed."
      "\n- (3)  Cones are decomposed within the Outer  Approximation domain only, i.e., the problem maintains the original  cones.",
      XPRS_PRECONEDECOMP, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PRECONEDECOMP

#ifdef XPRS_PRECONFIGURATION
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_preconfiguration XPRS_PRECONFIGURATION",
      "MIP Presolve: determines whether binary rows with only few repeating coefficients should be reformulated. The reformulation enumerates the extremal feasible configurations of a row and introduces new columns and rows to model the choice between these extremal configurations. This presolve operation can be disabled as part of the (advanced) IP reductions PRESOLVEOPS. "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatically determined."
      "\n- (0)  Disable configuration presolving.",
      XPRS_PRECONFIGURATION, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PRECONFIGURATION

#ifdef XPRS_PRECONVERTOBJTOCONS
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_preconvertobjtocons XPRS_PRECONVERTOBJTOCONS",
      "Presolve: convert a linear or quadratic objective function into an objective transfer constraint "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatically determined."
      "\n- (0)  Disable reformulation."
      "\n- (1)  Move only the quadratic part of the objective into a constraint."
      "\n- (2)  Move both the linear and quadratic parts of the objective into a constraint.",
      XPRS_PRECONVERTOBJTOCONS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PRECONVERTOBJTOCONS

#ifdef XPRS_PRECONVERTSEPARABLE
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_preconvertseparable XPRS_PRECONVERTSEPARABLE",
      "Presolve: reformulate problems with a non-diagonal quadratic objective and/or constraints as diagonal quadratic or second-order conic constraints. "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatically determined."
      "\n- (0)  Disable reformulation."
      "\n- (1)  Enable reformulation to diagonal quadratic constraints."
      "\n- (2)  Enable reformulation to diagonal quadratic constraints and reduction to second-order cones.",
      XPRS_PRECONVERTSEPARABLE, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PRECONVERTSEPARABLE

#ifdef XPRS_PREDOMCOL
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_predomcol XPRS_PREDOMCOL",
      "Presolve: Determines the level of dominated column removal reductions to perform when presolving a mixed integer problem. Only binary columns will be checked. "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatically determined."
      "\n- (0)  Disabled."
      "\n- (1)      Cautious strategy, limited effort looking for special structure.  "
      "\n- (2)      Same as 2 but checking all candidates.  "
      "\n- (3)      Includes 1 and 2 but also looks for more generic column domination.  ",
      XPRS_PREDOMCOL, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PREDOMCOL

#ifdef XPRS_PREDOMROW
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_predomrow XPRS_PREDOMROW",
      "Presolve: Determines the level of dominated row removal reductions to perform when presolving a problem. "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatically determined."
      "\n- (0)  Disabled."
      "\n- (1)  Cautious strategy."
      "\n- (2)  Medium strategy."
      "\n- (3)  Aggressive strategy. All candidate row combinations will be considered.",
      XPRS_PREDOMROW, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PREDOMROW

#ifdef XPRS_PREDUPROW
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_preduprow XPRS_PREDUPROW",
      "Presolve: Determines the type of duplicate rows to look for and eliminate when presolving a problem."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatically determined."
      "\n- (0)  Do not eliminate duplicate rows."
      "\n- (1)  Eliminate only rows that are identical in all variables."
      "\n- (2)  Same as option 1 plus eliminate duplicate rows with simple penalty variable expressions. (MIP only)."
      "\n- (3)  Same as option 2 plus eliminate duplicate rows with more complex penalty variable expressions. (MIP only).",
      XPRS_PREDUPROW, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PREDUPROW

#ifdef XPRS_PREELIMQUAD
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_preelimquad XPRS_PREELIMQUAD",
      "Presolve: Allows for elimination of quadratic variables via doubleton rows."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatically determined."
      "\n- (0)  Do not eliminate duplicate rows."
      "\n- (1)  Eliminate at least one quadratic variable for each doubleton row.",
      XPRS_PREELIMQUAD, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PREELIMQUAD

#ifdef XPRS_PREFOLDING
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_prefolding XPRS_PREFOLDING",
      "Presolve: Determines if a folding procedure should be used to aggregate continuous columns in an equitable partition. "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatically determined."
      "\n- (0)  Disabled."
      "\n- (1)  Enabled.",
      XPRS_PREFOLDING, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PREFOLDING

#ifdef XPRS_PREIMPLICATIONS
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_preimplications XPRS_PREIMPLICATIONS",
      "Presolve: Determines whether to use implication structures to remove redundant rows. If implication sequences are detected, they might also be used in probing."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatically determined."
      "\n- (0)  Do not use implications for sparsification."
      "\n- (1)  Use implications to remove reduandant rows.",
      XPRS_PREIMPLICATIONS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PREIMPLICATIONS

#ifdef XPRS_PRELINDEP
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_prelindep XPRS_PRELINDEP",
      "Presolve: Determines whether to check for and remove linearly dependent equality constraints when presolving a problem."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatically determined."
      "\n- (0)  Do not check for linearly dependent equality constraints."
      "\n- (1)  Check for and remove linearly dependent equality constraints.",
      XPRS_PRELINDEP, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PRELINDEP

#ifdef XPRS_PREOBJCUTDETECT
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_preobjcutdetect XPRS_PREOBJCUTDETECT",
      "Presolve: Determines whether to check for constraints that are parallel or near parallel to a linear objective function, and which can safely be removed. This reduction applies to MIPs only."
      "\n\n"
      "Values (default: 1):\n"
      "\n- (0)  Disable check and reductions."
      "\n- (1)  Enable check and reductions.",
      XPRS_PREOBJCUTDETECT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PREOBJCUTDETECT

#ifdef XPRS_PREPERMUTE
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_prepermute XPRS_PREPERMUTE",
      "This bit-vector control (see Section Bit-vector controls) specifies whether to randomly permute rows, columns and MIP entities when starting the presolve. With the default value 0, no permutation will take place."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  Permute rows."
      "\n- (1)  Permute columns."
      "\n- (2)  Permute MIP entities. This bit only affects MIP problems.",
      XPRS_PREPERMUTE, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PREPERMUTE

#ifdef XPRS_PREPERMUTESEED
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_prepermuteseed XPRS_PREPERMUTESEED",
      "This control sets the seed for the pseudo-random number generator for permuting the problem when starting the presolve. This control only has effects when PREPERMUTE is enabled."
      "\n\nDefault: 1",
      XPRS_PREPERMUTESEED, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PREPERMUTESEED

#ifdef XPRS_PREPROBING
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_preprobing XPRS_PREPROBING",
      "Presolve: Amount of probing to perform on binary variables during presolve. This is done by fixing a binary to each of its values in turn and analyzing the implications. "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Let the optimizer decide on the amount of probing."
      "\n- (0)  Disabled."
      "\n- (+1)  Light probing — only few implications will be examined."
      "\n- (+2)  Full probing — all implications for all binaries will be examined."
      "\n- (+3)  Full probing and repeat as long as the problem is significantly reduced.",
      XPRS_PREPROBING, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PREPROBING

#ifdef XPRS_PREPROTECTDUAL
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_preprotectdual XPRS_PREPROTECTDUAL",
      "Presolve: specifies whether the presolver should protect a given dual solution by maintaining the same level of dual feasibility. Enabling this control often results in a worse presolved model. This control only expected to be optionally enabled before calling XPRScrossoverlpsol. "
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  Disabled."
      "\n- (1)  Enabled. Protect the dual solution during presolve.",
      XPRS_PREPROTECTDUAL, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PREPROTECTDUAL

#ifdef XPRS_PREROOTEFFORT
    MPD( AddSolverOption_MergeDuplicates("heur:xprs_prerooteffort XPRS_PREROOTEFFORT",
      "Dial for the work spent during the Pre-root parallel heuristic phase. A positive value sets a suitable work limit that is dependent on problem-characteristics. Changing the value up/or down dials the work spent in this phase up or down. This control also enables/disables Pre-root parallel heuristics. "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-2)  Enable Pre-root parallel heuristics without a specific work limit for this phase. The phase will terminate if a different limit is hit, or if it runs out of heuristic work to do."
      "\n- (-1)  Enablement of Pre-root parallel heuristics is subject to HEUREMPHASIS."
      "\n- (0)  Disable Pre-root parallel heuristics."
      "\n- (x>0)  Enable Pre-root parallel heuristics with a work limit dependent on problem characteristics, using x as a factor to dial this work limit up or down.",
      XPRS_PREROOTEFFORT, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_PREROOTEFFORT

#ifdef XPRS_PREROOTTHREADS
    MPD( AddSolverOption_MergeDuplicates("heur:xprs_prerootthreads XPRS_PREROOTTHREADS",
      "Specifies an explicit number of threads that should be used for the Pre-root parallel heuristic phase. By default, this phase will use all threads available to the solver (as governed by the control THREADS). "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Use all available threads."
      "\n- (0)  Disable pre-root parallel heuristics."
      "\n- (n>0)  Use the specified number of threads, superseding the value of THREADS",
      XPRS_PREROOTTHREADS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PREROOTTHREADS

#ifdef XPRS_PREROOTWORKLIMIT
    MPD( AddSolverOption_MergeDuplicates("heur:xprs_prerootworklimit XPRS_PREROOTWORKLIMIT",
      "Set an explicit work limit in work units for the Pre-root parallel heuristic phase. Any positive value also enables this phase and runs it until the PREROOTWORKLIMIT units of work are hit. "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  No explicit work limit for the Pre-root parallel heuristic phase. If enabled, the work limit for this phase is controlled via PREROOTEFFORT."
      "\n- (0)  Disable Pre-root parallel heuristics."
      "\n- (x>0)  Enable Pre-root parallel heuristics with an explicit work limit of x work units. If set, this work limit has precedence over any work limit set by PREROOTEFFORT.",
      XPRS_PREROOTWORKLIMIT, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_PREROOTWORKLIMIT

#ifdef XPRS_PRESOLVE
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_presolve XPRS_PRESOLVE",
      "This control determines whether presolving should be performed prior to starting the main algorithm. Presolve attempts to simplify the problem by detecting and removing redundant constraints, tightening variable bounds, etc. In some cases, infeasibility may even be determined at this stage, or the optimal solution found."
      "\n\n"
      "Values (default: 1):\n"
      "\n- (-1)  Presolve applied, but a problem will not be declared infeasible if primal infeasibilities are detected. The problem will be solved by the LP optimization algorithm, returning an infeasible solution, which can sometimes be helpful."
      "\n- (0)  Presolve not applied."
      "\n- (1)  Presolve applied."
      "\n- (2)  Presolve applied, but redundant bounds are not removed. This can sometimes increase the efficiency of the barrier algorithm."
      "\n- (3)  Presolve is applied, and bounds detected to be redundant are always removed.",
      XPRS_PRESOLVE, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PRESOLVE

#ifdef XPRS_PRESOLVEMAXGROW
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_presolvemaxgrow XPRS_PRESOLVEMAXGROW",
      "Limit on how much the number of non-zero coefficients is allowed to grow during presolve, specified as a ratio of the number of non-zero coefficients in the original problem."
      "\n\nDefault: 0.1",
      XPRS_PRESOLVEMAXGROW, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_PRESOLVEMAXGROW

#ifdef XPRS_PRESOLVEOPS
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_presolveops XPRS_PRESOLVEOPS",
      "This bit-vector control (see Section Bit-vector controls) specifies the operations which are performed during the presolve."
      "\n\n"
      "Values (default: 511 (bits 0 — 8 incl. are set)):\n"
      "\n- (0)  Singleton column removal."
      "\n- (1)  Singleton row removal."
      "\n- (2)  Forcing row removal."
      "\n- (3)  Dual reductions."
      "\n- (4)  Redundant row removal."
      "\n- (5)  Duplicate column removal."
      "\n- (6)  Duplicate row removal."
      "\n- (7)  Strong dual reductions."
      "\n- (8)  Variable eliminations."
      "\n- (9)  No IP reductions."
      "\n- (10)  No domain changes for MIP entities (e.g., semi-continuous detection or shifting integers)."
      "\n- (11)  No advanced IP reductions."
      "\n- (12)  No eliminations on integers."
      "\n- (13)  No reductions based on solution enumeration."
      "\n- (14)  Linearly dependant row removal."
      "\n- (15)  No integer variable and SOS detection."
      "\n- (16)  No implied bounds."
      "\n- (17)  No clique presolve."
      "\n- (18)  No mod2 presolve.",
      XPRS_PRESOLVEOPS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PRESOLVEOPS

#ifdef XPRS_PRESOLVEPASSES
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_presolvepasses XPRS_PRESOLVEPASSES",
      "Number of reduction rounds to be performed in presolve"
      "\n\nDefault: 1",
      XPRS_PRESOLVEPASSES, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PRESOLVEPASSES

#ifdef XPRS_PRESORT
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_presort XPRS_PRESORT",
      "This bit-vector control (see Section Bit-vector controls) specifies whether to sort rows, columns and MIP entities by their names when starting the presolve. With the default value 0, no sorting will take place."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  Sort rows."
      "\n- (1)  Sort columns."
      "\n- (2)  Sort MIP entities. This bit only affects MIP problems.",
      XPRS_PRESORT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PRESORT

#ifdef XPRS_PRICINGALG
    MPD( AddSolverOption_MergeDuplicates("lp:xprs_pricingalg XPRS_PRICINGALG",
      "Simplex: This determines the primal simplex pricing method. It is used to select which variable enters the basis on each iteration. In general Devex pricing requires more time on each iteration, but may reduce the total number of iterations, whereas partial pricing saves time on each iteration, but may result in more iterations."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (-1)  Partial pricing."
      "\n- (0)  Determined automatically."
      "\n- (1)  Devex pricing."
      "\n- (2)  Steepest edge."
      "\n- (3)  Steepest edge with unit initial weights.",
      XPRS_PRICINGALG, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PRICINGALG

#ifdef XPRS_PRIMALOPS
    MPD( AddSolverOption_MergeDuplicates("lp:xprs_primalops XPRS_PRIMALOPS",
      "Primal simplex: allows fine tuning the variable selection in the primal simplex solver."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (0)  Use aggressive dj scaling."
      "\n- (1)  Conventional dj scaling."
      "\n- (2)  Use reluctant switching back to partial pricing."
      "\n- (3)  Use dynamic switching between cheap and expensive pricing strategies."
      "\n- (4)  Keep solving even after potential cycling is detected.",
      XPRS_PRIMALOPS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PRIMALOPS

#ifdef XPRS_PRIMALPERTURB
    MPD( AddSolverOption_MergeDuplicates("lp:xprs_primalperturb XPRS_PRIMALPERTURB",
      "The factor by which the problem will be perturbed prior to optimization by primal simplex. A value of 0.0 results in no perturbation prior to optimization.  Note the interconnection to the AUTOPERTURB control. If AUTOPERTURB is set to 1, the decision whether to perturb or not is left to the Optimizer. When the problem is automatically perturbed in primal simplex, however, the value of PRIMALPERTURB will be used for perturbation. "
      "\n\nDefault: -1 — determined automatically. ",
      XPRS_PRIMALPERTURB, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_PRIMALPERTURB

#ifdef XPRS_PRIMALUNSHIFT
    MPD( AddSolverOption_MergeDuplicates("lp:xprs_primalunshift XPRS_PRIMALUNSHIFT",
      "Determines whether primal is allowed to call dual to unshift."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  Allow the dual algorithm to be used to unshift."
      "\n- (1)  Don't allow the dual algorithm to be used to unshift.",
      XPRS_PRIMALUNSHIFT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PRIMALUNSHIFT

#ifdef XPRS_PSEUDOCOST
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_pseudocost XPRS_PSEUDOCOST",
      "Branch and Bound: The default pseudo cost used in estimation of the degradation associated with an unexplored node in the tree search.  A pseudo cost is associated with each integer decision variable and is an estimate of the amount by which the objective function will be worse if that variable is forced to an integral value. "
      "\n\nDefault: 0.01",
      XPRS_PSEUDOCOST, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_PSEUDOCOST

#ifdef XPRS_PWLDUALREDUCTIONS
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_pwldualreductions XPRS_PWLDUALREDUCTIONS",
      "This parameter specifies whether dual reductions should be applied to reduce the number of columns, rows and SOS-constraints added when transforming piecewise linear objectives and constraints to MIP structs. "
      "\n\n"
      "Values (default: 1):\n"
      "\n- (0)  Disabled. No dual reductions, add all columns, rows and SOS-constraints."
      "\n- (1)  Enabled. Only add neccessary columns, rows and sets, drop those implied by the objective sense.",
      XPRS_PWLDUALREDUCTIONS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PWLDUALREDUCTIONS

#ifdef XPRS_PWLNONCONVEXTRANSFORMATION
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_pwlnonconvextransformation XPRS_PWLNONCONVEXTRANSFORMATION",
      "This control specifies the reformulation method for piecewise linear constraints  at the beginning of the search.   Note that the chosen formulation will only be used if MIP entities  are necessary but not if presolve detected that a convex reformulation is possible. Furthermore, the  binary formulation will only be applied to piecewise linear constraints with bounded input variable,  otherwise the SOS2-formulation will be used. "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatic."
      "\n- (0)  Use a formulation based on SOS2-constraints."
      "\n- (1)  Use a formulation based on binary variables.",
      XPRS_PWLNONCONVEXTRANSFORMATION, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_PWLNONCONVEXTRANSFORMATION

#ifdef XPRS_QCCUTS
    MPD( AddSolverOption_MergeDuplicates("cut:xprs_qccuts XPRS_QCCUTS",
      "Branch and Bound: Limit on the number of rounds of  outer approximation cuts generated for the root node, when solving a mixed integer quadratic constrained or mixed integer second order conic problem with outer approximation. "
      "\n\nDefault: -1 — determined automatically.",
      XPRS_QCCUTS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_QCCUTS

#ifdef XPRS_QCROOTALG
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_qcrootalg XPRS_QCROOTALG",
      "This control determines which algorithm is to be used to solve the root of a mixed integer quadratic constrained or mixed integer second order cone problem, when outer approximation is used. "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Determined automatically."
      "\n- (0)  Use the barrier algorithm."
      "\n- (1)  Use the dual simplex on a relaxation of the problem constructed using outer approximation.",
      XPRS_QCROOTALG, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_QCROOTALG

#ifdef XPRS_QSIMPLEXOPS
    MPD( AddSolverOption_MergeDuplicates("qp:xprs_qsimplexops XPRS_QSIMPLEXOPS",
      "Controls the behavior of the quadratic simplex solvers via a bit-vector (see Section Bit-vector controls)."
      "\n\n"
      "Values (default: 0 ):\n"
      "\n- (0)  Force traditional primal first phase."
      "\n- (1)  Force BigM primal first phase."
      "\n- (2)  Force traditional dual first phase. "
      "\n- (3)  Force BigM dual first phase."
      "\n- (4)  Always use artificial bounds in dual."
      "\n- (5)  Use original problem basis only when warmstarting the KKT."
      "\n- (6)  Skip the primal bound flips for ranged primals (might cause more trouble than good if the bounds are very large)."
      "\n- (7)  Also do the single pivot crash."
      "\n- (8)  Do not apply aggressive perturbation in dual."
      "\n- (9)  Applies standard scaling to the KKT system."
      "\n- (10)  Do not fall back to using Barrier in case of numerical difficulties with quadratic simplex during a MIP solve."
      "\n- (11)  Use primal simplex to solve the phase 1 feasibility problem before applying quadratic primal simplex."
      "\n- (12)  Use dual simplex to solve the phase 1 feasibility problem before applying quadratic primal simplex."
      "\n- (13)  Use barrier algorithm to solve the phase 1 feasibility problem before applying quadratic primal simplex."
      "\n- (14)  Use partial pricing."
      "\n- (15)  Use full pricing."
      "\n- (16)  Perform cleanup if a superbasic solution is provided for warm-start.",
      XPRS_QSIMPLEXOPS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_QSIMPLEXOPS

#ifdef XPRS_QUADRATICUNSHIFT
    MPD( AddSolverOption_MergeDuplicates("qp:xprs_quadraticunshift XPRS_QUADRATICUNSHIFT",
      "Determines whether an extra solution purification step is called after a solution found by the quadratic simplex (either primal or dual)."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Determined automatically."
      "\n- (0)  No purification step."
      "\n- (1)  Always do the purification step.",
      XPRS_QUADRATICUNSHIFT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_QUADRATICUNSHIFT

#ifdef XPRS_RANDOMSEED
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_randomseed XPRS_RANDOMSEED",
      "Sets the initial seed to use for the pseudo-random number generator in the Optimizer. The sequence of random numbers is always reset using the seed when starting a new optimization run. "
      "\n\nDefault: 1",
      XPRS_RANDOMSEED, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_RANDOMSEED

#ifdef XPRS_REFACTOR
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_refactor XPRS_REFACTOR",
      "Indicates whether the optimization should restart using the current representation of the factorization in memory."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatic."
      "\n- (0)  Do not refactor on reoptimizing."
      "\n- (1)  Refactor on reoptimizing.",
      XPRS_REFACTOR, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_REFACTOR

#ifdef XPRS_REFINEOPS
    MPD( AddSolverOption_MergeDuplicates("sol:xprs_refineops XPRS_REFINEOPS",
      "This specifies when the solution refiner should be executed to reduce solution infeasibilities. The refiner will attempt to satisfy the target tolerances for all original linear constraints before presolve or scaling has been applied."
      "\n\n"
      "Values (default: 19 (bits 0, 1 and 4 are set)):\n"
      "\n- (0)  Run the solution refiner on an optimal solution of a continuous problem."
      "\n- (1)  Run the solution refiner when a new solution is found during a tree search. The refiner will be applied to the presolved solution before any post-solve operations are applied."
      "\n- (3)  Run the solution refiner on each node of the MIP search. "
      "\n- (4)  Run the solution refiner on an optimal solution before postsolve on a continuous problem. "
      "\n- (5)  Apply the iterative refiner to refine the solution. "
      "\n- (6)  Use higher precision in the iterative refinement. "
      "\n- (7)  If set, the iterative refiner will use the primal simplex algorithm. "
      "\n- (8)  If set, the iterative refiner will use the dual simplex algorithm. "
      "\n- (9)  Refine MIP solutions such that rounding them keeps the problem feasible when reoptimized. "
      "\n- (10)  Attempt to refine MIP solutions such that rounding them keeps the problem feasible when reoptimized, but accept integers solutions even if refinement fails. ",
      XPRS_REFINEOPS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_REFINEOPS

#ifdef XPRS_RELAXTREEMEMORYLIMIT
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_relaxtreememorylimit XPRS_RELAXTREEMEMORYLIMIT",
      "When the memory used by the branch and bound search tree exceeds the target specified by the TREEMEMORYLIMIT control, the optimizer will try to reduce this by writing nodes to the tree file.  In rare cases, usually where the solve has many millions of very small nodes, the tree structural data (which cannot be written to the tree file) will grow large enough to approach or exceed the tree's memory target.  When this happens, optimizer performance can degrade greatly as the solver makes heavy use of the tree file in preference to memory.  To prevent this, the solver will automatically relax the tree memory limit when it detects this case; the RELAXTREEMEMORYLIMIT control specifies the proportion of the previous memory limit by which to relax it.  Set RELAXTREEMEMORYLIMIT to 0.0 to force the Xpress Optimizer to never relax the tree memory limit in this way.  "
      "\n\nDefault: 0.1",
      XPRS_RELAXTREEMEMORYLIMIT, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_RELAXTREEMEMORYLIMIT

#ifdef XPRS_RELPIVOTTOL
    MPD( AddSolverOption_MergeDuplicates("sim:xprs_relpivottol XPRS_RELPIVOTTOL",
      "Simplex: At each iteration a pivot element is chosen within a given column of the matrix. The relative pivot tolerance, RELPIVOTTOL, is the size of the element chosen relative to the largest possible pivot element in the same column."
      "\n\nDefault: 1.0E-06",
      XPRS_RELPIVOTTOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_RELPIVOTTOL

#ifdef XPRS_REPAIRINDEFINITEQ
    MPD( AddSolverOption_MergeDuplicates("qp:xprs_repairindefiniteq XPRS_REPAIRINDEFINITEQ",
      "Controls if the optimizer should make indefinite quadratic matrices positive definite when it is possible."
      "\n\n"
      "Values (default: 1):\n"
      "\n- (0)  Repair if possible."
      "\n- (1)  Do not repair.",
      XPRS_REPAIRINDEFINITEQ, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_REPAIRINDEFINITEQ

#ifdef XPRS_REPAIRINFEASMAXTIME
    MPD( AddSolverOption_MergeDuplicates("inf:xprs_repairinfeasmaxtime XPRS_REPAIRINFEASMAXTIME",
      "This parameter is deprecated and will be removed in a future release. Overall time limit for the repairinfeas tool"
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  No time limit."
      "\n- (n>0)  If an integer solution has been found, stop MIP search after n seconds, otherwise continue until an integer solution is finally found."
      "\n- (n<0)  Stop in LP or MIP search after n seconds.",
      XPRS_REPAIRINFEASMAXTIME, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_REPAIRINFEASMAXTIME

#ifdef XPRS_REPAIRINFEASTIMELIMIT
    MPD( AddSolverOption_MergeDuplicates("inf:xprs_repairinfeastimelimit XPRS_REPAIRINFEASTIMELIMIT",
      "Overall time limit for the repairinfeas tool"
      "\n\n"
      "Values (default: 1e+20):\n"
      "\n- (>0)  Stop repairinfeas search after the given number of seconds.",
      XPRS_REPAIRINFEASTIMELIMIT, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_REPAIRINFEASTIMELIMIT

#ifdef XPRS_RESOURCESTRATEGY
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_resourcestrategy XPRS_RESOURCESTRATEGY",
      "Controls whether the optimizer is allowed to make nondeterministic decisions if memory is running low in an effort to preserve memory and finish the solve. Available memory (or container limits) are automatically detected but can also be changed by MAXMEMORYSOFT and MAXMEMORYHARD"
      "\n\n"
      "Values (default: 0):\n"
      "\n- (1)  Allow the optimizer to change the solve path if necessary to preserve memory when getting close to one of the memory limits.",
      XPRS_RESOURCESTRATEGY, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_RESOURCESTRATEGY

#ifdef XPRS_RLTCUTS
    MPD( AddSolverOption_MergeDuplicates("cut:xprs_rltcuts XPRS_RLTCUTS",
      "Determines whether RLT cuts should be separated in the Xpress Global Solver."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  The solver decides if RLT cuts are beneficial or not. This is the default setting."
      "\n- (0)  RLT cuts are disabled."
      "\n- (1)  RLT cuts are separated.",
      XPRS_RLTCUTS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_RLTCUTS

#ifdef XPRS_ROOTPRESOLVE
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_rootpresolve XPRS_ROOTPRESOLVE",
      "Determines if presolving should be performed on the problem after the tree search has finished with root cutting and heuristics."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Let the optimizer decide if the problem should be presolved again."
      "\n- (0)  Disabled."
      "\n- (+1)  Always presolve the root problem.",
      XPRS_ROOTPRESOLVE, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_ROOTPRESOLVE

#ifdef XPRS_SBBEST
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_sbbest XPRS_SBBEST",
      "Number of infeasible MIP entities to initialize pseudo costs for on each node."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  determined automatically."
      "\n- (0)  disable strong branching."
      "\n- (n>0)  perform strong branching on up to n entities at each node.",
      XPRS_SBBEST, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_SBBEST

#ifdef XPRS_SBEFFORT
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_sbeffort XPRS_SBEFFORT",
      "Adjusts the overall amount of effort when using strong branching to select an infeasible MIP entity to branch on."
      "\n\nDefault: 1.0",
      XPRS_SBEFFORT, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_SBEFFORT

#ifdef XPRS_SBESTIMATE
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_sbestimate XPRS_SBESTIMATE",
      "Branch and Bound: How to calculate pseudo costs from the local node when selecting an infeasible MIP entity to branch on. These pseudo costs are used in combination with local strong branching and history costs to select the branch candidate."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatically determined."
      "\n- (1-6)  Different variants of local pseudo costs.",
      XPRS_SBESTIMATE, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_SBESTIMATE

#ifdef XPRS_SBITERLIMIT
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_sbiterlimit XPRS_SBITERLIMIT",
      "Number of dual iterations to perform the strong branching for each entity."
      "\n\nDefault: -1 — determined automatically.",
      XPRS_SBITERLIMIT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_SBITERLIMIT

#ifdef XPRS_SBSELECT
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_sbselect XPRS_SBSELECT",
      "The size of the candidate list of MIP entities for strong branching."
      "\n\n"
      "Values (default: -2):\n"
      "\n- (-2)  Automatic (low effort)."
      "\n- (-1)  Automatic (high effort)."
      "\n- (n>=0)  Include n entities in the candidate list (but always at least SBBEST candidates).",
      XPRS_SBSELECT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_SBSELECT

#ifdef XPRS_SCALING
    MPD( AddSolverOption_MergeDuplicates("num:xprs_scaling XPRS_SCALING",
      "This bit-vector control (see Section Bit-vector controls) determines how the Optimizer will rescale a model internally before optimization. If set to 0, no scaling will take place."
      "\n\n"
      "Values (default: 163, meaning bits 0, 1, 5 and 7 are set):\n"
      "\n- (0)  Row scaling."
      "\n- (1)  Column scaling."
      "\n- (2)  Row scaling again."
      "\n- (3)  Maximum."
      "\n- (4)  Curtis-Reid."
      "\n- (5)  0: scale by geometric mean.1: scale by maximum element."
      "\n- (6)  Treat big-M rows as normal rows."
      "\n- (7)  Scale objective function for the simplex method."
      "\n- (8)  Exclude the quadratic part of constraint when calculating scaling factors."
      "\n- (9)  Scale before presolve."
      "\n- (10)  Do not scale rows up."
      "\n- (11)  Do not scale columns down."
      "\n- (12)  Do not apply automatic objective scaling."
      "\n- (13)  RHS scaling."
      "\n- (14)  Disable aggressive quadratic scaling."
      "\n- (15)  Enable explicit linear slack scaling.",
      XPRS_SCALING, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_SCALING

#ifdef XPRS_SDPCUTSTRATEGY
    MPD( AddSolverOption_MergeDuplicates("cut:xprs_sdpcutstrategy XPRS_SDPCUTSTRATEGY",
      " Level of SDP cutting planes separation: This specifies how aggressively SDP cuts should be separated. "
      "\n\n"
      "Values (default:  -1 ):\n"
      "\n- (-1)  Automatic - let the Optimizer decide."
      "\n- (0)  Separation of SDP cuts disabled."
      "\n- (1)  Conservative separation of SDP cuts."
      "\n- (2)  Moderate separation of SDP cuts."
      "\n- (3)  Aggressive separation of SDP cuts.",
      XPRS_SDPCUTSTRATEGY, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_SDPCUTSTRATEGY

#ifdef XPRS_SERIALIZEPREINTSOL
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_serializepreintsol XPRS_SERIALIZEPREINTSOL",
      "Setting SERIALIZEPREINTSOL to 1 will ensure that the preintsol callback is always fired in a deterministic order during a parallel MIP solve. This applies only when the control DETERMINISTIC is set to 1. "
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  The preintsol callbacks will be fired asynchronously from different threads."
      "\n- (1)  The preintsol callbacks will be fired in a deterministic order.",
      XPRS_SERIALIZEPREINTSOL, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_SERIALIZEPREINTSOL

#ifdef XPRS_SIFTING
    MPD( AddSolverOption_MergeDuplicates("lp:xprs_sifting XPRS_SIFTING",
      "Determines whether to enable sifting algorithm with the dual simplex method."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatically determined."
      "\n- (0)  Disable sifting."
      "\n- (1)  Enable sifting.",
      XPRS_SIFTING, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_SIFTING

#ifdef XPRS_SIFTPASSES
    MPD( AddSolverOption_MergeDuplicates("lp:xprs_siftpasses XPRS_SIFTPASSES",
      "Determines how quickly we allow to grow the worker problems during the sifting algorithm. Using larger values can increase the number of columns added to the worker problem which often results in increased solve times for the worker problems but the number of necessary sifting iterations may be reduced. "
      "\n\nDefault: 4",
      XPRS_SIFTPASSES, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_SIFTPASSES

#ifdef XPRS_SIFTPRESOLVEOPS
    MPD( AddSolverOption_MergeDuplicates("lp:xprs_siftpresolveops XPRS_SIFTPRESOLVEOPS",
      "Determines the presolve operations for solving the subproblems during the sifting algorithm."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Use the PRESOLVEOPS setting specified for the original problem."
      "\n- (>=0)  Use the value for the PRESOLVEOPS parameter for solving the subproblems during the sifting algorithm.",
      XPRS_SIFTPRESOLVEOPS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_SIFTPRESOLVEOPS

#ifdef XPRS_SIFTSWITCH
    MPD( AddSolverOption_MergeDuplicates("lp:xprs_siftswitch XPRS_SIFTSWITCH",
      "Determines which algorithm to use for solving the subproblems during sifting."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Dual simplex."
      "\n- (0)  Barrier."
      "\n- (>0)  Use the barrier algorithm while the number of dual infeasibilities is larger than this value, otherwise use dual simplex.",
      XPRS_SIFTSWITCH, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_SIFTSWITCH

#ifdef XPRS_SLEEPONTHREADWAIT
    MPD( AddSolverOption_MergeDuplicates("sys:xprs_sleeponthreadwait XPRS_SLEEPONTHREADWAIT",
      "This parameter is deprecated and will be removed in a future release. In previous versions this was used to determine if the threads should be put into a wait state when waiting for work."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatically determined depending on the CPU the Optimizer is running on."
      "\n- (0)  Keep the threads busy when waiting for work."
      "\n- (1)  Put the threads into a wait state when waiting for work.",
      XPRS_SLEEPONTHREADWAIT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_SLEEPONTHREADWAIT

#ifdef XPRS_SOLTIMELIMIT
    MPD( AddSolverOption_MergeDuplicates("lim:xprs_soltimelimit XPRS_SOLTIMELIMIT",
      "The maximum time in seconds that the Optimizer will run a MIP solve before it terminates, given that a solution has been found. As long as no solution has been found, this control will have no effect. "
      "\n\n"
      "Values (default: 1e+20):\n"
      "\n- (>0)  If an integer solution has been found, stop MIP search after the given number of seconds, otherwise continue until an integer solution is finally found.",
      XPRS_SOLTIMELIMIT, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_SOLTIMELIMIT

#ifdef XPRS_SOSREFTOL
    MPD( AddSolverOption_MergeDuplicates("tol:xprs_sosreftol XPRS_SOSREFTOL",
      "The minimum relative gap between the ordering values of elements in a special ordered set. The gap divided by the absolute value of the larger of the two adjacent values must be at least SOSREFTOL."
      "\n\nDefault: 1.0E-06",
      XPRS_SOSREFTOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_SOSREFTOL

#ifdef XPRS_SYMMETRY
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_symmetry XPRS_SYMMETRY",
      "Adjusts the overall amount of effort for symmetry detection."
      "\n\n"
      "Values (default: 1):\n"
      "\n- (0)  No symmetry detection."
      "\n- (1)  Conservative effort."
      "\n- (2)  Intensive symmetry search.",
      XPRS_SYMMETRY, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_SYMMETRY

#ifdef XPRS_SYMSELECT
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_symselect XPRS_SYMSELECT",
      "Adjusts the overall amount of effort for symmetry detection."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (0)  Search the whole matrix (otherwise the 0, 1 and -1 coefficients only)."
      "\n- (1)  Search all entities (otherwise binaries only).",
      XPRS_SYMSELECT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_SYMSELECT

#ifdef XPRS_THREADS
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_threads XPRS_THREADS",
      "The default number of threads used during optimization."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Determined automatically based on hardware configuration."
      "\n- (>0)  Number of threads to use.",
      XPRS_THREADS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_THREADS

#ifdef XPRS_TIMELIMIT
    MPD( AddSolverOption_MergeDuplicates("lim:xprs_timelimit XPRS_TIMELIMIT",
      "The maximum time in seconds that the Optimizer will run before it terminates, including the problem setup time and solution time. For MIP problems, this is the total time taken to solve all nodes. "
      "\n\n"
      "Values (default: 1e+20):\n"
      "\n- (>0)  Stop LP or MIP search after the given number of seconds.",
      XPRS_TIMELIMIT, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_TIMELIMIT

#ifdef XPRS_TRACE
    MPD( AddSolverOption_MergeDuplicates("pre:xprs_trace XPRS_TRACE",
      "Display the infeasibility diagnosis during presolve. If non-zero, an explanation of the logical deductions made by presolve to deduce infeasibility or unboundedness will be displayed on screen or sent to the message callback function."
      "\n\nDefault: 0",
      XPRS_TRACE, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_TRACE

#ifdef XPRS_TREECOMPRESSION
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_treecompression XPRS_TREECOMPRESSION",
      "When writing nodes to the gloal file, the optimizer can try to use data-compression techniques to reduce the size of the tree file on disk.  The TREECOMPRESSION control determines the strength of the data-compression algorithm used; higher values give superior data-compression at the affect of decreasing performance, while lower values compress quicker but not as effectively.  Where TREECOMPRESSION is set to 0, no data compression will be used on the tree file."
      "\n\nDefault: 2",
      XPRS_TREECOMPRESSION, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_TREECOMPRESSION

#ifdef XPRS_TREECOVERCUTS
    MPD( AddSolverOption_MergeDuplicates("cut:xprs_treecovercuts XPRS_TREECOVERCUTS",
      "Branch and Bound: The number of rounds of lifted cover inequalities generated at nodes other than the root node in the tree. Compare with the description for COVERCUTS.  A value of -1 indicates the number of rounds is determined automatically. "
      "\n\nDefault: -1",
      XPRS_TREECOVERCUTS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_TREECOVERCUTS

#ifdef XPRS_TREECUTSELECT
    MPD( AddSolverOption_MergeDuplicates("cut:xprs_treecutselect XPRS_TREECUTSELECT",
      "A bit-vector (see Section Bit-vector controls) providing detailed control of the cuts created during the tree search of a MIP solve. Use CUTSELECT to control cuts on the root node."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (5)  Clique cuts."
      "\n- (6)  Mixed Integer Rounding (MIR) cuts."
      "\n- (7)  Lifted cover cuts."
      "\n- (8)  Turn on row aggregation for MIR cuts."
      "\n- (11)  Flow path cuts."
      "\n- (12)  Implication cuts."
      "\n- (13)  Turn on automatic Lift and Project cutting strategy."
      "\n- (14)  Disable cutting from cut rows."
      "\n- (15)  Lifted GUB cover cuts."
      "\n- (16)  Zero-half cuts."
      "\n- (17)  Indicator constraint cuts."
      "\n- (18)  Strong Chvatal-Gomory cuts."
      "\n- (20)  Farkas cuts.",
      XPRS_TREECUTSELECT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_TREECUTSELECT

#ifdef XPRS_TREEDIAGNOSTICS
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_treediagnostics XPRS_TREEDIAGNOSTICS",
      "A bit-vector (see Section Bit-vector controls) providing control over how various tree-management-related messages get printed in the tree log file during the branch-and-bound search."
      "\n\n"
      "Values (default: 7):\n"
      "\n- (0)  Output regular summaries of current tree memory usage."
      "\n- (1)  Output messages whenever tree data is being written to tree file."
      "\n- (2)  Output progress messages while tree data is being written to the tree file, at an interval controlled by the TREEFILELOGINTERVAL control.",
      XPRS_TREEDIAGNOSTICS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_TREEDIAGNOSTICS

#ifdef XPRS_TREEFILELOGINTERVAL
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_treefileloginterval XPRS_TREEFILELOGINTERVAL",
      "This control sets the interval between progress messages output while writing tree data to the tree file, in seconds.  The solve is slowed greatly while data is being written to the tree file and this output allows the user to see how much progress is being made. "
      "\n\nDefault: 60",
      XPRS_TREEFILELOGINTERVAL, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_TREEFILELOGINTERVAL

#ifdef XPRS_TREEGOMCUTS
    MPD( AddSolverOption_MergeDuplicates("cut:xprs_treegomcuts XPRS_TREEGOMCUTS",
      "Branch and Bound: The number of rounds of Gomory cuts generated at nodes other than the first node in the tree. Compare with the description for GOMCUTS. A value of -1 indicates the number of rounds is determined automatically."
      "\n\nDefault: -1",
      XPRS_TREEGOMCUTS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_TREEGOMCUTS

#ifdef XPRS_TREEMEMORYLIMIT
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_treememorylimit XPRS_TREEMEMORYLIMIT",
      "A soft limit, in megabytes, for the amount of memory to use in storing the branch and bound search tree.  This doesn't include memory used for presolve, heuristics, solving the LP relaxation, etc. When set to 0 (the default), the optimizer will calculate a limit automatically based on the amount of free physical memory detected in the machine. When the memory used by the branch and bound tree exceeds this limit, the optimizer will try to reduce the memory usage by writing lower-rated sections of the tree to a file called the 'tree file'.  Though the solve can continue if it cannot bring the tree memory usage below the specified limit, performance will be inhibited and a message will be printed to the log."
      "\n\nDefault: 0 (calculate limit automatically)",
      XPRS_TREEMEMORYLIMIT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_TREEMEMORYLIMIT

#ifdef XPRS_TREEMEMORYSAVINGTARGET
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_treememorysavingtarget XPRS_TREEMEMORYSAVINGTARGET",
      "When the memory used by the branch-and-bound search tree exceeds the limit specified by the TREEMEMORYLIMIT control, the optimizer will try to save memory by writing lower-rated sections of the tree to the tree file.  The target amount of memory to save will be enough to bring memory usage back below the limit, plus enough extra to give the tree room to grow. The TREEMEMORYSAVINGTARGET control specifies the extra proportion of the tree's size to try to save; for example, if the tree memory limit is 1000Mb and TREEMEMORYSAVINGTARGET is 0.1, when the tree size exceeds 1000Mb the optimizer will try to reduce the tree size to 900Mb. Reducing the value of TREEMEMORYSAVINGTARGET will cause less extra nodes of the tree to be written to the tree file, but will result in the memory saving routine being triggered more often (as the tree will have less room in which to grow), which can reduce performance.  Increasing the value of TREEMEMORYSAVINGTARGET will cause additional, more highly-rated nodes, of the tree to be written to the tree file, which can cause performance issues if these nodes are required later in the solve."
      "\n\nDefault: 0.4",
      XPRS_TREEMEMORYSAVINGTARGET, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_TREEMEMORYSAVINGTARGET

#ifdef XPRS_TREEQCCUTS
    MPD( AddSolverOption_MergeDuplicates("qp:xprs_treeqccuts XPRS_TREEQCCUTS",
      "Branch and Bound: Limit on the number of rounds of  outer approximation cuts generated for nodes other than the root node, when solving a mixed integer quadratic constrained or mixed integer second order conic problem with outer approximation. "
      "\n\nDefault: -1 — determined automatically.",
      XPRS_TREEQCCUTS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_TREEQCCUTS

#ifdef XPRS_TUNERHISTORY
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_tunerhistory XPRS_TUNERHISTORY",
      "Tuner: Whether to reuse and append to previous tuner results of the same problem. "
      "\n\n"
      "Values (default: 2):\n"
      "\n- (0)  Discard any previous tuner results."
      "\n- (1)  Append new results to the previous tuner results, but do not reuse them."
      "\n- (2)  Reuse the previous results and append new results to it.",
      XPRS_TUNERHISTORY, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_TUNERHISTORY

#ifdef XPRS_TUNERMAXTIME
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_tunermaxtime XPRS_TUNERMAXTIME",
      "Tuner: The maximum time in seconds that the tuner will run before it terminates."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  No time limit."
      "\n- (>0)  Stop the tuner after the given number of seconds.",
      XPRS_TUNERMAXTIME, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_TUNERMAXTIME

#ifdef XPRS_TUNERMETHOD
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_tunermethod XPRS_TUNERMETHOD",
      "Tuner: Selects a factory tuner method. A tuner method consists of a list of controls with different settings that the tuner will evaluate and try to combine."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatically determined. The tuner will select the default method based on the problem type."
      "\n- (0)  Select the default LP tuner method."
      "\n- (1)  Select the default MIP tuner method."
      "\n- (2)  Select a more comprehensive MIP tuner method."
      "\n- (3)  Select a root-focus MIP tuner method."
      "\n- (4)  Select a tree-focus MIP tuner method."
      "\n- (5)  Select a simple MIP tuner method."
      "\n- (6)  Select the default SLP tuner method."
      "\n- (7)  Select the default MISLP tuner method."
      "\n- (8)  Select a MIP tuner method focussed on primal heuristics."
      "\n- (9)  Select the default Xpress Global tuner method.",
      XPRS_TUNERMETHOD, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_TUNERMETHOD

#ifdef XPRS_TUNERMETHODFILE
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_tunermethodfile XPRS_TUNERMETHODFILE",
      "Tuner: Defines a file from which the tuner can read user-defined tuner method."
      "\n\nDefault: (empty)",
      XPRS_TUNERMETHODFILE) );
#endif  // ifdef XPRS_TUNERMETHODFILE

#ifdef XPRS_TUNERMODE
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_tunermode XPRS_TUNERMODE",
      "Tuner: Whether to always enable the tuner or disable it."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  No effect."
      "\n- (0)  Always disable the tuner. XPRStune (TUNE) will have no effect."
      "\n- (1)  Always enable the tuner. XPRSmipoptimize (MIPOPTIMIZE), XPRSlpoptimize (LPOPTIMIZE), etc. will call the tuner before solving the problem.",
      XPRS_TUNERMODE, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_TUNERMODE

#ifdef XPRS_TUNEROUTPUT
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_tuneroutput XPRS_TUNEROUTPUT",
      "Tuner: Whether to output tuner results and logs to the file system."
      "\n\n"
      "Values (default: 1):\n"
      "\n- (0)  Don't output to the file system."
      "\n- (1)  Output results and logs to the file system.",
      XPRS_TUNEROUTPUT, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_TUNEROUTPUT

#ifdef XPRS_TUNEROUTPUTPATH
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_tuneroutputpath XPRS_TUNEROUTPUTPATH",
      "Tuner: Defines a root path to which the tuner writes the result file and logs."
      "\n\nDefault: tuneroutput",
      XPRS_TUNEROUTPUTPATH) );
#endif  // ifdef XPRS_TUNEROUTPUTPATH

#ifdef XPRS_TUNERPERMUTE
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_tunerpermute XPRS_TUNERPERMUTE",
      "Tuner: Defines the number of permutations to solve for each control setting."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  Solve the original problem only for each setting."
      "\n- (n>0)  Solve the original problem and n permuted problems for each setting.",
      XPRS_TUNERPERMUTE, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_TUNERPERMUTE

#ifdef XPRS_TUNERSESSIONNAME
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_tunersessionname XPRS_TUNERSESSIONNAME",
      "Tuner: Defines a session name for the tuner."
      "\n\nDefault: (empty)",
      XPRS_TUNERSESSIONNAME) );
#endif  // ifdef XPRS_TUNERSESSIONNAME

#ifdef XPRS_TUNERTARGET
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_tunertarget XPRS_TUNERTARGET",
      "Tuner: Defines the tuner target -- what should be evaluated when comparing two runs with different control settings."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatically determined. The tuner will choose the default target based on problem type."
      "\n- (0)  Solution time then gap. (MIP/MISLP default)"
      "\n- (1)  Solution time then best bound."
      "\n- (2)  Solution time then best integer solution."
      "\n- (3)  The primal dual integral."
      "\n- (4)  Time only. (LP/SLP default)"
      "\n- (5)  SLP objective only. (SLP/MISLP choice)"
      "\n- (6)  SLP validation number only. (SLP/MISLP choice)"
      "\n- (7)  Gap only."
      "\n- (8)  Best bound only."
      "\n- (9)  Best integer solution only."
      "\n- (10)  Best primal integral. (Only for individual instances, not for problem sets)",
      XPRS_TUNERTARGET, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_TUNERTARGET

#ifdef XPRS_TUNERTHREADS
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_tunerthreads XPRS_TUNERTHREADS",
      "Tuner: the number of threads used by the tuner."
      "\n\n"
      "Values (default: 1):\n"
      "\n- (-1)  Choose automaticlly."
      "\n- (1)  The tuner will run in sequential."
      "\n- (n>1)  The tuner will run in parallel with n threads.",
      XPRS_TUNERTHREADS, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_TUNERTHREADS

#ifdef XPRS_TUNERVERBOSE
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_tunerverbose XPRS_TUNERVERBOSE",
      "Tuner: whether the tuner should prints detailed information for each run."
      "\n\n"
      "Values (default: 1):\n"
      "\n- (1)  Print extra information."
      "\n- (0)  Print less information.",
      XPRS_TUNERVERBOSE, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_TUNERVERBOSE

#ifdef XPRS_USERSOLHEURISTIC
    MPD( AddSolverOption_MergeDuplicates("heur:xprs_usersolheuristic XPRS_USERSOLHEURISTIC",
      " Determines how much effort to put into running a local search heuristic to find a feasible integer solution from a partial or infeasible user solution. "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatically determined."
      "\n- (0)  Search heuristic disabled."
      "\n- (1)  Light effort."
      "\n- (2)  Moderate effort."
      "\n- (3)  High effort.",
      XPRS_USERSOLHEURISTIC, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_USERSOLHEURISTIC

#ifdef XPRS_VARSELECTION
    MPD( AddSolverOption_MergeDuplicates("mip:xprs_varselection XPRS_VARSELECTION",
      "Branch and Bound: This determines the formula used to calculate the estimate of each integer variable, and thus which integer variable is selected to be branched on at a given node. The variable selected to be branched on is the one with the maximum estimate."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Determined automatically."
      "\n- (1)  The minimum of the 'up' and 'down' pseudo costs."
      "\n- (2)  The 'up' pseudo cost plus the 'down' pseudo cost."
      "\n- (3)  The maximum of the 'up' and 'down' pseudo costs, plus twice the minimum of the 'up' and 'down' pseudo costs."
      "\n- (4)  The maximum of the 'up' and 'down' pseudo costs."
      "\n- (5)  The 'down' pseudo cost."
      "\n- (6)  The 'up' pseudo cost."
      "\n- (7)  A weighted combination of the 'up' and 'down' pseudo costs, where the weights depend on how fractional the variable is."
      "\n- (8)  The product of the 'up' and 'down' pseudo costs.",
      XPRS_VARSELECTION, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_VARSELECTION

#ifdef XPRS_VERSION
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_version XPRS_VERSION",
      "The Optimizer version number, e.g. 1301 meaning release 13.01."
      "\n\nDefault: Software version dependent",
      XPRS_VERSION, INT_MIN, INT_MAX) );
#endif  // ifdef XPRS_VERSION

#ifdef XPRS_WORKLIMIT
    MPD( AddSolverOption_MergeDuplicates("tech:xprs_worklimit XPRS_WORKLIMIT",
      "The maximum work (measured in work units) that the Optimizer will run before it terminates. WORK is accumulated during the search and ever increasing. In contrast to TIME, WORK is independent of the hardware and platform on which the search is conducted. The WORKLIMIT serves as a deterministic stopping criterion. When it is reached, it leaves the optimizer in a reproducible state. "
      "\n\n"
      "Values (default: 1e+20):\n"
      "\n- (>0)  Stop LP or MIP search when the given number of work units is reached.",
      XPRS_WORKLIMIT, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XPRS_WORKLIMIT

  }  // AddOptimizerOptions()

};  // class CompiledOptimizerOptions

}  // namespace mp
