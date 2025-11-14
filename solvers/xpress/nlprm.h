#include <climits>
#include <cfloat>

#include "mp/common.h"
#include "mp/error.h"
#include "mp/backend-std.h"


namespace mp {

/// A mix-in class to add Xpress parameters.
/// Translated from '../mp/solvers/xpress/nlprm.h'
/// on Fri Nov 14 16:36:46 2025
///
template <class Impl>
class CompiledNonlinearOptions {
public:
  /// Add up to 254 'Nonlinear' parameters
  void AddNonlinearOptions() {

#ifdef XKTR_PARAM_ALGORITHM
    MPD( AddSolverOption_MergeDuplicates("xktr:param_algorithm XKTR_PARAM_ALGORITHM KNITRO_PARAM_ALGORITHM",
      "Indicates which algorithm to use to solve nonlinear problems"
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)   (auto) let Knitro automatically choose an algorithm, based on the problem characteristics."
      "\n- (1)   (direct) use the Interior/Direct algorithm."
      "\n- (2)   (cg) use the Interior/CG algorithm."
      "\n- (3)   (active) use the Active Set algorithm."
      "\n- (4)   (sqp) use the SQP algorithm."
      "\n- (5)   (multi) run all algorithms, perhaps in parallel."
      "\n- (6)   (al) use the Augmented Lagrangian algorithm.",
      XKTR_PARAM_ALGORITHM, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_ALGORITHM

#ifdef XKTR_PARAM_BAR_DIRECTINTERVAL
    MPD( AddSolverOption_MergeDuplicates("xktr:param_bar_directinterval XKTR_PARAM_BAR_DIRECTINTERVAL KNITRO_PARAM_BAR_DIRECTINTERVAL",
      "Controls the maximum number of consecutive conjugate gradient (CG) steps before Knitro will try to enforce that a step is taken using direct linear algebra. "
      "\n\nDefault: 10",
      XKTR_PARAM_BAR_DIRECTINTERVAL, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_BAR_DIRECTINTERVAL

#ifdef XKTR_PARAM_BAR_FEASIBLE
    MPD( AddSolverOption_MergeDuplicates("xktr:param_bar_feasible XKTR_PARAM_BAR_FEASIBLE KNITRO_PARAM_BAR_FEASIBLE",
      "Specifies whether special emphasis is placed on getting and staying feasible in the interior-point algorithms."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)   (no) No special emphasis on feasibility."
      "\n- (1)   (stay) Iterates must satisfy inequality constraints once they become sufficiently feasible."
      "\n- (2)   (get) Special emphasis is placed on getting feasible before trying to optimize."
      "\n- (3)   (get_stay) Implement both options 1 and 2 above.",
      XKTR_PARAM_BAR_FEASIBLE, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_BAR_FEASIBLE

#ifdef XKTR_PARAM_BAR_FEASMODETOL
    MPD( AddSolverOption_MergeDuplicates("xktr:param_bar_feasmodetol XKTR_PARAM_BAR_FEASMODETOL KNITRO_PARAM_BAR_FEASMODETOL",
      "Specifies the tolerance in equation that determines whether Knitro will force subsequent iterates to remain feasible. "
      "\n\nDefault: 1.0e-4",
      XKTR_PARAM_BAR_FEASMODETOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XKTR_PARAM_BAR_FEASMODETOL

#ifdef XKTR_PARAM_BAR_INITMU
    MPD( AddSolverOption_MergeDuplicates("xktr:param_bar_initmu XKTR_PARAM_BAR_INITMU KNITRO_PARAM_BAR_INITMU",
      "Specifies the initial value for the barrier parameter : μ used with the barrier algorithms. This option has no effect on the Active Set algorithm. "
      "\n\nDefault: 1.0e-1",
      XKTR_PARAM_BAR_INITMU, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XKTR_PARAM_BAR_INITMU

#ifdef XKTR_PARAM_BAR_INITPT
    MPD( AddSolverOption_MergeDuplicates("xktr:param_bar_initpt XKTR_PARAM_BAR_INITPT KNITRO_PARAM_BAR_INITPT",
      "Indicates whether an initial point strategy is used with barrier algorithms."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)   (auto) Let Knitro automatically choose the strategy."
      "\n- (1)   (yes) Shift the initial slacks and multipliers to improve barrier algorithm performance."
      "\n- (2)   (no) Do no alter the initial slacks and multipliers.",
      XKTR_PARAM_BAR_INITPT, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_BAR_INITPT

#ifdef XKTR_PARAM_BAR_MAXBACKTRACK
    MPD( AddSolverOption_MergeDuplicates("xktr:param_bar_maxbacktrack XKTR_PARAM_BAR_MAXBACKTRACK KNITRO_PARAM_BAR_MAXBACKTRACK",
      "Indicates the maximum allowable number of backtracks during the linesearch of the Interior/Direct algorithm before reverting to a CG step. "
      "\n\nDefault: 3",
      XKTR_PARAM_BAR_MAXBACKTRACK, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_BAR_MAXBACKTRACK

#ifdef XKTR_PARAM_BAR_MAXCROSSIT
    MPD( AddSolverOption_MergeDuplicates("xktr:param_bar_maxcrossit XKTR_PARAM_BAR_MAXCROSSIT KNITRO_PARAM_BAR_MAXCROSSIT",
      "Specifies the maximum number of crossover iterations before termination. "
      "\n\nDefault: 0",
      XKTR_PARAM_BAR_MAXCROSSIT, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_BAR_MAXCROSSIT

#ifdef XKTR_PARAM_BAR_MAXREFACTOR
    MPD( AddSolverOption_MergeDuplicates("xktr:param_bar_maxrefactor XKTR_PARAM_BAR_MAXREFACTOR KNITRO_PARAM_BAR_MAXREFACTOR",
      "Indicates the maximum number of refactorizations of the KKT system per iteration of the Interior/Direct algorithm before reverting to a CG step. "
      "\n\nDefault: -1",
      XKTR_PARAM_BAR_MAXREFACTOR, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_BAR_MAXREFACTOR

#ifdef XKTR_PARAM_BAR_MURULE
    MPD( AddSolverOption_MergeDuplicates("xktr:param_bar_murule XKTR_PARAM_BAR_MURULE KNITRO_PARAM_BAR_MURULE",
      "Indicates which strategy to use for modifying the barrier parameter mu in the barrier algorithms."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)    (auto) Let Knitro automatically choose the strategy."
      "\n- (1)    (monotone) Monotonically decrease the barrier parameter. Available for both barrier algorithms."
      "\n- (2)    (adaptive) Use an adaptive rule based on the complementarity gap to determine the value of the barrier parameter. Available for both barrier algorithms."
      "\n- (3)    (probing) Use a probing (affine-scaling) step to dynamically determine the barrier parameter. Available only for the Interior/Direct algorithm."
      "\n- (4)    (dampmpc) Use a Mehrotra predictor-corrector type rule to determine the barrier parameter, with safeguards on the corrector step. Available only for the Interior/Direct algorithm."
      "\n- (5)    (fullmpc) Use a Mehrotra predictor-corrector type rule to determine the barrier parameter, without safeguards on the corrector step. Available only for the Interior/Direct algorithm."
      "\n- (6)    (quality) Minimize a quality function at each iteration to determine the barrier parameter. Available only for the Interior/Direct algorithm.",
      XKTR_PARAM_BAR_MURULE, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_BAR_MURULE

#ifdef XKTR_PARAM_BAR_PENCONS
    MPD( AddSolverOption_MergeDuplicates("xktr:param_bar_pencons XKTR_PARAM_BAR_PENCONS KNITRO_PARAM_BAR_PENCONS",
      "Indicates whether a penalty approach is applied to the constraints."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)    (auto) Let Knitro automatically choose the strategy."
      "\n- (1)    (none) No constraints are penalized."
      "\n- (2)    (all) A penalty approach is applied to all general constraints.",
      XKTR_PARAM_BAR_PENCONS, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_BAR_PENCONS

#ifdef XKTR_PARAM_BAR_PENRULE
    MPD( AddSolverOption_MergeDuplicates("xktr:param_bar_penrule XKTR_PARAM_BAR_PENRULE KNITRO_PARAM_BAR_PENRULE",
      "Indicates which penalty parameter strategy to use for determining whether or not to accept a trial iterate."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)   (auto) Let Knitro automatically choose the strategy."
      "\n- (1)   (single) Use a single penalty parameter in the merit function to weight feasibility versus optimality."
      "\n- (2)   (flex) Use a more tolerant and flexible step acceptance procedure based on a range of penalty parameter values.",
      XKTR_PARAM_BAR_PENRULE, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_BAR_PENRULE

#ifdef XKTR_PARAM_BAR_SWITCHRULE
    MPD( AddSolverOption_MergeDuplicates("xktr:param_bar_switchrule XKTR_PARAM_BAR_SWITCHRULE KNITRO_PARAM_BAR_SWITCHRULE",
      "Indicates whether or not the barrier algorithms will allow switching from an optimality phase to a pure feasibility phase. "
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)   (auto) Let Knitro determine the switching procedure."
      "\n- (1)   (never) Never switch to feasibility phase."
      "\n- (2)   (level1) Allow switches to feasibility phase."
      "\n- (3)   (level2) Use a more aggressive switching rule.",
      XKTR_PARAM_BAR_SWITCHRULE, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_BAR_SWITCHRULE

#ifdef XKTR_PARAM_DELTA
    MPD( AddSolverOption_MergeDuplicates("xktr:param_delta XKTR_PARAM_DELTA KNITRO_PARAM_DELTA",
      "Specifies the initial trust region radius scaling factor used to determine the initial trust region size. "
      "\n\nDefault: 1.0e0",
      XKTR_PARAM_DELTA, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XKTR_PARAM_DELTA

#ifdef XKTR_PARAM_FEASTOL
    MPD( AddSolverOption_MergeDuplicates("xktr:param_feastol XKTR_PARAM_FEASTOL KNITRO_PARAM_FEASTOL",
      "Specifies the final relative stopping tolerance for the feasibility error."
      "\n\nDefault: 1.0e-6",
      XKTR_PARAM_FEASTOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XKTR_PARAM_FEASTOL

#ifdef XKTR_PARAM_FEASTOLABS
    MPD( AddSolverOption_MergeDuplicates("xktr:param_feastolabs XKTR_PARAM_FEASTOLABS KNITRO_PARAM_FEASTOLABS",
      "Specifies the final absolute stopping tolerance for the feasibility error."
      "\n\nDefault: 0.0e0",
      XKTR_PARAM_FEASTOLABS, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XKTR_PARAM_FEASTOLABS

#ifdef XKTR_PARAM_GRADOPT
    MPD( AddSolverOption_MergeDuplicates("xktr:param_gradopt XKTR_PARAM_GRADOPT KNITRO_PARAM_GRADOPT",
      "Specifies how to compute the gradients of the objective and constraint functions. "
      "\n\n"
      "Values (default: 1):\n"
      "\n- (1)   (exact) User provides a routine for computing the exact gradients."
      "\n- (2)   (forward) Knitro computes gradients by forward finite-differences."
      "\n- (3)   (central) Knitro computes gradients by central finite differences.",
      XKTR_PARAM_GRADOPT, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_GRADOPT

#ifdef XKTR_PARAM_HESSOPT
    MPD( AddSolverOption_MergeDuplicates("xktr:param_hessopt XKTR_PARAM_HESSOPT KNITRO_PARAM_HESSOPT",
      "Specifies how to compute the (approximate) Hessian of the Lagrangian. "
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)   (auto) Let Knitro make an automatic choice."
      "\n- (1)   (exact) User provides a routine for computing the exact Hessian."
      "\n- (2)   (bfgs) Knitro computes a (dense) quasi-Newton BFGS Hessian."
      "\n- (3)   (sr1) Knitro computes a (dense) quasi-Newton SR1 Hessian."
      "\n- (4)   (finite_diff) Knitro computes Hessian-vector products using finite-differences."
      "\n- (5)   (product) User provides a routine to compute the Hessian-vector products."
      "\n- (6)   (lbfgs) Knitro computes a limited-memory quasi-Newton BFGS Hessian (its size is determined by the option lmsize)."
      "\n- (7)   (gauss_newton) Knitro computes a Gauss-Newton approximation of the hessian (available for least-squares only, and default value for least-squares)",
      XKTR_PARAM_HESSOPT, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_HESSOPT

#ifdef XKTR_PARAM_HONORBNDS
    MPD( AddSolverOption_MergeDuplicates("xktr:param_honorbnds XKTR_PARAM_HONORBNDS KNITRO_PARAM_HONORBNDS",
      "Indicates whether or not to enforce satisfaction of simple variable bounds throughout the optimization.  "
      "\n\n"
      "Values (default: 2):\n"
      "\n- (0)   (no) Knitro does not require that the bounds on the variables be satisfied at intermediate iterates."
      "\n- (1)   (always) Knitro enforces that the initial point and all subsequent solution estimates satisfy the bounds on the variables."
      "\n- (2)   (initpt) Knitro enforces that the initial point satisfies the bounds on the variables.",
      XKTR_PARAM_HONORBNDS, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_HONORBNDS

#ifdef XKTR_PARAM_INFEASTOL
    MPD( AddSolverOption_MergeDuplicates("xktr:param_infeastol XKTR_PARAM_INFEASTOL KNITRO_PARAM_INFEASTOL",
      "Specifies the (relative) tolerance used for declaring infeasibility of a model."
      "\n\nDefault: 1.0e-8",
      XKTR_PARAM_INFEASTOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XKTR_PARAM_INFEASTOL

#ifdef XKTR_PARAM_LMSIZE
    MPD( AddSolverOption_MergeDuplicates("xktr:param_lmsize XKTR_PARAM_LMSIZE KNITRO_PARAM_LMSIZE",
      "Specifies the number of limited memory pairs stored when approximating the Hessian using the limited-memory quasi-Newton BFGS option."
      "\n\nDefault: 10",
      XKTR_PARAM_LMSIZE, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_LMSIZE

#ifdef XKTR_PARAM_MAXCGIT
    MPD( AddSolverOption_MergeDuplicates("xktr:param_maxcgit XKTR_PARAM_MAXCGIT KNITRO_PARAM_MAXCGIT",
      "Specifies the number of limited memory pairs stored when approximating the Hessian using the limited-memory quasi-Newton BFGS option."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)   Let Knitro automatically choose a value based on the problem size."
      "\n- (n)   At most n>0 CG iterations may be performed during one minor iteration of Knitro.",
      XKTR_PARAM_MAXCGIT, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_MAXCGIT

#ifdef XKTR_PARAM_MAXIT
    MPD( AddSolverOption_MergeDuplicates("xktr:param_maxit XKTR_PARAM_MAXIT KNITRO_PARAM_MAXIT",
      "Specifies the maximum number of iterations before termination."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)   Let Knitro automatically choose a value based on the problem type. Currently Knitro sets this value to 10000 for LPs/NLPs and 3000 for MIP problems."
      "\n- (n)   At most n>0 iterations may be performed before terminating.",
      XKTR_PARAM_MAXIT, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_MAXIT

#ifdef XKTR_PARAM_MIP_BRANCHRULE
    MPD( AddSolverOption_MergeDuplicates("xktr:param_mip_branchrule XKTR_PARAM_MIP_BRANCHRULE KNITRO_PARAM_MIP_BRANCHRULE",
      "Specifies which branching rule to use for MIP branch and bound procedure."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)   (auto) Let Knitro automatically choose the branching rule."
      "\n- (1)   (most_frac) Use most fractional (most infeasible) branching."
      "\n- (2)   (pseudcost) Use pseudo-cost branching."
      "\n- (3)   (strong) Use strong branching (see options XKTR_PARAM_MIP_STRONG_CANDLIM, XKTR_PARAM_MIP_STRONG_LEVEL,  XKTR_PARAM_MIP_STRONG_MAXIT for further control of strong branching procedure).",
      XKTR_PARAM_MIP_BRANCHRULE, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_MIP_BRANCHRULE

#ifdef XKTR_PARAM_MIP_GUB_BRANCH
    MPD( AddSolverOption_MergeDuplicates("xktr:param_mip_gub_branch XKTR_PARAM_MIP_GUB_BRANCH KNITRO_PARAM_MIP_GUB_BRANCH",
      "Specifies whether or not to branch on generalized upper bounds (GUBs)."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)   (no) Do not branch on GUBs."
      "\n- (1)   (yes) Allow branching on GUBs.",
      XKTR_PARAM_MIP_GUB_BRANCH, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_MIP_GUB_BRANCH

#ifdef XKTR_PARAM_MIP_HEURISTIC
    MPD( AddSolverOption_MergeDuplicates("xktr:param_mip_heuristic XKTR_PARAM_MIP_HEURISTIC KNITRO_PARAM_MIP_HEURISTIC",
      "Specifies which MIP heuristic search approach to apply to try to find an initial integer feasible point."
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)   (auto) Let Knitro choose the heuristic to apply (if any)."
      "\n- (1)   (none) No heuristic search applied."
      "\n- (2)   (feaspump) Apply feasibility pump heuristic."
      "\n- (3)   (mpec) Apply heuristic based on MPEC formulation.",
      XKTR_PARAM_MIP_HEURISTIC, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_MIP_HEURISTIC

#ifdef XKTR_PARAM_MIP_HEURISTIC_MAXIT
    MPD( AddSolverOption_MergeDuplicates("xktr:param_mip_heuristic_maxit XKTR_PARAM_MIP_HEURISTIC_MAXIT KNITRO_PARAM_MIP_HEURISTIC_MAXIT",
      "Specifies the maximum number of iterations to allow for MIP heuristic, if one is enabled."
      "\n\nDefault: 100",
      XKTR_PARAM_MIP_HEURISTIC_MAXIT, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_MIP_HEURISTIC_MAXIT

#ifdef XKTR_PARAM_MIP_IMPLICATNS
    MPD( AddSolverOption_MergeDuplicates("xktr:param_mip_implicatns XKTR_PARAM_MIP_IMPLICATNS KNITRO_PARAM_MIP_IMPLICATNS",
      "Specifies whether or not to add constraints to the MIP derived from logical implications."
      "\n\n"
      "Values (default: 1):\n"
      "\n- (0)   (no) Do not add constraints from logical implications."
      "\n- (1)   (yes) Knitro adds constraints from logical implications.",
      XKTR_PARAM_MIP_IMPLICATNS, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_MIP_IMPLICATNS

#ifdef XKTR_PARAM_MIP_INTEGERTOL
    MPD( AddSolverOption_MergeDuplicates("xktr:param_mip_integertol XKTR_PARAM_MIP_INTEGERTOL KNITRO_PARAM_INTEGERTOL",
      "This value specifies the threshold for deciding whether or not a variable is determined to be an integer."
      "\n\nDefault: 1.0e-8",
      XKTR_PARAM_MIP_INTEGERTOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XKTR_PARAM_MIP_INTEGERTOL

#ifdef XKTR_PARAM_MIP_INTGAPABS
    MPD( AddSolverOption_MergeDuplicates("xktr:param_mip_intgapabs XKTR_PARAM_MIP_INTGAPABS KNITRO_PARAM_INTGAPABS",
      "The absolute integrality gap stop tolerance for MIP."
      "\n\nDefault: 1.0e-6",
      XKTR_PARAM_MIP_INTGAPABS, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XKTR_PARAM_MIP_INTGAPABS

#ifdef XKTR_PARAM_MIP_INTGAPREL
    MPD( AddSolverOption_MergeDuplicates("xktr:param_mip_intgaprel XKTR_PARAM_MIP_INTGAPREL KNITRO_PARAM_INTGAPREL",
      "The relative integrality gap stop tolerance for MIP."
      "\n\nDefault: 1.0e-6",
      XKTR_PARAM_MIP_INTGAPREL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XKTR_PARAM_MIP_INTGAPREL

#ifdef XKTR_PARAM_MIP_KNAPSACK
    MPD( AddSolverOption_MergeDuplicates("xktr:param_mip_knapsack XKTR_PARAM_MIP_KNAPSACK KNITRO_PARAM_MIP_KNAPSACK",
      "Specifies rules for adding MIP knapsack cuts."
      "\n\n"
      "Values (default: 1):\n"
      "\n- (0)   (none) Do not add knapsack cuts."
      "\n- (1)   (ineqs) Add cuts derived from inequalities only."
      "\n- (2)   (ineqs_eqs) Add cuts derived from both inequalities and equalities.",
      XKTR_PARAM_MIP_KNAPSACK, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_MIP_KNAPSACK

#ifdef XKTR_PARAM_MIP_LPALG
    MPD( AddSolverOption_MergeDuplicates("xktr:param_mip_lpalg XKTR_PARAM_MIP_LPALG KNITRO_PARAM_MIP_LPALG",
      "Specifies which algorithm to use for any linear programming (LP) subproblem solves that may occur in the MIP branch and bound procedure. "
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)   (auto) Let Knitro automatically choose an algorithm, based on the problem characteristics."
      "\n- (1)   (direct) Use the Interior/Direct (barrier) algorithm."
      "\n- (2)   (cg) Use the Interior/CG (barrier) algorithm."
      "\n- (3)   (active) Use the Active Set (simplex) algorithm.",
      XKTR_PARAM_MIP_LPALG, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_MIP_LPALG

#ifdef XKTR_PARAM_MIP_MAXNODES
    MPD( AddSolverOption_MergeDuplicates("xktr:param_mip_maxnodes XKTR_PARAM_MIP_MAXNODES KNITRO_PARAM_MIP_MAXNODES",
      "Specifies the maximum number of nodes explored."
      "\n\nDefault: 100000",
      XKTR_PARAM_MIP_MAXNODES, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_MIP_MAXNODES

#ifdef XKTR_PARAM_MIP_MAXSOLVES
    MPD( AddSolverOption_MergeDuplicates("xktr:param_mip_maxsolves XKTR_PARAM_MIP_MAXSOLVES KNITRO_PARAM_MIP_MAXSOLVES",
      "Specifies the maximum number of subproblem solves allowed (0 means no limit)."
      "\n\nDefault: 200000",
      XKTR_PARAM_MIP_MAXSOLVES, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_MIP_MAXSOLVES

#ifdef XKTR_PARAM_MIP_METHOD
    MPD( AddSolverOption_MergeDuplicates("xktr:param_mip_method XKTR_PARAM_MIP_METHOD KNITRO_PARAM_MIP_METHOD",
      "Specifies which MIP method to use. "
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)   (auto) Let Knitro automatically choose the method."
      "\n- (1)   (BB) Use the standard branch and bound method."
      "\n- (2)   (HQG) Use the hybrid Quesada-Grossman method (for convex, nonlinear problems only).",
      XKTR_PARAM_MIP_METHOD, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_MIP_METHOD

#ifdef XKTR_PARAM_MIP_OUTINTERVAL
    MPD( AddSolverOption_MergeDuplicates("xktr:param_mip_outinterval XKTR_PARAM_MIP_OUTINTERVAL KNITRO_PARAM_MIP_OUTINTERVAL",
      "Specifies node printing interval for XKTR_PARAM_MIP_OUTLEVEL when XKTR_PARAM_MIP_OUTLEVEL > 0. "
      "\n\n"
      "Values (default: 10):\n"
      "\n- (0)   Print output every node."
      "\n- (2)   Print output every 2nd node."
      "\n- (N)   Print output every Nth node.",
      XKTR_PARAM_MIP_OUTINTERVAL, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_MIP_OUTINTERVAL

#ifdef XKTR_PARAM_MIP_OUTLEVEL
    MPD( AddSolverOption_MergeDuplicates("xktr:param_mip_outlevel XKTR_PARAM_MIP_OUTLEVEL KNITRO_PARAM_MIP_OUTLEVEL",
      "Specifies how much MIP information to print. "
      "\n\n"
      "Values (default: 1):\n"
      "\n- (0)   (none) Do not print any MIP node information."
      "\n- (1)   (iters) Print one line of output for every node.",
      XKTR_PARAM_MIP_OUTLEVEL, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_MIP_OUTLEVEL

#ifdef XKTR_PARAM_MIP_PSEUDOINIT
    MPD( AddSolverOption_MergeDuplicates("xktr:param_mip_pseudoinit XKTR_PARAM_MIP_PSEUDOINIT KNITRO_PARAM_MIP_PSEUDOINIT",
      "Specifies the method used to initialize pseudo-costs corresponding to variables that have not yet been branched on in the MIP method. "
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)   Let Knitro automatically choose the method."
      "\n- (1)   Initialize using the average value of computed pseudo-costs."
      "\n- (2)   Initialize using strong branching.",
      XKTR_PARAM_MIP_PSEUDOINIT, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_MIP_PSEUDOINIT

#ifdef XKTR_PARAM_MIP_ROOTALG
    MPD( AddSolverOption_MergeDuplicates("xktr:param_mip_rootalg XKTR_PARAM_MIP_ROOTALG KNITRO_PARAM_MIP_ROOTALG",
      "Specifies which algorithm to use for the root node solve in MIP (same options as XKTR_PARAM_ALGORITHM user option). "
      "\n\nDefault: 0",
      XKTR_PARAM_MIP_ROOTALG, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_MIP_ROOTALG

#ifdef XKTR_PARAM_MIP_ROUNDING
    MPD( AddSolverOption_MergeDuplicates("xktr:param_mip_rounding XKTR_PARAM_MIP_ROUNDING KNITRO_PARAM_MIP_ROUNDING",
      "Specifies the MIP rounding rule to apply. "
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)   (auto) Let Knitro choose the rounding rule."
      "\n- (1)   (none) Do not round if a node is infeasible."
      "\n- (2)   (heur_only) Round using a fast heuristic only."
      "\n- (3)   (nlp_sometimes) Round and solve a subproblem if likely to succeed."
      "\n- (4)   (nlp_always) Always round and solve a subproblem.",
      XKTR_PARAM_MIP_ROUNDING, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_MIP_ROUNDING

#ifdef XKTR_PARAM_MIP_SELECTRULE
    MPD( AddSolverOption_MergeDuplicates("xktr:param_mip_selectrule XKTR_PARAM_MIP_SELECTRULE KNITRO_PARAM_MIP_SELECTRULE",
      "Specifies the MIP select rule for choosing the next node in the branch and bound tree. "
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)   (auto) Let Knitro choose the node selection rule."
      "\n- (1)   (depth_first) Search the tree using a depth first procedure."
      "\n- (2)   (best_bound) Select the node with the best relaxation bound."
      "\n- (3)   (combo_1) Use depth first unless pruned, then best bound.",
      XKTR_PARAM_MIP_SELECTRULE, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_MIP_SELECTRULE

#ifdef XKTR_PARAM_MIP_STRONG_CANDLIM
    MPD( AddSolverOption_MergeDuplicates("xktr:param_mip_strong_candlim XKTR_PARAM_MIP_STRONG_CANDLIM KNITRO_PARAM_MIP_STRONG_CANDLIM",
      "Specifies the maximum number of candidates to explore for MIP strong branching. "
      "\n\nDefault: 10",
      XKTR_PARAM_MIP_STRONG_CANDLIM, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_MIP_STRONG_CANDLIM

#ifdef XKTR_PARAM_MIP_STRONG_LEVEL
    MPD( AddSolverOption_MergeDuplicates("xktr:param_mip_strong_level XKTR_PARAM_MIP_STRONG_LEVEL KNITRO_PARAM_MIP_STRONG_LEVEL",
      "Specifies the maximum number of tree levels on which to perform MIP strong branching. "
      "\n\nDefault: 10",
      XKTR_PARAM_MIP_STRONG_LEVEL, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_MIP_STRONG_LEVEL

#ifdef XKTR_PARAM_MIP_STRONG_MAXIT
    MPD( AddSolverOption_MergeDuplicates("xktr:param_mip_strong_maxit XKTR_PARAM_MIP_STRONG_MAXIT KNITRO_PARAM_MIP_STRONG_MAXIT",
      "Specifies the maximum number of iterations to allow for MIP strong branching solves. "
      "\n\nDefault: 1000",
      XKTR_PARAM_MIP_STRONG_MAXIT, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_MIP_STRONG_MAXIT

#ifdef XKTR_PARAM_MIP_TERMINATE
    MPD( AddSolverOption_MergeDuplicates("xktr:param_mip_terminate XKTR_PARAM_MIP_TERMINATE KNITRO_PARAM_MIP_TERMINATE",
      "Specifies conditions for terminating the MIP algorithm. "
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)   (optimal) Terminate at optimum."
      "\n- (1)   (feasible) Terminate at first integer feasible point.",
      XKTR_PARAM_MIP_TERMINATE, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_MIP_TERMINATE

#ifdef XKTR_PARAM_OBJRANGE
    MPD( AddSolverOption_MergeDuplicates("xktr:param_objrange XKTR_PARAM_OBJRANGE KNITRO_PARAM_OBJRANGE",
      "Specifies the extreme limits of the objective function for purposes of determining unboundedness."
      "\n\nDefault: 1.0e20",
      XKTR_PARAM_OBJRANGE, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XKTR_PARAM_OBJRANGE

#ifdef XKTR_PARAM_OPTTOL
    MPD( AddSolverOption_MergeDuplicates("xktr:param_opttol XKTR_PARAM_OPTTOL KNITRO_PARAM_OPTTOL",
      "Specifies the final relative stopping tolerance for the KKT (optimality) error."
      "\n\nDefault: 1.0e-6",
      XKTR_PARAM_OPTTOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XKTR_PARAM_OPTTOL

#ifdef XKTR_PARAM_OPTTOLABS
    MPD( AddSolverOption_MergeDuplicates("xktr:param_opttolabs XKTR_PARAM_OPTTOLABS KNITRO_PARAM_OPTTOLABS",
      "Specifies the final absolute stopping tolerance for the KKT (optimality) error."
      "\n\nDefault: 0.0e0",
      XKTR_PARAM_OPTTOLABS, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XKTR_PARAM_OPTTOLABS

#ifdef XKTR_PARAM_OUTLEV
    MPD( AddSolverOption_MergeDuplicates("xktr:param_outlev XKTR_PARAM_OUTLEV KNITRO_PARAM_OUTLEV",
      "Controls the level of output produced by Knitro. "
      "\n\n"
      "Values (default: 2):\n"
      "\n- (0)   (none) Printing of all output is suppressed."
      "\n- (1)   (summary) Print only summary information."
      "\n- (2)   (iter_10) Print basic information every 10 iterations."
      "\n- (3)   (iter) Print basic information at each iteration."
      "\n- (4)   (iter_verbose) Print basic information and the function count at each iteration."
      "\n- (5)   (iter_x) Print all the above, and the values of the solution vector x."
      "\n- (6)   (all) Print all the above, and the values of the constraints c at x and the Lagrange multipliers lambda.",
      XKTR_PARAM_OUTLEV, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_OUTLEV

#ifdef XKTR_PARAM_PRESOLVE
    MPD( AddSolverOption_MergeDuplicates("xktr:param_presolve XKTR_PARAM_PRESOLVE KNITRO_PARAM_PRESOLVE",
      "Determine whether or not to use the Knitro presolver to try to simplify the model by removing variables or constraints. Specifies conditions for terminating the MIP algorithm. "
      "\n\n"
      "Values (default: 1):\n"
      "\n- (0)   (none) Do not use Knitro presolver."
      "\n- (1)   (basic) Use the Knitro basic presolver.",
      XKTR_PARAM_PRESOLVE, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_PRESOLVE

#ifdef XKTR_PARAM_PRESOLVE_TOL
    MPD( AddSolverOption_MergeDuplicates("xktr:param_presolve_tol XKTR_PARAM_PRESOLVE_TOL KNITRO_PARAM_PRESOLVE_TOL",
      "Determines the tolerance used by the Knitro presolver to remove variables and constraints from the model."
      "\n\nDefault: 1.0e-6",
      XKTR_PARAM_PRESOLVE_TOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XKTR_PARAM_PRESOLVE_TOL

#ifdef XKTR_PARAM_SCALE
    MPD( AddSolverOption_MergeDuplicates("xktr:param_scale XKTR_PARAM_SCALE KNITRO_PARAM_SCALE",
      "Performs a scaling of the objective and constraint functions based on their values at the initial point. "
      "\n\n"
      "Values (default: 1):\n"
      "\n- (0)   (no) No scaling is performed."
      "\n- (1)   (yes) Knitro is allowed to scale the objective function and constraints.",
      XKTR_PARAM_SCALE, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_SCALE

#ifdef XKTR_PARAM_SOC
    MPD( AddSolverOption_MergeDuplicates("xktr:param_soc XKTR_PARAM_SOC KNITRO_PARAM_SOC",
      "Specifies whether or not to try second order corrections (SOC). "
      "\n\n"
      "Values (default: 1):\n"
      "\n- (0)   (no) No second order correction steps are attempted."
      "\n- (1)   (maybe) Second order correction steps may be attempted on some iterations."
      "\n- (2)   (yes) Second order correction steps are always attempted if the original step is rejected and there are nonlinear constraints. ",
      XKTR_PARAM_SOC, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_SOC

#ifdef XKTR_PARAM_SOLTYPE
    MPD( AddSolverOption_MergeDuplicates("xktr:param_soltype XKTR_PARAM_SOLTYPE KNITRO_PARAM_SOLTYPE",
      "This option specifies the solution returned by Knitro. Generally, the solution converged to by Knitro is a locally optimal solution that corresponds to the best feasible solution found. However, on rare occasions, Knitro may enounter a feasible solution during the optimization process that has a better objective value than the final solution converged to by Knitro. Setting soltype = 1 in this case will return this iterate. "
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)   (final) Always return the final solution to which Knitro converges."
      "\n- (1)   (bestfeas) Always return the best feasible solution encountered during the optimization.",
      XKTR_PARAM_SOLTYPE, INT_MIN, INT_MAX) );
#endif  // ifdef XKTR_PARAM_SOLTYPE

#ifdef XKTR_PARAM_XTOL
    MPD( AddSolverOption_MergeDuplicates("xktr:param_xtol XKTR_PARAM_XTOL KNITRO_PARAM_XTOL",
      "The optimization process will terminate if the relative change in all components of the solution point estimate is less than xtol. "
      "\n\nDefault: 1.0e-15",
      XKTR_PARAM_XTOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XKTR_PARAM_XTOL

#ifdef XSLP_ALGORITHM
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_algorithm XSLP_ALGORITHM SLPALGORITHM",
      "Bit map describing the SLP algorithm(s) to be used"
      "\n\n"
      "Values (default: 166 (sets bits 1, 2, 5, 7)):\n"
      "\n- (0)  Do not apply step bounds."
      "\n- (1)  Apply step bounds to SLP delta vectors only when required."
      "\n- (2)  Estimate step bounds from early SLP iterations."
      "\n- (3)  Use dynamic damping."
      "\n- (4)  Do not update values which are converged within strict tolerance."
      "\n- (5)  Retain previous value when cascading if determining row is zero."
      "\n- (6)  Reset XSLP_DELTA_Z to zero when converged and continue SLP."
      "\n- (7)  Quick convergence check."
      "\n- (8)  Escalate penalties."
      "\n- (9)  Use the primal simplex algorithm when all error vectors become inactive."
      "\n- (11)  Continue optimizing after penalty cost reaches maximum."
      "\n- (12)  Accept a solution which has converged even if there are still significant active penalty error vectors."
      "\n- (13)  Skip the solution polishing step if the LP postsolve returns a slightly infeasible, but claimed optimal solution."
      "\n- (14)  Step bounds are updated to accomodate cascaded values (otherwise cascaded values are pushed to respect step bounds)."
      "\n- (15)  Apply clamping when converged on extended criteria only with some variables having active step bounds."
      "\n- (16)  Apply clamping when converged on extended criteria only.",
      XSLP_ALGORITHM, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_ALGORITHM

#ifdef XSLP_ANALYZE
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_analyze XSLP_ANALYZE SLPANALYZE",
      "Bit map activating additional options supporting model / solution path analysis"
      "\n\n"
      "Values (default: 0):\n"
      "\n- (3)  Include an extended iteration summary."
      "\n- (4)  Run infeasibility analysis on infeasible iterations."
      "\n- (6)  Write the linearizations to disk at every XSLP_AUTOSAVE iterations."
      "\n- (7)  Write the initial basis of the linearizations to disk at every XSLP_AUTOSAVE iterations."
      "\n- (8)  Create an XSLP save file at every XSLP_AUTOSAVE iterations.",
      XSLP_ANALYZE, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_ANALYZE

#ifdef XSLP_ATOL_A
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_atol_a XSLP_ATOL_A SLPATOL_A",
      "Absolute delta convergence tolerance"
      "\n\nDefault: -1.0",
      XSLP_ATOL_A, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_ATOL_A

#ifdef XSLP_ATOL_R
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_atol_r XSLP_ATOL_R SLPATOL_R",
      "Relative delta convergence tolerance"
      "\n\nDefault: -1.0",
      XSLP_ATOL_R, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_ATOL_R

#ifdef XSLP_AUGMENTATION
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_augmentation XSLP_AUGMENTATION SLPAUGMENTATION",
      "Bit map describing the SLP augmentation method(s) to be used"
      "\n\n"
      "Values (default: 12 (sets bits 2 and 3)):\n"
      "\n- (0)  Minimum augmentation."
      "\n- (1)  Even handed augmentation."
      "\n- (2)  Penalty error vectors on all non-linear equality constraints."
      "\n- (3)  Penalty error vectors on all non-linear inequality constraints."
      "\n- (4)  Penalty vectors to exceed step bounds."
      "\n- (5)  Use arithmetic means to estimate penalty weights."
      "\n- (6)  Estimate step bounds from values of row coefficients."
      "\n- (7)  Estimate step bounds from absolute values of row coefficients."
      "\n- (8)  Row-based step bounds."
      "\n- (9)  Penalty error vectors on all constraints."
      "\n- (10)  Intial values do not imply an SLP variable."
      "\n- (12)  Avoid running an LP around fixed initial values trying to get feasible.",
      XSLP_AUGMENTATION, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_AUGMENTATION

#ifdef XSLP_AUTOSAVE
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_autosave XSLP_AUTOSAVE SLPAUTOSAVE",
      "Frequency with which to save the model"
      "\n\nDefault: 0",
      XSLP_AUTOSAVE, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_AUTOSAVE

#ifdef XSLP_BARCROSSOVERSTART
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_barcrossoverstart XSLP_BARCROSSOVERSTART SLPBARCROSSOVERSTART",
      "Default crossover activation behaviour for barrier start"
      "\n\nDefault: 0",
      XSLP_BARCROSSOVERSTART, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_BARCROSSOVERSTART

#ifdef XSLP_BARLIMIT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_barlimit XSLP_BARLIMIT SLPBARLIMIT",
      "Number of initial SLP iterations using the barrier method"
      "\n\nDefault: 0",
      XSLP_BARLIMIT, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_BARLIMIT

#ifdef XSLP_BARSTALLINGLIMIT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_barstallinglimit XSLP_BARSTALLINGLIMIT SLPBARSTALLINGLIMIT",
      "Number of iterations to allow numerical failures in barrier before switching to dual"
      "\n\nDefault: 3",
      XSLP_BARSTALLINGLIMIT, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_BARSTALLINGLIMIT

#ifdef XSLP_BARSTALLINGOBJLIMIT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_barstallingobjlimit XSLP_BARSTALLINGOBJLIMIT SLPBARSTALLINGOBJLIMIT",
      "Number of iterations over which to measure the objective change for barrier iterations with no crossover"
      "\n\nDefault: 3",
      XSLP_BARSTALLINGOBJLIMIT, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_BARSTALLINGOBJLIMIT

#ifdef XSLP_BARSTALLINGTOL
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_barstallingtol XSLP_BARSTALLINGTOL SLPBARSTALLINGTOL",
      "Required change in the objective when progress is measured in barrier iterations without crossover"
      "\n\nDefault: 0.05",
      XSLP_BARSTALLINGTOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_BARSTALLINGTOL

#ifdef XSLP_BARSTARTOPS
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_barstartops XSLP_BARSTARTOPS SLPBARSTARTOPS",
      "Controls behaviour when the barrier is used to solve the linearizations"
      "\n\n"
      "Values (default: -1):\n"
      "\n- (0)  Check objective progress when no crossover is applied."
      "\n- (1)  Fall back to dual simplex if too many numerical problems are reported by the barrier."
      "\n- (2)  If a non-vertex converged solution found by barrier without crossover can be returned as a final solution.",
      XSLP_BARSTARTOPS, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_BARSTARTOPS

#ifdef XSLP_BOUNDTHRESHOLD
    MPD( AddSolverOption_MergeDuplicates("pre:xslp_boundthreshold XSLP_BOUNDTHRESHOLD SLPBOUNDTHRESHOLD",
      "The maximum size of a bound that can be introduced by nonlinear presolve."
      "\n\nDefault: 1.0e+10",
      XSLP_BOUNDTHRESHOLD, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_BOUNDTHRESHOLD

#ifdef XSLP_CALCTHREADS
    MPD( AddSolverOption_MergeDuplicates("tech:xslp_calcthreads XSLP_CALCTHREADS NLPCALCTHREADS",
      "Number of threads used for formula and derivatives evaluations"
      "\n\nDefault: -1 (determined by XSLP_THREADS)",
      XSLP_CALCTHREADS, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_CALCTHREADS

#ifdef XSLP_CASCADE
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_cascade XSLP_CASCADE SLPCASCADE",
      "Bit map describing the cascading to be used"
      "\n\n"
      "Values (default: 257):\n"
      "\n- (0)  Apply cascading to all variables with determining rows."
      "\n- (1)  Apply cascading to SLP variables which appear in coefficients and which would change by more than XPRS_FEASTOL."
      "\n- (2)  Apply cascading to all SLP variables which appear in coefficients."
      "\n- (3)  Apply cascading to SLP variables which are structural and which would change by more than XPRS_FEASTOL."
      "\n- (4)  Apply cascading to all SLP variables which are structural."
      "\n- (5)  Create secondary order groupping DR rows with instantiated user functions together in the order."
      "\n- (6)  In cases where the determining column is below XSLP_DRCOLTOL, fix at the previous rather than current value."
      "\n- (7)  In cases where the determining column is below XSLP_DRCOLTOL, fix within a range XSLP_DRFIXRANGE of previous value."
      "\n- (8)  Automatically determine whether to apply cascading.",
      XSLP_CASCADE, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_CASCADE

#ifdef XSLP_CASCADENLIMIT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_cascadenlimit XSLP_CASCADENLIMIT SLPCASCADENLIMIT",
      "Maximum number of iterations for cascading with non-linear determining rows"
      "\n\nDefault: 10",
      XSLP_CASCADENLIMIT, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_CASCADENLIMIT

#ifdef XSLP_CASCADETOL_PA
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_cascadetol_pa XSLP_CASCADETOL_PA SLPCASCADETOL_PA",
      "Absolute cascading print tolerance"
      "\n\nDefault: 0.01",
      XSLP_CASCADETOL_PA, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_CASCADETOL_PA

#ifdef XSLP_CASCADETOL_PR
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_cascadetol_pr XSLP_CASCADETOL_PR SLPCASCADETOL_PR",
      "Relative cascading print tolerance"
      "\n\nDefault: 0.01",
      XSLP_CASCADETOL_PR, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_CASCADETOL_PR

#ifdef XSLP_CDTOL_A
    MPD( AddSolverOption_MergeDuplicates("tol:xslp_cdtol_a XSLP_CDTOL_A SLPCDTOL_A",
      "Absolute tolerance for deducing constant derivatives"
      "\n\nDefault: 1.0e-08",
      XSLP_CDTOL_A, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_CDTOL_A

#ifdef XSLP_CDTOL_R
    MPD( AddSolverOption_MergeDuplicates("tol:xslp_cdtol_r XSLP_CDTOL_R SLPCDTOL_R",
      "Relative tolerance for deducing constant derivatives"
      "\n\nDefault: 1.0e-08",
      XSLP_CDTOL_R, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_CDTOL_R

#ifdef XSLP_CLAMPSHRINK
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_clampshrink XSLP_CLAMPSHRINK SLPCLAMPSHRINK",
      "Shrink ratio used to impose strict convergence on variables converged in extended criteria only"
      "\n\nDefault: 0.3",
      XSLP_CLAMPSHRINK, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_CLAMPSHRINK

#ifdef XSLP_CLAMPVALIDATIONTOL_A
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_clampvalidationtol_a XSLP_CLAMPVALIDATIONTOL_A SLPCLAMPVALIDATIONTOL_A",
      "Absolute validation tolerance for applying XSLP_CLAMPSHRINK"
      "\n\nDefault: 1.0e-06",
      XSLP_CLAMPVALIDATIONTOL_A, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_CLAMPVALIDATIONTOL_A

#ifdef XSLP_CLAMPVALIDATIONTOL_R
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_clampvalidationtol_r XSLP_CLAMPVALIDATIONTOL_R SLPCLAMPVALIDATIONTOL_R",
      "Relative validation tolerance for applying XSLP_CLAMPSHRINK"
      "\n\nDefault: 1.0e-06",
      XSLP_CLAMPVALIDATIONTOL_R, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_CLAMPVALIDATIONTOL_R

#ifdef XSLP_CONTROL
    MPD( AddSolverOption_MergeDuplicates("bit:xslp_control XSLP_CONTROL",
      "Bit map describing which Xpress NonLinear functions also activate the corresponding Optimizer Library function"
      "\n\n"
      "Values (default: 0 (no bits set)):\n"
      "\n- (0)  Xpress NonLinear problem management functions do NOT invoke the corresponding Optimizer Library function for the underlying linear problem."
      "\n- (1)  XSLPcopycontrols does NOT invoke XPRScopycontrols."
      "\n- (2)  XSLPcopycallbacks does NOT invoke XPRScopycallbacks."
      "\n- (3)  XSLPcopyprob does NOT invoke XPRScopyprob."
      "\n- (4)  XSLPsetdefaults does NOT invoke XPRSsetdefaults."
      "\n- (5)  XSLPsave does NOT invoke XPRSsave."
      "\n- (6)  XSLPrestore does NOT invoke XPRSrestore.",
      XSLP_CONTROL, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_CONTROL

#ifdef XSLP_CONVERGENCEOPS
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_convergenceops XSLP_CONVERGENCEOPS SLPCONVERGENCEOPS",
      "Bit map describing which convergence tests should be carried out"
      "\n\n"
      "Values (default: 39935 (bits 0-9, 11-12, and 15 are set)):\n"
      "\n- (0)  Execute the closure tolerance checks."
      "\n- (1)  Execute the delta tolerance checks."
      "\n- (2)  Execute the matrix tolerance checks."
      "\n- (3)  Execute the impact tolerance checks."
      "\n- (4)  Execute the slack impact tolerance checks."
      "\n- (5)  Check for user provided convergence."
      "\n- (6)  Execute the objective range checks."
      "\n- (7)  Execute the objective range + constraint activity check."
      "\n- (8)  Execute the objective range + active step bound check."
      "\n- (9)  Execute the convergence continuation check."
      "\n- (10)  Take scaling of individual variables / rows into account."
      "\n- (11)  Execute the validation target convergence checks."
      "\n- (12)  Execute the first order optimality target convergence checks."
      "\n- (13)  Allow convex quadratic problems to converge on extended criteria."
      "\n- (15)  Do not declare converged if still sufficient improvement in objective.",
      XSLP_CONVERGENCEOPS, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_CONVERGENCEOPS

#ifdef XSLP_CTOL
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_ctol XSLP_CTOL SLPCTOL",
      "Closure convergence tolerance"
      "\n\nDefault: -1.0",
      XSLP_CTOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_CTOL

#ifdef XSLP_CUTSTRATEGY
    MPD( AddSolverOption_MergeDuplicates("mislp:xslp_cutstrategy XSLP_CUTSTRATEGY SLPCUTSTRATEGY",
      "Determines whihc cuts to apply in the MISLP search when the default SLP-in-MIP strategy is used."
      "\n\nDefault: 0",
      XSLP_CUTSTRATEGY, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_CUTSTRATEGY

#ifdef XSLP_DAMP
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_damp XSLP_DAMP SLPDAMP",
      "Damping factor for updating values of variables"
      "\n\nDefault: 1",
      XSLP_DAMP, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_DAMP

#ifdef XSLP_DAMPEXPAND
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_dampexpand XSLP_DAMPEXPAND SLPDAMPEXPAND",
      "Multiplier to increase damping factor during dynamic damping"
      "\n\nDefault: 1",
      XSLP_DAMPEXPAND, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_DAMPEXPAND

#ifdef XSLP_DAMPMAX
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_dampmax XSLP_DAMPMAX SLPDAMPMAX",
      "Maximum value for the damping factor of a variable during dynamic damping"
      "\n\nDefault: 1",
      XSLP_DAMPMAX, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_DAMPMAX

#ifdef XSLP_DAMPMIN
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_dampmin XSLP_DAMPMIN SLPDAMPMIN",
      "Minimum value for the damping factor of a variable during dynamic damping"
      "\n\nDefault: 1",
      XSLP_DAMPMIN, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_DAMPMIN

#ifdef XSLP_DAMPSHRINK
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_dampshrink XSLP_DAMPSHRINK SLPDAMPSHRINK",
      "Multiplier to decrease damping factor during dynamic damping"
      "\n\nDefault: 1",
      XSLP_DAMPSHRINK, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_DAMPSHRINK

#ifdef XSLP_DAMPSTART
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_dampstart XSLP_DAMPSTART SLPDAMPSTART",
      "SLP iteration at which damping is activated"
      "\n\nDefault: 0",
      XSLP_DAMPSTART, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_DAMPSTART

#ifdef XSLP_DEFAULTIV
    MPD( AddSolverOption_MergeDuplicates("dat:xslp_defaultiv XSLP_DEFAULTIV NLPDEFAULTIV",
      "Default initial value for an SLP variable if none is explicitly given"
      "\n\nDefault: 100",
      XSLP_DEFAULTIV, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_DEFAULTIV

#ifdef XSLP_DEFAULTSTEPBOUND
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_defaultstepbound XSLP_DEFAULTSTEPBOUND SLPDEFAULTSTEPBOUND",
      "Minimum initial value for the step bound of an SLP variable if none is explicitly given"
      "\n\nDefault: 16",
      XSLP_DEFAULTSTEPBOUND, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_DEFAULTSTEPBOUND

#ifdef XSLP_DELAYUPDATEROWS
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_delayupdaterows XSLP_DELAYUPDATEROWS SLPDELAYUPDATEROWS",
      "Number of SLP iterations before update rows are fully activated"
      "\n\nDefault: 2",
      XSLP_DELAYUPDATEROWS, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_DELAYUPDATEROWS

#ifdef XSLP_DELTACOST
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_deltacost XSLP_DELTACOST SLPDELTACOST",
      "Initial penalty cost multiplier for penalty delta vectors"
      "\n\nDefault: 200",
      XSLP_DELTACOST, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_DELTACOST

#ifdef XSLP_DELTACOSTFACTOR
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_deltacostfactor XSLP_DELTACOSTFACTOR SLPDELTACOSTFACTOR",
      "Factor for increasing cost multiplier on total penalty delta vectors"
      "\n\nDefault: 1.3",
      XSLP_DELTACOSTFACTOR, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_DELTACOSTFACTOR

#ifdef XSLP_DELTAFORMAT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_deltaformat XSLP_DELTAFORMAT SLPDELTAFORMAT",
      "Formatting string for creation of names for SLP delta vectors"
      "\n\nDefault: pD_%s where p is a unique prefix for names in the current problem",
      XSLP_DELTAFORMAT) );
#endif  // ifdef XSLP_DELTAFORMAT

#ifdef XSLP_DELTAMAXCOST
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_deltamaxcost XSLP_DELTAMAXCOST SLPDELTAMAXCOST",
      "Maximum penalty cost multiplier for penalty delta vectors"
      "\n\nDefault: XPRS_PLUSINFINITY",
      XSLP_DELTAMAXCOST, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_DELTAMAXCOST

#ifdef XSLP_DELTAOFFSET
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_deltaoffset XSLP_DELTAOFFSET SLPDELTAOFFSET",
      "Position of first character of SLP variable name used to create name of delta vector"
      "\n\nDefault: 0",
      XSLP_DELTAOFFSET, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_DELTAOFFSET

#ifdef XSLP_DELTAZLIMIT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_deltazlimit XSLP_DELTAZLIMIT SLPDELTAZLIMIT",
      "Number of SLP iterations during which to apply XSLP_DELTA_Z"
      "\n\nDefault: 0",
      XSLP_DELTAZLIMIT, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_DELTAZLIMIT

#ifdef XSLP_DELTA_A
    MPD( AddSolverOption_MergeDuplicates("diff:xslp_delta_a XSLP_DELTA_A SLPDELTA_A",
      "Absolute perturbation of values for calculating numerical derivatives"
      "\n\nDefault: 0.001",
      XSLP_DELTA_A, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_DELTA_A

#ifdef XSLP_DELTA_INFINITY
    MPD( AddSolverOption_MergeDuplicates("diff:xslp_delta_infinity XSLP_DELTA_INFINITY SLPDELTA_INFINITY",
      "Maximum value for partial derivatives"
      "\n\nDefault: 1.0e+15",
      XSLP_DELTA_INFINITY, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_DELTA_INFINITY

#ifdef XSLP_DELTA_R
    MPD( AddSolverOption_MergeDuplicates("diff:xslp_delta_r XSLP_DELTA_R SLPDELTA_R",
      "Relative perturbation of values for calculating numerical derivatives"
      "\n\nDefault: 0.001",
      XSLP_DELTA_R, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_DELTA_R

#ifdef XSLP_DELTA_X
    MPD( AddSolverOption_MergeDuplicates("diff:xslp_delta_x XSLP_DELTA_X SLPDELTA_X",
      "Minimum absolute value of delta coefficients to be retained"
      "\n\nDefault: 1.0e-6",
      XSLP_DELTA_X, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_DELTA_X

#ifdef XSLP_DELTA_Z
    MPD( AddSolverOption_MergeDuplicates("diff:xslp_delta_z XSLP_DELTA_Z SLPDELTA_Z",
      "Tolerance used when calculating derivatives"
      "\n\nDefault: 0.00001",
      XSLP_DELTA_Z, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_DELTA_Z

#ifdef XSLP_DELTA_ZERO
    MPD( AddSolverOption_MergeDuplicates("diff:xslp_delta_zero XSLP_DELTA_ZERO SLPDELTA_ZERO",
      "Absolute zero acceptance tolerance used when calculating derivatives"
      "\n\nDefault: -1.0 (not applied)",
      XSLP_DELTA_ZERO, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_DELTA_ZERO

#ifdef XSLP_DERIVATIVES
    MPD( AddSolverOption_MergeDuplicates("diff:xslp_derivatives XSLP_DERIVATIVES NLPDERIVATIVES",
      "Bitmap describing the method of calculating derivatives"
      "\n\n"
      "Values (default: 1):\n"
      "\n- (0)  analytic derivatives where possible"
      "\n- (1)  avoid embedding numerical derivatives of instantiated functions into analytic derivatives",
      XSLP_DERIVATIVES, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_DERIVATIVES

#ifdef XSLP_DETERMINISTIC
    MPD( AddSolverOption_MergeDuplicates("mislp:xslp_deterministic XSLP_DETERMINISTIC NLPDETERMINISTIC",
      "Determines if the parallel features of SLP should be guaranteed to be deterministic"
      "\n\nDefault: 1",
      XSLP_DETERMINISTIC, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_DETERMINISTIC

#ifdef XSLP_DJTOL
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_djtol XSLP_DJTOL SLPDJTOL",
      "Tolerance on DJ value for determining if a variable is at its step bound"
      "\n\nDefault: 1.0e-6",
      XSLP_DJTOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_DJTOL

#ifdef XSLP_DRCOLDJTOL
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_drcoldjtol XSLP_DRCOLDJTOL SLPDRCOLDJTOL",
      "Reduced cost tolerance on the delta variable when fixing due to the determining column being below XSLP_DRCOLTOL."
      "\n\nDefault: 0.0",
      XSLP_DRCOLDJTOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_DRCOLDJTOL

#ifdef XSLP_DRCOLTOL
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_drcoltol XSLP_DRCOLTOL SLPDRCOLTOL",
      "The minimum absolute magnitude of a determining column, for which the determined variable is still regarded as well defined"
      "\n\nDefault: 1.0e-6",
      XSLP_DRCOLTOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_DRCOLTOL

#ifdef XSLP_DRFIXRANGE
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_drfixrange XSLP_DRFIXRANGE SLPDRFIXRANGE",
      "The range around the previous value where variables are fixed in cascading if the determining column is below XSLP_DRCOLTOL."
      "\n\nDefault: 0.1",
      XSLP_DRFIXRANGE, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_DRFIXRANGE

#ifdef XSLP_ECFCHECK
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_ecfcheck XSLP_ECFCHECK SLPECFCHECK",
      "Check feasibility at the point of linearization for extended convergence criteria"
      "\n\n"
      "Values (default: 1):\n"
      "\n- (0)  no check (extended criteria are always used);"
      "\n- (1)  check until one infeasible constraint is found;"
      "\n- (2)  check all constraints.",
      XSLP_ECFCHECK, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_ECFCHECK

#ifdef XSLP_ECFTOL_A
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_ecftol_a XSLP_ECFTOL_A SLPECFTOL_A",
      "Absolute tolerance on testing feasibility at the point of linearization"
      "\n\nDefault: -1.0",
      XSLP_ECFTOL_A, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_ECFTOL_A

#ifdef XSLP_ECFTOL_R
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_ecftol_r XSLP_ECFTOL_R SLPECFTOL_R",
      "Relative tolerance on testing feasibility at the point of linearization"
      "\n\nDefault: -1.0",
      XSLP_ECFTOL_R, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_ECFTOL_R

#ifdef XSLP_ECHOXPRSMESSAGES
    MPD( AddSolverOption_MergeDuplicates("log:xslp_echoxprsmessages XSLP_ECHOXPRSMESSAGES",
      "Controls if the XSLP message callback should relay messages from the XPRS library."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  automatic: if an XSLP message callback is not set, then messages from the nonlinear solver are sent to the XPRS message callback; if an XSLP message callback is set, then messages are not echoed."
      "\n- (0)  the XPRS and XSLP message callbacks are treated as independent."
      "\n- (1)  messages from the XPRS message callback are sent to the XSLP message callback."
      "\n- (2)  messages from the nonlinear solver are sent to the XPRS message callback.",
      XSLP_ECHOXPRSMESSAGES, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_ECHOXPRSMESSAGES

#ifdef XSLP_ENFORCECOSTSHRINK
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_enforcecostshrink XSLP_ENFORCECOSTSHRINK SLPENFORCECOSTSHRINK",
      "Factor by which to decrease the current penalty multiplier when enforcing rows. "
      "\n\nDefault: 0.00001",
      XSLP_ENFORCECOSTSHRINK, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_ENFORCECOSTSHRINK

#ifdef XSLP_ENFORCEMAXCOST
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_enforcemaxcost XSLP_ENFORCEMAXCOST SLPENFORCEMAXCOST",
      "Maximum penalty cost in the objective before enforcing most violating rows"
      "\n\nDefault: 1.0e+11",
      XSLP_ENFORCEMAXCOST, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_ENFORCEMAXCOST

#ifdef XSLP_ERRORCOST
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_errorcost XSLP_ERRORCOST SLPERRORCOST",
      "Initial penalty cost multiplier for penalty error vectors"
      "\n\nDefault: 200",
      XSLP_ERRORCOST, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_ERRORCOST

#ifdef XSLP_ERRORCOSTFACTOR
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_errorcostfactor XSLP_ERRORCOSTFACTOR SLPERRORCOSTFACTOR",
      "Factor for increasing cost multiplier on total penalty error vectors"
      "\n\nDefault: 1.3",
      XSLP_ERRORCOSTFACTOR, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_ERRORCOSTFACTOR

#ifdef XSLP_ERRORMAXCOST
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_errormaxcost XSLP_ERRORMAXCOST SLPERRORMAXCOST",
      "Maximum penalty cost multiplier for penalty error vectors"
      "\n\nDefault: XPRS_PLUSINFINITY",
      XSLP_ERRORMAXCOST, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_ERRORMAXCOST

#ifdef XSLP_ERROROFFSET
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_erroroffset XSLP_ERROROFFSET SLPERROROFFSET",
      "Position of first character of constraint name used to create name of penalty error vectors"
      "\n\nDefault: 0",
      XSLP_ERROROFFSET, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_ERROROFFSET

#ifdef XSLP_ERRORTOL_A
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_errortol_a XSLP_ERRORTOL_A SLPERRORTOL_A",
      "Absolute tolerance for error vectors"
      "\n\nDefault: 0.00001",
      XSLP_ERRORTOL_A, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_ERRORTOL_A

#ifdef XSLP_ERRORTOL_P
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_errortol_p XSLP_ERRORTOL_P SLPERRORTOL_P",
      "Absolute tolerance for printing error vectors"
      "\n\nDefault: 0.0001",
      XSLP_ERRORTOL_P, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_ERRORTOL_P

#ifdef XSLP_ESCALATION
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_escalation XSLP_ESCALATION SLPESCALATION",
      "Factor for increasing cost multiplier on individual penalty error vectors"
      "\n\nDefault: 1.25",
      XSLP_ESCALATION, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_ESCALATION

#ifdef XSLP_ETOL_A
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_etol_a XSLP_ETOL_A SLPETOL_A",
      "Absolute tolerance on penalty vectors"
      "\n\nDefault: 0.0001",
      XSLP_ETOL_A, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_ETOL_A

#ifdef XSLP_ETOL_R
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_etol_r XSLP_ETOL_R SLPETOL_R",
      "Relative tolerance on penalty vectors"
      "\n\nDefault: 0.0001",
      XSLP_ETOL_R, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_ETOL_R

#ifdef XSLP_EVALUATE
    MPD( AddSolverOption_MergeDuplicates("func:xslp_evaluate XSLP_EVALUATE NLPEVALUATE",
      "Evaluation strategy for user functions"
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  use derivatives where possible;"
      "\n- (1)  always re-evaluate.",
      XSLP_EVALUATE, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_EVALUATE

#ifdef XSLP_EVTOL_A
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_evtol_a XSLP_EVTOL_A SLPEVTOL_A",
      "Absolute tolerance on total penalty costs"
      "\n\nDefault: -1.0",
      XSLP_EVTOL_A, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_EVTOL_A

#ifdef XSLP_EVTOL_R
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_evtol_r XSLP_EVTOL_R SLPEVTOL_R",
      "Relative tolerance on total penalty costs"
      "\n\nDefault: -1.0",
      XSLP_EVTOL_R, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_EVTOL_R

#ifdef XSLP_EXPAND
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_expand XSLP_EXPAND SLPEXPAND",
      "Multiplier to increase a step bound"
      "\n\nDefault: 2",
      XSLP_EXPAND, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_EXPAND

#ifdef XSLP_FEASTOLTARGET
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_feastoltarget XSLP_FEASTOLTARGET SLPFEASTOLTARGET",
      "When set, this defines a target feasibility tolerance to which the linearizations are solved to"
      "\n\nDefault: 0 (ignored, not set)",
      XSLP_FEASTOLTARGET, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_FEASTOLTARGET

#ifdef XSLP_FILTER
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_filter XSLP_FILTER SLPFILTER",
      "Bit map for controlling solution updates"
      "\n\n"
      "Values (default: 3 (bit 0,1)):\n"
      "\n- (0)  retain best solution according to the merit function."
      "\n- (1)  check cascaded solutions against improvements in the merit function."
      "\n- (2)  force minimum step sizes in line search."
      "\n- (3)  accept the trust region step is the line search returns a zero step size.",
      XSLP_FILTER, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_FILTER

#ifdef XSLP_FINDIV
    MPD( AddSolverOption_MergeDuplicates("heur:xslp_findiv XSLP_FINDIV NLPFINDIV",
      "Option for running a heuristic to find a feasible initial point"
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatic (default)."
      "\n- (0)  Disable the heuristic."
      "\n- (1)  Enable the heuristic.",
      XSLP_FINDIV, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_FINDIV

#ifdef XSLP_FUNCEVAL
    MPD( AddSolverOption_MergeDuplicates("func:xslp_funceval XSLP_FUNCEVAL NLPFUNCEVAL",
      "Bit map for determining the method of evaluating user functions and their derivatives"
      "\n\n"
      "Values (default: 0):\n"
      "\n- (3)  evaluate function whenever independent variables change."
      "\n- (4)  evaluate function when independent variables change outside tolerances."
      "\n- (5)  application of bits 3-4:  0 = functions which do not have a defined re-evaluation mode;1 = all functions."
      "\n- (6)  tangential derivatives."
      "\n- (7)  forward derivatives"
      "\n- (8)  application of bits 6-7:  0 = functions which do not have a defined derivative mode;1 = all functions.",
      XSLP_FUNCEVAL, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_FUNCEVAL

#ifdef XSLP_GRANULARITY
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_granularity XSLP_GRANULARITY SLPGRANULARITY",
      "Base for calculating penalty costs"
      "\n\nDefault: 4",
      XSLP_GRANULARITY, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_GRANULARITY

#ifdef XSLP_GRIDHEURSELECT
    MPD( AddSolverOption_MergeDuplicates("heur:xslp_gridheurselect XSLP_GRIDHEURSELECT SLPGRIDHEURSELECT",
      "Bit map selectin which heuristics to run if the problem has variable with an integer delta"
      "\n\n"
      "Values (default: 6):\n"
      "\n- (0)  Enumeration: try all combinations."
      "\n- (1)  Simple search heuristics."
      "\n- (2)  Simulated annealing.",
      XSLP_GRIDHEURSELECT, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_GRIDHEURSELECT

#ifdef XSLP_HESSIAN
    MPD( AddSolverOption_MergeDuplicates("diff:xslp_hessian XSLP_HESSIAN NLPHESSIAN",
      "Second order differentiation mode when using analytical derivatives "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1,0)  automatic selection"
      "\n- (1)  numerical derivatives (finite difference)"
      "\n- (2)  symbolic differentiation"
      "\n- (3)  automatic differentiation",
      XSLP_HESSIAN, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_HESSIAN

#ifdef XSLP_HEURSTRATEGY
    MPD( AddSolverOption_MergeDuplicates("heur:xslp_heurstrategy XSLP_HEURSTRATEGY SLPHEURSTRATEGY",
      "Branch and Bound: This specifies the MINLP heuristic strategy. On some problems it is worth trying more comprehensive heuristic strategies by setting HEURSTRATEGY to 2 or 3."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatic selection of heuristic strategy (depending on XPRS_HEUREMPHASIS)."
      "\n- (0)  No heuristics."
      "\n- (1)  Basic heuristic strategy."
      "\n- (2)  Enhanced heuristic strategy."
      "\n- (3)  Extensive heuristic strategy."
      "\n- (4)  Run all heuristics without effort limits.",
      XSLP_HEURSTRATEGY, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_HEURSTRATEGY

#ifdef XSLP_INFEASLIMIT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_infeaslimit XSLP_INFEASLIMIT SLPINFEASLIMIT",
      "The maximum number of consecutive infeasible SLP iterations which can occur before Xpress-SLP terminates"
      "\n\nDefault: 3",
      XSLP_INFEASLIMIT, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_INFEASLIMIT

#ifdef XSLP_INFINITY
    MPD( AddSolverOption_MergeDuplicates("num:xslp_infinity XSLP_INFINITY NLPINFINITY",
      "Value returned by a divide-by-zero in a formula"
      "\n\nDefault: 1.0e+10",
      XSLP_INFINITY, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_INFINITY

#ifdef XSLP_ITERFALLBACKOPS
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_iterfallbackops XSLP_ITERFALLBACKOPS SLPITERFALLBACKOPS",
      "Alternative LP level control values for numerically challengeing problems"
      "\n\nDefault: none",
      XSLP_ITERFALLBACKOPS) );
#endif  // ifdef XSLP_ITERFALLBACKOPS

#ifdef XSLP_ITERLIMIT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_iterlimit XSLP_ITERLIMIT SLPITERLIMIT",
      "The maximum number of SLP iterations"
      "\n\nDefault: 1000",
      XSLP_ITERLIMIT, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_ITERLIMIT

#ifdef XSLP_ITOL_A
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_itol_a XSLP_ITOL_A SLPITOL_A",
      "Absolute impact convergence tolerance"
      "\n\nDefault: -1.0",
      XSLP_ITOL_A, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_ITOL_A

#ifdef XSLP_ITOL_R
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_itol_r XSLP_ITOL_R SLPITOL_R",
      "Relative impact convergence tolerance"
      "\n\nDefault: -1.0",
      XSLP_ITOL_R, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_ITOL_R

#ifdef XSLP_IVNAME
    MPD( AddSolverOption_MergeDuplicates("tech:xslp_ivname XSLP_IVNAME NLPIVNAME",
      "Name of the set of initial values to be used"
      "\n\nDefault: none",
      XSLP_IVNAME) );
#endif  // ifdef XSLP_IVNAME

#ifdef XSLP_JACOBIAN
    MPD( AddSolverOption_MergeDuplicates("diff:xslp_jacobian XSLP_JACOBIAN NLPJACOBIAN",
      "First order differentiation mode when using analytical derivatives"
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1,0)  automatic selection"
      "\n- (1)  numerical derivatives (finite difference)"
      "\n- (2)  symbolic differentiation"
      "\n- (3)  automatic differentiation",
      XSLP_JACOBIAN, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_JACOBIAN

#ifdef XSLP_KEEPEQUALSCOLUMN
    MPD( AddSolverOption_MergeDuplicates("tech:xslp_keepequalscolumn XSLP_KEEPEQUALSCOLUMN NLPKEEPEQUALSCOLUMN",
      "When set to a nonzero value, the MPS reader will keep the equals column in the problem"
      "\n\nDefault: 0",
      XSLP_KEEPEQUALSCOLUMN, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_KEEPEQUALSCOLUMN

#ifdef XSLP_LINQUADBR
    MPD( AddSolverOption_MergeDuplicates("pre:xslp_linquadbr XSLP_LINQUADBR NLPLINQUADBR",
      "Use linear and quadratic constraints and objective function to further reduce bounds on all variables"
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  automatic selection"
      "\n- (0)  disable"
      "\n- (1)  enable",
      XSLP_LINQUADBR, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_LINQUADBR

#ifdef XSLP_LOG
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_log XSLP_LOG NLPLOG",
      "Level of printing during SLP iterations"
      "\n\n"
      "Values (default: 0):\n"
      "\n- (-1)  none"
      "\n- (0)  minimal"
      "\n- (1)  normal: iteration, penalty vectors"
      "\n- (2)  omit from convergence log any variables which have converged"
      "\n- (3)  omit from convergence log any variables which have already converged (except variables on step bounds)"
      "\n- (4)  include all variables in convergence log"
      "\n- (5)  include user function call communications in the log",
      XSLP_LOG, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_LOG

#ifdef XSLP_LSITERLIMIT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_lsiterlimit XSLP_LSITERLIMIT SLPLSITERLIMIT",
      "Number of iterations in the line search"
      "\n\nDefault: 0",
      XSLP_LSITERLIMIT, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_LSITERLIMIT

#ifdef XSLP_LSPATTERNLIMIT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_lspatternlimit XSLP_LSPATTERNLIMIT SLPLSPATTERNLIMIT",
      "Number of iterations in the pattern search preceding the line search"
      "\n\nDefault: 0",
      XSLP_LSPATTERNLIMIT, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_LSPATTERNLIMIT

#ifdef XSLP_LSSTART
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_lsstart XSLP_LSSTART SLPLSSTART",
      "Iteration in which to active the line search"
      "\n\nDefault: 8",
      XSLP_LSSTART, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_LSSTART

#ifdef XSLP_LSZEROLIMIT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_lszerolimit XSLP_LSZEROLIMIT SLPLSZEROLIMIT",
      "Maximum number of zero length line search steps before line search is deactivated"
      "\n\nDefault: 5",
      XSLP_LSZEROLIMIT, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_LSZEROLIMIT

#ifdef XSLP_MATRIXTOL
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_matrixtol XSLP_MATRIXTOL SLPMATRIXTOL",
      " Nonzero tolerance for dropping coefficients from the linearization."
      "\n\nDefault: 0.0",
      XSLP_MATRIXTOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_MATRIXTOL

#ifdef XSLP_MAXTIME
    MPD( AddSolverOption_MergeDuplicates("lim:xslp_maxtime XSLP_MAXTIME NLPMAXTIME",
      "The maximum time in seconds that the SLP optimization will run before it terminates"
      "\n\nDefault: 0",
      XSLP_MAXTIME, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_MAXTIME

#ifdef XSLP_MAXWEIGHT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_maxweight XSLP_MAXWEIGHT SLPMAXWEIGHT",
      "Maximum penalty weight for delta or error vectors"
      "\n\nDefault: 100",
      XSLP_MAXWEIGHT, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_MAXWEIGHT

#ifdef XSLP_MEMORYFACTOR
    MPD( AddSolverOption_MergeDuplicates("tech:xslp_memoryfactor XSLP_MEMORYFACTOR",
      "Factor for expanding size of dynamic arrays in memory"
      "\n\nDefault: 1.6",
      XSLP_MEMORYFACTOR, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_MEMORYFACTOR

#ifdef XSLP_MERITLAMBDA
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_meritlambda XSLP_MERITLAMBDA NLPMERITLAMBDA",
      "Factor by which the net objective is taken into account in the merit function"
      "\n\nDefault: 0.0",
      XSLP_MERITLAMBDA, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_MERITLAMBDA

#ifdef XSLP_MINSBFACTOR
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_minsbfactor XSLP_MINSBFACTOR SLPMINSBFACTOR",
      "Factor by which step bounds can be decreased beneath XSLP_ATOL_A"
      "\n\nDefault: 1.0",
      XSLP_MINSBFACTOR, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_MINSBFACTOR

#ifdef XSLP_MINUSDELTAFORMAT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_minusdeltaformat XSLP_MINUSDELTAFORMAT SLPMINUSDELTAFORMAT",
      "Formatting string for creation of names for SLP negative penalty delta vectors"
      "\n\nDefault: pD-%s where p is a unique prefix for names in the current problem",
      XSLP_MINUSDELTAFORMAT) );
#endif  // ifdef XSLP_MINUSDELTAFORMAT

#ifdef XSLP_MINUSERRORFORMAT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_minuserrorformat XSLP_MINUSERRORFORMAT SLPMINUSERRORFORMAT",
      "Formatting string for creation of names for SLP negative penalty error vectors"
      "\n\nDefault: pE-%s where p is a unique prefix for names in the current problem",
      XSLP_MINUSERRORFORMAT) );
#endif  // ifdef XSLP_MINUSERRORFORMAT

#ifdef XSLP_MINWEIGHT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_minweight XSLP_MINWEIGHT SLPMINWEIGHT",
      "Minimum penalty weight for delta or error vectors"
      "\n\nDefault: 0.01",
      XSLP_MINWEIGHT, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_MINWEIGHT

#ifdef XSLP_MIPALGORITHM
    MPD( AddSolverOption_MergeDuplicates("mislp:xslp_mipalgorithm XSLP_MIPALGORITHM SLPMIPALGORITHM",
      "Bitmap describing the MISLP algorithms to be used"
      "\n\n"
      "Values (default: 17 (bits 0 and4 are set)):\n"
      "\n- (0)  Solve initial SLP to convergence."
      "\n- (2)  Relax step bounds according to XSLP_MIPRELAXSTEPBOUNDS after initial node."
      "\n- (3)  Fix step bounds according to XSLP_MIPFIXSTEPBOUNDS after initial node."
      "\n- (4)  Relax step bounds according to XSLP_MIPRELAXSTEPBOUNDS at each node."
      "\n- (5)  Fix step bounds according to XSLP_MIPFIXSTEPBOUNDS at each node."
      "\n- (6)  Limit iterations at each node to XSLP_MIPITERLIMIT."
      "\n- (7)  Relax step bounds according to XSLP_MIPRELAXSTEPBOUNDS after MIP solution is found."
      "\n- (8)  Fix step bounds according to XSLP_MIPFIXSTEPBOUNDS after MIP solution is found."
      "\n- (9)  Use MIP at each SLP iteration instead of SLP at each node."
      "\n- (10)  Use MIP on converged SLP solution and then SLP on the resulting MIP solution.",
      XSLP_MIPALGORITHM, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_MIPALGORITHM

#ifdef XSLP_MIPCUTOFFCOUNT
    MPD( AddSolverOption_MergeDuplicates("mislp:xslp_mipcutoffcount XSLP_MIPCUTOFFCOUNT SLPMIPCUTOFFCOUNT",
      "Number of SLP iterations to check when considering a node for cutting off"
      "\n\nDefault: 5",
      XSLP_MIPCUTOFFCOUNT, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_MIPCUTOFFCOUNT

#ifdef XSLP_MIPCUTOFFLIMIT
    MPD( AddSolverOption_MergeDuplicates("mislp:xslp_mipcutofflimit XSLP_MIPCUTOFFLIMIT SLPMIPCUTOFFLIMIT",
      "Number of SLP iterations to check when considering a node for cutting off"
      "\n\nDefault: 10",
      XSLP_MIPCUTOFFLIMIT, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_MIPCUTOFFLIMIT

#ifdef XSLP_MIPCUTOFF_A
    MPD( AddSolverOption_MergeDuplicates("mislp:xslp_mipcutoff_a XSLP_MIPCUTOFF_A SLPMIPCUTOFF_A",
      "Absolute objective function cutoff for MIP termination"
      "\n\nDefault: 0.00001",
      XSLP_MIPCUTOFF_A, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_MIPCUTOFF_A

#ifdef XSLP_MIPCUTOFF_R
    MPD( AddSolverOption_MergeDuplicates("mislp:xslp_mipcutoff_r XSLP_MIPCUTOFF_R SLPMIPCUTOFF_R",
      "Absolute objective function cutoff for MIP termination"
      "\n\nDefault: 0.00001",
      XSLP_MIPCUTOFF_R, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_MIPCUTOFF_R

#ifdef XSLP_MIPDEFAULTALGORITHM
    MPD( AddSolverOption_MergeDuplicates("mislp:xslp_mipdefaultalgorithm XSLP_MIPDEFAULTALGORITHM SLPMIPDEFAULTALGORITHM",
      "Default algorithm to be used during the tree search in MISLP"
      "\n\nDefault: 3",
      XSLP_MIPDEFAULTALGORITHM, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_MIPDEFAULTALGORITHM

#ifdef XSLP_MIPERRORTOL_A
    MPD( AddSolverOption_MergeDuplicates("mislp:xslp_miperrortol_a XSLP_MIPERRORTOL_A SLPMIPERRORTOL_A",
      "Absolute penalty error cost tolerance for MIP cut-off"
      "\n\nDefault: 0 (inactive)",
      XSLP_MIPERRORTOL_A, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_MIPERRORTOL_A

#ifdef XSLP_MIPERRORTOL_R
    MPD( AddSolverOption_MergeDuplicates("mislp:xslp_miperrortol_r XSLP_MIPERRORTOL_R SLPMIPERRORTOL_R",
      "Relative penalty error cost tolerance for MIP cut-off"
      "\n\nDefault: 0 (inactive)",
      XSLP_MIPERRORTOL_R, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_MIPERRORTOL_R

#ifdef XSLP_MIPFIXSTEPBOUNDS
    MPD( AddSolverOption_MergeDuplicates("mislp:xslp_mipfixstepbounds XSLP_MIPFIXSTEPBOUNDS SLPMIPFIXSTEPBOUNDS",
      "Bitmap describing the step-bound fixing strategy during MISLP"
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  Fix step bounds on structural SLP variables which are not in coefficients."
      "\n- (1)  Fix step bounds on all structural SLP variables."
      "\n- (2)  Fix step bounds on SLP variables appearing only in coefficients."
      "\n- (3)  Fix step bounds on SLP variables appearing in coefficients.",
      XSLP_MIPFIXSTEPBOUNDS, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_MIPFIXSTEPBOUNDS

#ifdef XSLP_MIPITERLIMIT
    MPD( AddSolverOption_MergeDuplicates("mislp:xslp_mipiterlimit XSLP_MIPITERLIMIT SLPMIPITERLIMIT",
      "Maximum number of SLP iterations at each node"
      "\n\nDefault: 0",
      XSLP_MIPITERLIMIT, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_MIPITERLIMIT

#ifdef XSLP_MIPLOG
    MPD( AddSolverOption_MergeDuplicates("mislp:xslp_miplog XSLP_MIPLOG SLPMIPLOG",
      "Frequency with which MIP status is printed"
      "\n\nDefault: 0 (deterministic logging)",
      XSLP_MIPLOG, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_MIPLOG

#ifdef XSLP_MIPOCOUNT
    MPD( AddSolverOption_MergeDuplicates("mislp:xslp_mipocount XSLP_MIPOCOUNT SLPMIPOCOUNT",
      "Number of SLP iterations at each node over which to measure objective function variation"
      "\n\nDefault: 5",
      XSLP_MIPOCOUNT, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_MIPOCOUNT

#ifdef XSLP_MIPOTOL_A
    MPD( AddSolverOption_MergeDuplicates("mislp:xslp_mipotol_a XSLP_MIPOTOL_A SLPMIPOTOL_A",
      "Absolute objective function tolerance for MIP termination"
      "\n\nDefault: 0.00001",
      XSLP_MIPOTOL_A, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_MIPOTOL_A

#ifdef XSLP_MIPOTOL_R
    MPD( AddSolverOption_MergeDuplicates("mislp:xslp_mipotol_r XSLP_MIPOTOL_R SLPMIPOTOL_R",
      "Relative objective function tolerance for MIP termination"
      "\n\nDefault: 0.00001",
      XSLP_MIPOTOL_R, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_MIPOTOL_R

#ifdef XSLP_MIPRELAXSTEPBOUNDS
    MPD( AddSolverOption_MergeDuplicates("mislp:xslp_miprelaxstepbounds XSLP_MIPRELAXSTEPBOUNDS SLPMIPRELAXSTEPBOUNDS",
      "Bitmap describing the step-bound relaxation strategy during MISLP"
      "\n\n"
      "Values (default: 15 (relax all types)):\n"
      "\n- (0)  Relax step bounds on structural SLP variables which are not in coefficients."
      "\n- (1)  Relax step bounds on all structural SLP variables."
      "\n- (2)  Relax step bounds on SLP variables appearing only in coefficients."
      "\n- (3)  Relax step bounds on SLP variables appearing in coefficients.",
      XSLP_MIPRELAXSTEPBOUNDS, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_MIPRELAXSTEPBOUNDS

#ifdef XSLP_MSMAXBOUNDRANGE
    MPD( AddSolverOption_MergeDuplicates("alg:xslp_msmaxboundrange XSLP_MSMAXBOUNDRANGE MSMAXBOUNDRANGE",
      "Defines the maximum range inside which initial points are generated by multistart presets"
      "\n\nDefault: 1000",
      XSLP_MSMAXBOUNDRANGE, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_MSMAXBOUNDRANGE

#ifdef XSLP_MTOL_A
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_mtol_a XSLP_MTOL_A SLPMTOL_A",
      "Absolute effective matrix element convergence tolerance"
      "\n\nDefault: -1.0",
      XSLP_MTOL_A, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_MTOL_A

#ifdef XSLP_MTOL_R
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_mtol_r XSLP_MTOL_R SLPMTOL_R",
      "Relative effective matrix element convergence tolerance"
      "\n\nDefault: -1.0",
      XSLP_MTOL_R, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_MTOL_R

#ifdef XSLP_MULTISTART
    MPD( AddSolverOption_MergeDuplicates("alg:xslp_multistart XSLP_MULTISTART MULTISTART",
      "The multistart main control. Defines if the multistart search is to be initiated, or if only the baseline model is to be solved. "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Depends on if any multistart jobs have been added."
      "\n- (0)  Multistart is off."
      "\n- (1)  Multistart is on.",
      XSLP_MULTISTART, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_MULTISTART

#ifdef XSLP_MULTISTART_LOG
    MPD( AddSolverOption_MergeDuplicates("alg:xslp_multistart_log XSLP_MULTISTART_LOG MULTISTART_LOG",
      "The level of logging during the multistart run. "
      "\n\nDefault: 0",
      XSLP_MULTISTART_LOG, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_MULTISTART_LOG

#ifdef XSLP_MULTISTART_MAXSOLVES
    MPD( AddSolverOption_MergeDuplicates("alg:xslp_multistart_maxsolves XSLP_MULTISTART_MAXSOLVES MULTISTART_MAXSOLVES",
      "The maximum number of jobs to create during the multistart search. "
      "\n\nDefault: -1 (no upper limit)",
      XSLP_MULTISTART_MAXSOLVES, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_MULTISTART_MAXSOLVES

#ifdef XSLP_MULTISTART_MAXTIME
    MPD( AddSolverOption_MergeDuplicates("alg:xslp_multistart_maxtime XSLP_MULTISTART_MAXTIME MULTISTART_MAXTIME",
      "The maximum total time to be spent in the mutlistart search. "
      "\n\nDefault: 0 (no upper limit)",
      XSLP_MULTISTART_MAXTIME, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_MULTISTART_MAXTIME

#ifdef XSLP_MULTISTART_POOLSIZE
    MPD( AddSolverOption_MergeDuplicates("alg:xslp_multistart_poolsize XSLP_MULTISTART_POOLSIZE MULTISTART_POOLSIZE",
      "The maximum number of problem objects allowed to pool up before synchronization in the deterministic multistart. "
      "\n\nDefault: 2",
      XSLP_MULTISTART_POOLSIZE, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_MULTISTART_POOLSIZE

#ifdef XSLP_MULTISTART_SEED
    MPD( AddSolverOption_MergeDuplicates("alg:xslp_multistart_seed XSLP_MULTISTART_SEED MULTISTART_SEED",
      "Random seed used for the automatic generation of initial point when loading multistart presets "
      "\n\nDefault: 0",
      XSLP_MULTISTART_SEED, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_MULTISTART_SEED

#ifdef XSLP_MULTISTART_THREADS
    MPD( AddSolverOption_MergeDuplicates("tech:xslp_multistart_threads XSLP_MULTISTART_THREADS MULTISTART_THREADS",
      "The maximum number of threads to be used in multistart"
      "\n\nDefault: -1 (determined by XSLP_THREADS)",
      XSLP_MULTISTART_THREADS, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_MULTISTART_THREADS

#ifdef XSLP_MVTOL
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_mvtol XSLP_MVTOL SLPMVTOL",
      "Marginal value tolerance for determining if a constraint is slack"
      "\n\nDefault: -1.0",
      XSLP_MVTOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_MVTOL

#ifdef XSLP_NLPSOLVER
    MPD( AddSolverOption_MergeDuplicates("alg:xslp_nlpsolver XSLP_NLPSOLVER NLPSOLVER",
      "Controls whether to call FICO Xpress Global or one of the local solvers"
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  If the license allows and there are no user functions or multistart jobs, FICO Xpress Global will be called, otherwise a local solver."
      "\n- (1)  The algorithm selected by XSLP_SOLVER will be used to find a locally optimal solution"
      "\n- (2)  FICO Xpress Global will be used to find a globally optimal solution",
      XSLP_NLPSOLVER, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_NLPSOLVER

#ifdef XSLP_OBJTHRESHOLD
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_objthreshold XSLP_OBJTHRESHOLD SLPOBJTHRESHOLD",
      "Assumed maximum value of the objective function in absolute value."
      "\n\nDefault: 1.0e+15",
      XSLP_OBJTHRESHOLD, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_OBJTHRESHOLD

#ifdef XSLP_OBJTOPENALTYCOST
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_objtopenaltycost XSLP_OBJTOPENALTYCOST SLPOBJTOPENALTYCOST",
      "Factor to estimate initial penalty costs from objective function"
      "\n\nDefault: 0",
      XSLP_OBJTOPENALTYCOST, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_OBJTOPENALTYCOST

#ifdef XSLP_OCOUNT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_ocount XSLP_OCOUNT SLPOCOUNT",
      "Number of SLP iterations over which to measure objective function variation for static objective (2) convergence criterion"
      "\n\nDefault: 5",
      XSLP_OCOUNT, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_OCOUNT

#ifdef XSLP_OPTIMALITYTOLTARGET
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_optimalitytoltarget XSLP_OPTIMALITYTOLTARGET SLPOPTIMALITYTOLTARGET",
      "When set, this defines a target optimality tolerance to which the linearizations are solved to"
      "\n\nDefault: 0 (ignored, not set)",
      XSLP_OPTIMALITYTOLTARGET, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_OPTIMALITYTOLTARGET

#ifdef XSLP_OTOL_A
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_otol_a XSLP_OTOL_A SLPOTOL_A",
      "Absolute static objective (2) convergence tolerance"
      "\n\nDefault: -1.0",
      XSLP_OTOL_A, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_OTOL_A

#ifdef XSLP_OTOL_R
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_otol_r XSLP_OTOL_R SLPOTOL_R",
      "Relative static objective (2) convergence tolerance"
      "\n\nDefault: -1.0",
      XSLP_OTOL_R, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_OTOL_R

#ifdef XSLP_PENALTYCOLFORMAT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_penaltycolformat XSLP_PENALTYCOLFORMAT SLPPENALTYCOLFORMAT",
      "Formatting string for creation of the names of the SLP penalty transfer vectors"
      "\n\nDefault: pPC_%s where p is a unique prefix for names in the current problem",
      XSLP_PENALTYCOLFORMAT) );
#endif  // ifdef XSLP_PENALTYCOLFORMAT

#ifdef XSLP_PENALTYINFOSTART
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_penaltyinfostart XSLP_PENALTYINFOSTART SLPPENALTYINFOSTART",
      "Iteration from which to record row penalty information"
      "\n\nDefault: 3",
      XSLP_PENALTYINFOSTART, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_PENALTYINFOSTART

#ifdef XSLP_PENALTYROWFORMAT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_penaltyrowformat XSLP_PENALTYROWFORMAT SLPPENALTYROWFORMAT",
      "Formatting string for creation of the names of the SLP penalty rows"
      "\n\nDefault: pPR_%s where p is a unique prefix for names in the current problem",
      XSLP_PENALTYROWFORMAT) );
#endif  // ifdef XSLP_PENALTYROWFORMAT

#ifdef XSLP_PLUSDELTAFORMAT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_plusdeltaformat XSLP_PLUSDELTAFORMAT SLPPLUSDELTAFORMAT",
      "Formatting string for creation of names for SLP positive penalty delta vectors"
      "\n\nDefault: pD+%s where p is a unique prefix for names in the current problem",
      XSLP_PLUSDELTAFORMAT) );
#endif  // ifdef XSLP_PLUSDELTAFORMAT

#ifdef XSLP_PLUSERRORFORMAT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_pluserrorformat XSLP_PLUSERRORFORMAT SLPPLUSERRORFORMAT",
      "Formatting string for creation of names for SLP positive penalty error vectors"
      "\n\nDefault: pE+%s where p is a unique prefix for names in the current problem",
      XSLP_PLUSERRORFORMAT) );
#endif  // ifdef XSLP_PLUSERRORFORMAT

#ifdef XSLP_POSTSOLVE
    MPD( AddSolverOption_MergeDuplicates("pre:xslp_postsolve XSLP_POSTSOLVE NLPPOSTSOLVE",
      "This control determines whether postsolving should be performed automatically"
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Postsolve if the problem could be solved to optimality/infeasibility."
      "\n- (0)  Do not automatically postsolve."
      "\n- (1)  Postsolve automatically.",
      XSLP_POSTSOLVE, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_POSTSOLVE

#ifdef XSLP_PRESOLVE
    MPD( AddSolverOption_MergeDuplicates("pre:xslp_presolve XSLP_PRESOLVE NLPPRESOLVE",
      "This control determines whether presolving should be performed on the nonlinear problem prior to starting the main algorithm"
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Disable nonlinear presolve if and only if Optimizer presolve is disabled."
      "\n- (0)  Disable nonlinear presolve."
      "\n- (1)  Activate nonlinear presolve."
      "\n- (2)  Low memory presolve. Original problem is not restored by postsolve and dual solution may not be completely postsolved.",
      XSLP_PRESOLVE, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_PRESOLVE

#ifdef XSLP_PRESOLVELEVEL
    MPD( AddSolverOption_MergeDuplicates("pre:xslp_presolvelevel XSLP_PRESOLVELEVEL NLPPRESOLVELEVEL",
      "This control determines the level of changes presolve may carry out on the problem and whether column/row indices may change"
      "\n\n"
      "Values (default: XSLP_PRESOLVELEVEL_FULL):\n"
      "\n- (1)  Individual rows only presolve, no dropped columns/rows or index changes, no nonlinear transformations (XSLP_PRESOLVELEVEL_LOCALIZED)."
      "\n- (2)  All linear presolve that does not drop columns/rows, no index changes, no nonlinear transformations (XSLP_PRESOLVELEVEL_BASIC)."
      "\n- (3)  Full linear presolve including dropping columns/rows and index changes, no nonlinear transformations (XSLP_PRESOLVELEVEL_LINEAR)."
      "\n- (4)  Full presolve (XSLP_PRESOLVELEVEL_FULL).",
      XSLP_PRESOLVELEVEL, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_PRESOLVELEVEL

#ifdef XSLP_PRESOLVEOPS
    MPD( AddSolverOption_MergeDuplicates("pre:xslp_presolveops XSLP_PRESOLVEOPS NLPPRESOLVEOPS",
      "Bitmap indicating the SLP presolve actions to be taken"
      "\n\n"
      "Values (default: 2104):\n"
      "\n- (0)  Generic SLP presolve."
      "\n- (1)  Explicitly fix columns identified as fixed to zero."
      "\n- (2)  Explicitly fix all columns identified as fixed."
      "\n- (3)  SLP bound tightening."
      "\n- (4)  MISLP bound tightening."
      "\n- (5)  Bound tightening based on function domains."
      "\n- (8)  Do not presolve coefficients."
      "\n- (9)  Do not remove delta variables."
      "\n- (10)  Avoid reductions that can not be dual postsolved."
      "\n- (11)  Allow eliminations on determined variables."
      "\n- (12)  Avoid performing linear reductions at the nlp level."
      "\n- (13)  Avoid simplifying nonlinear expressions.",
      XSLP_PRESOLVEOPS, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_PRESOLVEOPS

#ifdef XSLP_PRESOLVEZERO
    MPD( AddSolverOption_MergeDuplicates("pre:xslp_presolvezero XSLP_PRESOLVEZERO NLPPRESOLVEZERO",
      "Minimum absolute value for a variable which is identified as nonzero during SLP presolve"
      "\n\nDefault: 1.0E-09",
      XSLP_PRESOLVEZERO, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_PRESOLVEZERO

#ifdef XSLP_PRESOLVE_ELIMTOL
    MPD( AddSolverOption_MergeDuplicates("pre:xslp_presolve_elimtol XSLP_PRESOLVE_ELIMTOL NLPPRESOLVE_ELIMTOL",
      "Tolerance for nonlinear eliminations during SLP presolve"
      "\n\nDefault: 0.001",
      XSLP_PRESOLVE_ELIMTOL, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_PRESOLVE_ELIMTOL

#ifdef XSLP_PRIMALINTEGRALALPHA
    MPD( AddSolverOption_MergeDuplicates("log:xslp_primalintegralalpha XSLP_PRIMALINTEGRALALPHA NLPPRIMALINTEGRALALPHA",
      "Decay term for primal integral computation"
      "\n\nDefault: 0",
      XSLP_PRIMALINTEGRALALPHA, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_PRIMALINTEGRALALPHA

#ifdef XSLP_PRIMALINTEGRALREF
    MPD( AddSolverOption_MergeDuplicates("log:xslp_primalintegralref XSLP_PRIMALINTEGRALREF NLPPRIMALINTEGRALREF",
      "Reference solution value to take into account when calculating the primal integral"
      "\n\nDefault: XPRS_PLUSINFINITY",
      XSLP_PRIMALINTEGRALREF, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_PRIMALINTEGRALREF

#ifdef XSLP_PROBING
    MPD( AddSolverOption_MergeDuplicates("pre:xslp_probing XSLP_PROBING NLPPROBING",
      "This control determines whether probing on a subset of variables should be performed prior to starting the main algorithm. Probing runs multiple times bound reduction in order to further tighten the bounding box."
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  Automatic."
      "\n- (0)  Disable SLP probing."
      "\n- (1)  Activate SLP probing only on binary variables."
      "\n- (2)  Activate SLP probing only on binary or unbounded integer variables."
      "\n- (3)  Activate SLP probing only on binary or integer variables."
      "\n- (4)  Activate SLP probing only on binary, integer variables, and unbounded continuous variables."
      "\n- (5)  Activate SLP probing on any variable.",
      XSLP_PROBING, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_PROBING

#ifdef XSLP_REFORMULATE
    MPD( AddSolverOption_MergeDuplicates("pre:xslp_reformulate XSLP_REFORMULATE NLPREFORMULATE",
      "Controls the problem reformulations carried out before augmentation. This allows SLP to take advantage of dedicated algorithms for special problem classes. "
      "\n\n"
      "Values (default: 511 (bits 0 — 8 incl. are set)):\n"
      "\n- (0)  Solve convex quadratic objectives using the XPRS library ."
      "\n- (1)  Convert non-convex quadratic objectives to SLP constructs ."
      "\n- (2)  Solve convex quadratic constraints using the XPRS library."
      "\n- (3)  Convert non-convex QCQP constraints to SLP constructs."
      "\n- (4)  Keep second order cones in the XPRS problem to keep them in the linearizations."
      "\n- (5)  Convexity of a quadratic only problem may be checked by calling the optimizer to solve the instance."
      "\n- (6)  Convert pievewise linear functions to MIP constructs."
      "\n- (7)  Convert ABS functions to MIP constraints if the full problem can be made not nonlinear."
      "\n- (8)  Convert MIN and MAX functions to MIP expressions if the full problem can be made not nonlinear."
      "\n- (9)  Always convert ABS expressions."
      "\n- (10)  Always convert MIN and MAX expressions.",
      XSLP_REFORMULATE, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_REFORMULATE

#ifdef XSLP_SAMECOUNT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_samecount XSLP_SAMECOUNT SLPSAMECOUNT",
      "Number of steps reaching the step bound in the same direction before step bounds are increased"
      "\n\nDefault: 3",
      XSLP_SAMECOUNT, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_SAMECOUNT

#ifdef XSLP_SAMEDAMP
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_samedamp XSLP_SAMEDAMP SLPSAMEDAMP",
      "Number of steps in same direction before damping factor is increased"
      "\n\nDefault: 3",
      XSLP_SAMEDAMP, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_SAMEDAMP

#ifdef XSLP_SBLOROWFORMAT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_sblorowformat XSLP_SBLOROWFORMAT SLPSBLOROWFORMAT",
      "Formatting string for creation of names for SLP lower step bound rows"
      "\n\nDefault: pSB-%s where p is a unique prefix for names in the current problem",
      XSLP_SBLOROWFORMAT) );
#endif  // ifdef XSLP_SBLOROWFORMAT

#ifdef XSLP_SBNAME
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_sbname XSLP_SBNAME SLPSBNAME",
      "Name of the set of initial step bounds to be used"
      "\n\nDefault: none",
      XSLP_SBNAME) );
#endif  // ifdef XSLP_SBNAME

#ifdef XSLP_SBROWOFFSET
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_sbrowoffset XSLP_SBROWOFFSET SLPSBROWOFFSET",
      "Position of first character of SLP variable name used to create name of SLP lower and upper step bound rows"
      "\n\nDefault: 0",
      XSLP_SBROWOFFSET, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_SBROWOFFSET

#ifdef XSLP_SBSTART
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_sbstart XSLP_SBSTART SLPSBSTART",
      "SLP iteration after which step bounds are first applied"
      "\n\nDefault: 8",
      XSLP_SBSTART, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_SBSTART

#ifdef XSLP_SBUPROWFORMAT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_sbuprowformat XSLP_SBUPROWFORMAT SLPSBUPROWFORMAT",
      "Formatting string for creation of names for SLP upper step bound rows"
      "\n\nDefault: pSB+%s where p is a unique prefix for names in the current problem",
      XSLP_SBUPROWFORMAT) );
#endif  // ifdef XSLP_SBUPROWFORMAT

#ifdef XSLP_SCALE
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_scale XSLP_SCALE SLPSCALE",
      "When to re-scale the SLP problem"
      "\n\n"
      "Values (default: 1):\n"
      "\n- (0)  No re-scaling."
      "\n- (1)  Re-scale every SLP iteration up to XSLP_SCALECOUNT iterations after the end of barrier optimization."
      "\n- (2)  Re-scale every SLP iteration up to XSLP_SCALECOUNT iterations in total."
      "\n- (3)  Re-scale every SLP iteration until primal simplex is automatically invoked."
      "\n- (4)  Re-scale every SLP iteration."
      "\n- (5)  Re-scale every  XSLP_SCALECOUNT SLP iterations."
      "\n- (6)  Re-scale every  XSLP_SCALECOUNT SLP iterations after the end of barrier optimization.",
      XSLP_SCALE, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_SCALE

#ifdef XSLP_SCALECOUNT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_scalecount XSLP_SCALECOUNT SLPSCALECOUNT",
      "Iteration limit used in determining when to re-scale the SLP matrix"
      "\n\nDefault: 0",
      XSLP_SCALECOUNT, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_SCALECOUNT

#ifdef XSLP_SHRINK
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_shrink XSLP_SHRINK SLPSHRINK",
      "Multiplier to reduce a step bound"
      "\n\nDefault: 0.5",
      XSLP_SHRINK, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_SHRINK

#ifdef XSLP_SHRINKBIAS
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_shrinkbias XSLP_SHRINKBIAS SLPSHRINKBIAS",
      "Defines an overwrite / adjustment of step bounds for improving iterations"
      "\n\nDefault: 0 (ignored, not set)",
      XSLP_SHRINKBIAS, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_SHRINKBIAS

#ifdef XSLP_SLPLOG
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_slplog XSLP_SLPLOG SLPLOG",
      "Frequency with which SLP status is printed"
      "\n\nDefault: 1",
      XSLP_SLPLOG, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_SLPLOG

#ifdef XSLP_SOLVER
    MPD( AddSolverOption_MergeDuplicates("alg:xslp_solver XSLP_SOLVER LOCALSOLVER",
      "Selects the library to use for local solves"
      "\n\n"
      "Values (default: -1):\n"
      "\n- (-1)  automatic selection, based on model characteristics and solver availability"
      "\n- (0)  use Xpress-SLP (always available)"
      "\n- (1)  use Knitro if available"
      "\n- (2)  use Xpress-Optimizer if possible (convex quadratic problems only)",
      XSLP_SOLVER, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_SOLVER

#ifdef XSLP_STOL_A
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_stol_a XSLP_STOL_A SLPSTOL_A",
      "Absolute slack convergence tolerance"
      "\n\nDefault: -1.0",
      XSLP_STOL_A, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_STOL_A

#ifdef XSLP_STOL_R
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_stol_r XSLP_STOL_R SLPSTOL_R",
      "Relative slack convergence tolerance"
      "\n\nDefault: -1.0",
      XSLP_STOL_R, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_STOL_R

#ifdef XSLP_STOPOUTOFRANGE
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_stopoutofrange XSLP_STOPOUTOFRANGE NLPSTOPOUTOFRANGE",
      "Stop optimization and return error code if internal function argument is out of range"
      "\n\nDefault: 0",
      XSLP_STOPOUTOFRANGE, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_STOPOUTOFRANGE

#ifdef XSLP_THREADS
    MPD( AddSolverOption_MergeDuplicates("tech:xslp_threads XSLP_THREADS NLPTHREADS",
      "Default number of threads to be used"
      "\n\nDefault: -1 (use XPRS_THREADS value)",
      XSLP_THREADS, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_THREADS

#ifdef XSLP_THREADSAFEUSERFUNC
    MPD( AddSolverOption_MergeDuplicates("func:xslp_threadsafeuserfunc XSLP_THREADSAFEUSERFUNC NLPTHREADSAFEUSERFUNC",
      "Defines if user functions are allowed to be called in parallel"
      "\n\n"
      "Values (default: 0 (no parallel user function calls)):\n"
      "\n- (0)  user function are not thread safe, and will not be called in parallel"
      "\n- (1)  user functions are thread safe, and may be called in parallel",
      XSLP_THREADSAFEUSERFUNC, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_THREADSAFEUSERFUNC

#ifdef XSLP_TOLNAME
    MPD( AddSolverOption_MergeDuplicates("tech:xslp_tolname XSLP_TOLNAME SLPTOLNAME",
      "Name of the set of tolerance sets to be used"
      "\n\nDefault: none",
      XSLP_TOLNAME) );
#endif  // ifdef XSLP_TOLNAME

#ifdef XSLP_TRACEMASK
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_tracemask XSLP_TRACEMASK SLPTRACEMASK",
      "Mask of variable or row names that are to be traced through the SLP iterates"
      "\n\nDefault: none (no tracing)",
      XSLP_TRACEMASK) );
#endif  // ifdef XSLP_TRACEMASK

#ifdef XSLP_TRACEMASKOPS
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_tracemaskops XSLP_TRACEMASKOPS SLPTRACEMASKOPS",
      "Controls the information printed for XSLP_TRACEMASK. The order in which the information is printed is determined by the order of bits in XSLP_TRACEMASKOPS. "
      "\n\n"
      "Values (default: -1):\n"
      "\n- (0)  The variable name is used as a mask, not as an exact fit."
      "\n- (1)  Use mask to trace rows."
      "\n- (2)  Use mask to trace columns."
      "\n- (3)  Use mask to trace cascaded SLP variables."
      "\n- (4)  Show row / column category."
      "\n- (5)  Trace slack values."
      "\n- (6)  Trace dual values."
      "\n- (7)  Trace row penalty multiplier."
      "\n- (8)  Trace variable values (as returned by the lineariation)."
      "\n- (9)  Trace reduced costs."
      "\n- (10)  Trace slp value (value used in linearization and cascaded)."
      "\n- (11)  Trace step bounds."
      "\n- (12)  Trace convergence status."
      "\n- (13)  Trace line search.",
      XSLP_TRACEMASKOPS, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_TRACEMASKOPS

#ifdef XSLP_UNFINISHEDLIMIT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_unfinishedlimit XSLP_UNFINISHEDLIMIT SLPUNFINISHEDLIMIT",
      "The number of consecutive SLP iterations that may have an unfinished status before the solve is terminated."
      "\n\nDefault: 3",
      XSLP_UNFINISHEDLIMIT, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_UNFINISHEDLIMIT

#ifdef XSLP_UPDATEFORMAT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_updateformat XSLP_UPDATEFORMAT SLPUPDATEFORMAT",
      "Formatting string for creation of names for SLP update rows"
      "\n\nDefault: pU_%s where p is a unique prefix for names in the current problem",
      XSLP_UPDATEFORMAT) );
#endif  // ifdef XSLP_UPDATEFORMAT

#ifdef XSLP_UPDATEOFFSET
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_updateoffset XSLP_UPDATEOFFSET SLPUPDATEOFFSET",
      "Position of first character of SLP variable name used to create name of SLP update row"
      "\n\nDefault: 0",
      XSLP_UPDATEOFFSET, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_UPDATEOFFSET

#ifdef XSLP_VALIDATIONFACTOR
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_validationfactor XSLP_VALIDATIONFACTOR NLPVALIDATIONFACTOR",
      "Minimum improvement in validation targets to continue iterating"
      "\n\nDefault: 0.001",
      XSLP_VALIDATIONFACTOR, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_VALIDATIONFACTOR

#ifdef XSLP_VALIDATIONTARGET_K
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_validationtarget_k XSLP_VALIDATIONTARGET_K NLPVALIDATIONTARGET_K",
      "Optimality target tolerance"
      "\n\nDefault: 1e-6",
      XSLP_VALIDATIONTARGET_K, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_VALIDATIONTARGET_K

#ifdef XSLP_VALIDATIONTARGET_R
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_validationtarget_r XSLP_VALIDATIONTARGET_R NLPVALIDATIONTARGET_R",
      "Feasiblity target tolerance"
      "\n\nDefault: 1e-6",
      XSLP_VALIDATIONTARGET_R, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_VALIDATIONTARGET_R

#ifdef XSLP_VALIDATIONTOL_A
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_validationtol_a XSLP_VALIDATIONTOL_A NLPVALIDATIONTOL_A",
      "Absolute tolerance for the XSLPvalidate procedure"
      "\n\nDefault: 0.00001",
      XSLP_VALIDATIONTOL_A, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_VALIDATIONTOL_A

#ifdef XSLP_VALIDATIONTOL_K
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_validationtol_k XSLP_VALIDATIONTOL_K NLPVALIDATIONTOL_K",
      "Relative tolerance for the XSLPvalidatekkt procedure"
      "\n\nDefault: 0.00001",
      XSLP_VALIDATIONTOL_K, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_VALIDATIONTOL_K

#ifdef XSLP_VALIDATIONTOL_R
    MPD( AddSolverOption_MergeDuplicates("tol:xslp_validationtol_r XSLP_VALIDATIONTOL_R NLPVALIDATIONTOL_R",
      "Relative tolerance for the XSLPvalidate procedure"
      "\n\nDefault: 0.00001",
      XSLP_VALIDATIONTOL_R, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_VALIDATIONTOL_R

#ifdef XSLP_VCOUNT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_vcount XSLP_VCOUNT SLPVCOUNT",
      "Number of SLP iterations over which to measure static objective (3) convergence"
      "\n\nDefault: 0",
      XSLP_VCOUNT, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_VCOUNT

#ifdef XSLP_VLIMIT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_vlimit XSLP_VLIMIT SLPVLIMIT",
      "Number of SLP iterations after which static objective (3) convergence testing starts"
      "\n\nDefault: 0",
      XSLP_VLIMIT, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_VLIMIT

#ifdef XSLP_VTOL_A
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_vtol_a XSLP_VTOL_A SLPVTOL_A",
      "Absolute static objective (3) convergence tolerance"
      "\n\nDefault: -1.0",
      XSLP_VTOL_A, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_VTOL_A

#ifdef XSLP_VTOL_R
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_vtol_r XSLP_VTOL_R SLPVTOL_R",
      "Relative static objective (3) convergence tolerance"
      "\n\nDefault: -1.0",
      XSLP_VTOL_R, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_VTOL_R

#ifdef XSLP_WCOUNT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_wcount XSLP_WCOUNT SLPWCOUNT",
      "Number of SLP iterations over which to measure the objective for the extended convergence continuation criterion"
      "\n\nDefault: 0",
      XSLP_WCOUNT, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_WCOUNT

#ifdef XSLP_WTOL_A
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_wtol_a XSLP_WTOL_A SLPWTOL_A",
      "Absolute extended convergence continuation tolerance"
      "\n\nDefault: -1.0",
      XSLP_WTOL_A, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_WTOL_A

#ifdef XSLP_WTOL_R
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_wtol_r XSLP_WTOL_R SLPWTOL_R",
      "Relative extended convergence continuation tolerance"
      "\n\nDefault: -1.0",
      XSLP_WTOL_R, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_WTOL_R

#ifdef XSLP_XCOUNT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_xcount XSLP_XCOUNT SLPXCOUNT",
      "Number of SLP iterations over which to measure static objective (1) convergence"
      "\n\nDefault: 5",
      XSLP_XCOUNT, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_XCOUNT

#ifdef XSLP_XLIMIT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_xlimit XSLP_XLIMIT SLPXLIMIT",
      "Number of SLP iterations up to which static objective (1) convergence testing is performed"
      "\n\nDefault: 100",
      XSLP_XLIMIT, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_XLIMIT

#ifdef XSLP_XTOL_A
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_xtol_a XSLP_XTOL_A SLPXTOL_A",
      "Absolute static objective function (1) tolerance"
      "\n\nDefault: -1.0",
      XSLP_XTOL_A, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_XTOL_A

#ifdef XSLP_XTOL_R
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_xtol_r XSLP_XTOL_R SLPXTOL_R",
      "Relative static objective function (1) tolerance"
      "\n\nDefault: -1.0",
      XSLP_XTOL_R, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_XTOL_R

#ifdef XSLP_ZERO
    MPD( AddSolverOption_MergeDuplicates("tol:xslp_zero XSLP_ZERO NLPZERO",
      "Absolute tolerance"
      "\n\nDefault: 1.0E-15",
      XSLP_ZERO, -DBL_MAX, DBL_MAX) );
#endif  // ifdef XSLP_ZERO

#ifdef XSLP_ZEROCRITERION
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_zerocriterion XSLP_ZEROCRITERION SLPZEROCRITERION",
      "Bitmap determining the behavior of the placeholder deletion procedure"
      "\n\n"
      "Values (default: 0):\n"
      "\n- (0)  (=1) Remove placeholders in nonbasic SLP variables"
      "\n- (1)  (=2) Remove placeholders in nonbasic delta variables"
      "\n- (2)  (=4) Remove placeholders in a basic SLP variable if its update row is                            nonbasic"
      "\n- (3)  (=8) Remove placeholders in a basic delta variable if its update row                            is nonbasic and the corresponding SLP variable is nonbasic"
      "\n- (4)  (=16) Remove placeholders in a basic delta variable if the determining                            row for the corresponding SLP variable is nonbasic"
      "\n- (5)  (=32) Print information about zero placeholders",
      XSLP_ZEROCRITERION, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_ZEROCRITERION

#ifdef XSLP_ZEROCRITERIONCOUNT
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_zerocriterioncount XSLP_ZEROCRITERIONCOUNT SLPZEROCRITERIONCOUNT",
      "Number of consecutive times a placeholder entry is zero before being considered for deletion"
      "\n\nDefault: 0",
      XSLP_ZEROCRITERIONCOUNT, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_ZEROCRITERIONCOUNT

#ifdef XSLP_ZEROCRITERIONSTART
    MPD( AddSolverOption_MergeDuplicates("slp:xslp_zerocriterionstart XSLP_ZEROCRITERIONSTART SLPZEROCRITERIONSTART",
      "SLP iteration at which criteria for deletion of placeholder entries are first activated."
      "\n\nDefault: 0",
      XSLP_ZEROCRITERIONSTART, INT_MIN, INT_MAX) );
#endif  // ifdef XSLP_ZEROCRITERIONSTART

  }  // AddNonlinearOptions()

};  // class CompiledNonlinearOptions

}  // namespace mp
