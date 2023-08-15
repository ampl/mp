Summary of recent updates to the AMPL MP Library
================================================

## 20230815
- Fixed a bug causing repeated names for
  auxiliary variables and constraints.
- Option values can be assigned without '='.
- Changed default tolerance for strict comparisons
  to 0 (option cmp:eps, #102.)
- Fixed a bug where equivalent conditional
  comparisons were not unified.


## 20230728
- Option 'tech:writesolution' #218.
- Option 'writeprob' ('tech:writemodel') ASL-compatible.
- Hint when 'writeprob' fails: use 'writesol'.


## 20230726
- Fixed inequalities of integer expressions with
  non-integer constants, see test_int_non_int.mod.
- Backend std feature WRITE_SOLUTION.
- Fixed parsing quoted string options.


## 20230724
- Option [solver_]auxfiles rc; transfers names
	of variables and constraints into the model;
	(solver)_options 'cvt:names=0-3' controls names.


## 20230714
- Print warnings in non-verbose mode too.
- 'barrier' equivalent to 'barrier=1' for integer options.


## 20230621
- Fix quadratic objective with repeated subexpressions.


## 20230616
- Smaller reformulations for conditional comparisons.
- Option *cvt:names* sets whether to read AMPL
  variable names or to provide generic names.


## 20230531
- Cones: recognize (affine_expr) >= y * exp(z/y)
  as exponential cone.
- Cones: recognize xy >= 1 as rotated SOC.
- Wrong solver options are gracefully reported via
  solve_message.


## 20230515
- *Recognize exponential conic constraints*.
  Exponential cones are recognized and passed to the
  solver, if supported.
  

## 20230424
- *Pass variable names* if read from a `col` file with the 
  same name of the `nl` file being read.
- *Fixed #203*: starting solution is now not passed to the 
  solver if empty.
  

## 20230321
- *Recognize second-order cones*
  Recognize SOCP constraints from algebra and pass them
  natively, or transform to quadratics.


## 20230207
- *Handle boolean constants* in ProblemFlattener.


## 20221228
- *More to #163*: ignore SOS with zero weights.


## 20221222
- *Fixed #195*: case-insensitive option synonyms.

- *Fixed #194*
   Report correct objno for feasibility problems in .sol file,
   so that AMPL can print "Objective = find feasible solution".


## 20221211
- *==> else*
   Implemented implication with 'else': *constr1* ==> *constr2* [else *constr3*]   

- *PLApproxRelTol, PLApproxDomain*
   Parameters to control piecewise-linear approximation.
   cvt:plapprox:reltol default value changed from 1e-5 to 0.01.

## 20221012
- *Piecewise-linear approximation of quadratics*
    Automatic for linear solvers.
    For convex QP solvers, set the following options:
    cvt:quadobj=0 cvt:quadcon=0 to linearize nonconvex objective(s)
    and/or constraints.
    Recognizing x^2 for stronger univariate approximation.

## 20220928
- *Piecewise-linear approximation of univariate nonlinear functions*
    Approximation of exp, a^x, x^a, log, log10, trigonometric and hyperbolic
    functions.

- *Default value of big-M*
    For linearization of logical constraints on variables without finite bounds,
    option cvt:mip:bigM can provide a default big-M bound.
    

## 20220725
- *Solution file export* 
    On Windows now creates files with LF only to avoid issues when exporting
    suffixes to AMPL.

    Multiple solutions export file format amended.


## 20220720
- *Propagating suffixes via expression trees into flat constraints*
    Partially implemented #184. x-gurobi accepts options
    'funcpieces...' and corresponding suffixes which are passed into
    GRBaddgenconstrExp etc.

    Subexpressions: note that if a subexpression is contained in several
    constraints, for contradicting suffix values the maximum is taken.

- *Option 'cvt:writegraph'*
   Exporting the flattening / conversion graph in JSON Lines format (WIP).


## 20220617
- *PowConstraint is reduced to quadratics in some cases*
    For constant non-negative integer exponent and base variables
    with negative lower bound, PowConstraint
    is reduced to quadratics (possibly with auxiliary variables).
    Reason: Gurobi's GRBaddgenconstrPow
    not accepting negative bases.

- *Context for algebraic constraints*
    Context is now propagated for algebraic constraints.
    For example, 3x + max(y, z) <= 6 will result in 3 linear
    constraints. (Earlier this was done for logical constraints
    and objectives).


## 20220526
- *Special ordered sets*
    Fixed: SOS are now recognized even if the suffix '.ref' 
    value is integer


## 20220511
- *Complementarity constraints: also quadratics*
    Complementarity constraints now handle quadratics.

- *Branch develop is used for new code*
    The active development branch is now *develop*.

- *Convert quadratic range constraints to QuadCon(LE/EQ/GE)*
    Gurobi and COPT do not support quadratic range constraints.
    (Gurobi's linear ranges are not feature-complete).
    Conversion of linear range constraints into one-side rhs
    constraints has been generalized for any algebraic ones.


## 20220408
- *Complementarity constraints: 1st go*
    Conversion to MIP of complementarity constraints
    (no quadratics but functional subexpressions ok).

- *DivConstraint and DivConverter_MIP*
    ModelFlattener now receives the Div expression
    and MIPFlatConverter handles it via quadratics.

- *Fix NL input variable order*

- *Reduce default strict comparison tolerance*
    Change *cvt:mip:eps* default value to 1e-4.

- *Build on MacOS 12.3, in particular on Apple M1*
    Fixed linking on MacOS 12.3 and FindCPLEX.cmake.
    For Apple M1, manually set -DCMAKE_OSX_ARCHITECTURES="x86_64"
    in CMake when building with CPLEX 22.1 because
    it contains only Intel libraries.

- *Expression maps*
    FlatConverter eliminates subexpressions of all types.
    A subexpression means here a duplicate expression, such as
    abs(x+2) occurring several times in the model (here x+2
    is a nested subexpression).

- *AMPLS C API*
    C API allowing access to underlying solver API.
    Replaces the previous Solver C API (solver-c.cc).
    Toy driver `gurobi_ampls` exemplifies API usage.


## 20220216
- *Improved warnings (#161, #163)*:
    In verbose mode, FlatConverter / Backend print warning summary
    before and after solving
    
- *Fixed solver option parsing in Windows (#160)*

- *Reworked Backend / ModelAPI class hierarchy (#162)*:
    In particular, also generalized old MP hierarchy
    (Solver / ProblemBuilder / AppSolutionHandler / SolverNLHandler)

- *Allowing SOS constraints with repeated weights (#163)*:
    Although Gurobi states SOS weights should be unique, it accepts them repeated.
    This happens when AMPL linearizes a PL function with redundant (repeated) slopes.
    It seems better to use PL functions natively (*option pl_linearize 0;*).
