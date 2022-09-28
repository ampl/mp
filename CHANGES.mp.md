Summary of recent updates to the AMPL MP Library
================================================

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
