Summary of recent updates to the AMPL MP Library
================================================

## 20220403
- *DivConstraint and DivConverter_MIP*
    ModelFlattener now receives the Div expression
    and MIPFlatConverter handles it via quadratics.

- *Fix NL input variable order*

- *Reduce default strict comparison tolerance*
    Change *cvt:mip:eps* default value to 1e-4.

- *Build on MacOS 12.3, in particular on Apple M1*
    Fixed linking on MacOS 12.3 and FindCPLEX.cmake.
    For Apple M1, manually set -DCMAKE_OSX_ARCHITECTURES="x86_64" in CMake when building
    with CPLEX 22.1 because it contains only Intel libraries.

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
 for all
