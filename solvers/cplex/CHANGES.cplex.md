Summary of recent updates to CPLEX for AMPL
==============================================


## 20250814
- Changes in MP
  - Improved preprocessing of logical
    and combinatorial expressions
    (options cvt:pre:unnest, cvt:pre:sort).
  - Option cvt:pre:boundlogarg (default 0) to bound
    arguments of logarithm nonnegative. Previously
    always done, sometimes deteriorating performance
    of nonlinear solvers.
- Options alg:numericalemphasis, tech:memoryemphasis.


## 20250801
- Changes in MP
  - Tolerances set by options pre:feastol,
    pre:feastolrel both need to be violated
    to produce a warning on contradicting
    variable/constraint bounds. Previously
    the preprocessor failed on any violation,
    without letting the solver try.
  - Options cvt:compl, cvt:compl:eps control
    complementarity reformulations.


## 20250617
- Option lp:opttol (opttol).
- Changes in MP
  - Multi-objective emulator: added support for
    objective-specific options via objective suffixes
    beginning with *option_*
  - Option *cvt:unnest*: bits 2 and 4 switch on
    inlining of linear and quadratic subexpressions
    produced during reformulations (by default on).
  - Options *cvt:pre:ctx2ineq*, *cvt:pre:ctx2count*
    to control context propagation into conditional
    comparisons #267.


## 20250429
- Fix a bug in parsing of quadratic expressions,
  which could wrongly parse products of unequal
  linear expressions, such as (x-3)*(x-z-5).
- Options alg:dual and pre:dual are coordinated for
	compatibility with the legacy ASL driver.


## 20250424
- Changes in MP
  - Option cvt:qp2pass (default even faster parsing
    of quadratics)
- Compatibility updates
  - Option netopt accepts numeric argument.
    For keyword-only, use option alg:network.
    Note the new default (compared to cplexasl):
    netopt=0.
  - Options netdisplay, netfeasibility, netfind,
    netoptimality, and netpricing
  - Options primal, dual (both dummy)
- Fixed output logging incl. options mipdisplay,
  mipinterval, lpdisplay, bardisplay, modisplay,
  netdisplay


## 20250329
- Changes in MP:
  - Option cvt:multoutcard to limit the size of
    out-multiplied QP expressions. Can improve speed
    on large models.
  - Improved parsing of quadratic expressions.


## 20250308
- Changes in MP.


## 20241212
- Releasing with CPLEX version 22.1.2


## 20241206
- Added support for sensitivity analysis (keyword `sensitivity`)
- Added keywords `mip:nodefile`, `tech:workdir`, `tech:workfilelim`


## 20240828
- Options `primalopt`, `dualopt`, `baropt`, `siftopt`, `netopt`, `bendersopt` are now
  flags (e.g. they have to be set without a value)


## 20240823
- CPLEX MP driver is now the default. To use the previous ASL-based driver set:
  `option solver cplexasl;`
- Added many options available in the ASL driver
- Notable changes with ASL driver:
   - Keyword 'basis_cond' is now 'kappa' and follows the standard MP implementation
   - Multiobjective optimization follows the MP standard implementation
   - Choosing solution algorithm for the LP problems or the initial MIP subproblem 
     (formerly `mipstartalg`) now follows the same standard as other MP drivers: 
     see `alg:method` in the `-=` output
   - Former option `timing` is now implemented through the corresponding
     MP feature. To print the number of available cores use
     `tech:numcores`.
   - Use MP feature `tech:writemodelonly` instead of `writeprob` + `nosolve`.
     Note that `nosolve` still applies when specifying `writemipstart`
- Conversion of basis status for constraints depending on sense

## 20240801
- Fix a problem that occured when reporting IIS
- Removed some noise from console output


## 20240724
- Option *acc:_all*
	- Useful to disable all reformulations (acc:_all=2),
		or force linearization (acc:_all=0).
- Option *cvt:prod*     
  - Controls reformulation of binary products into logical 
    constraints.
- Faster input of quadratic expressions.


## 20240617
- *Multi-objective emulator*
	- obj:multi=2 forces emulation, even if MO natively supported.
	- Fixed a bug in the objective degradation suffixes
		.objasbtol, .objreltol.


## 20240604
- Presolve division by constant, resulting in fewer constraints
- Fix no-solution case in multi-objective emulator


## 20240529
- *Multi-objective emulator*
	- All flat MP solvers support multi-objective mode (obj:multi=1),
		either natively, or via emulation.
	- Suffixes .objpriority, .objweight, .objabstol, .objreltol.
	- [BREAKING] Default intuitive handling of .objweight,
		see option obj:multi:weight, even when natively supported.


## 20240429
- [BREAKING] Merged `report_times` and `timing`; they 
  are now aliases, set the value to 1 to have basic info,
  to 2 to have more detailed info.

  
## 20240320
- *SOS constraints*.
  - Fixed handling of SOS2 constraints created by AMPL
    as reformulations of PL expressions (`option
    pl_linearize 1`, default; set to 0 to use the solver's
    native PL functions or MP linearization.)
  - Disallow repeated weights for SOS constraints
    (suffixes `.sosno`/`.ref`.)
- *Option `report_times`* 
- *Unused `acc:` options*.
  - The constraint acceptance options `acc:...`
    for non-handled constraints are ignored
    (previously triggered error.)


## 20240115
- *Solve result codes*.
  - List codes by running (solver) -!
  - [BREAKING] Standardized codes. Major changes:
    - 100-199 (solved?) means solution candidate
      provided, but can be suboptimal/infeasible
    - 300-349 means unbounded problem but
      feasible solution returned
    - 400-449 means limit/interrupt but feasible
  - [BREAKING] sol:chk:fail returns code 150 (solved?)
- Improved translation of *SOCP constraints*.
  - Options cvt:socp, cvt:socp2qc.
- Compact solution check warnings
- Fixed presolve of the power function #226.
- Added fixed model for `mip:basis`
- Added options:
  - `alg:method` and flags to choose the solution algorithm
  - `mip:nodemethod` to choose the algorithm used for MIP node problems
  - `lp:solutiontype` to choose if to use crossover to get a basis basis solution after barrier


## 20231122
- Public beta; implemented most [features](https://mp.ampl.com/features-guide.html) in MP:
  model export, warm start, input and output basis, feasibility relaxation,
  multiple solutions, kappa, unbounded rays, IIS, return MIP gap, 
  return dual bound
- Native Model support:
    - continuous and integer variables
    - multiple objectives, quadratic objectives
    - linear and quadratic constraints, indicator constraints and PL functions
    - special ordered sets (type 1 and 2)


## 20231117
- MP update: fixed graceful exit on Ctrl-C from AMPL in Linux
  and fixed issue with reading text-format NL files


## 20230815
- Fixed a bug causing repeated names for
  auxiliary variables and constraints.
- Option values can be assigned without '='.
- Changed default tolerance for strict comparisons
  to 0 (option cmp:eps, #102.)
- Fixed a bug where equivalent conditional
  comparisons were not unified.

 
