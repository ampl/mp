Summary of recent updates to the AMPL MP Library
================================================

## Unreleased

- New option *cvt:pre:inteps* (synonyms *pre:inteps*, *roundintbnd*):
  integrality tolerance, default 1e-6. A value within it of an integer
  is taken as that integer for bounds of integer variables (otherwise
  they are rounded inwards), right-hand sides of comparisons of
  integer-valued expressions, compared constants and *numberof*
  values, also in logical expressions. Fix: such comparisons could give
  a wrong optimum, depending on how they were written: `x <= 0.3/0.1`
  (2.9999999999999996) allowed integer `x = 3` as a constraint, but not
  inside a disjunction; `x == 0.3/0.1` was false. Negative: exact.
  Values from 0.5 are rejected.
- Comparisons inside logical expressions are decided by the converter
  from the bounds of their expressions with option *cvt:pre:feastol*
  (absolute, in the units of the expression) and a margin for rounding
  errors of the compared constant (new option *cvt:pre:ulps*, default
  16 units in the last place; 0: none). A non-strict comparison holds
  within the tolerance: `x >= 4` with `x <= 3.9999995` is no longer
  decided false, and `(x + 1)/3 >= 2 ==> y >= 3` at `x <= 5` is no
  longer decided false by a rounding error (which made the returned
  solution violate the model). A strict comparison is decided true only
  if it holds exactly, and false as the negation of the non-strict one:
  `x < 1e-7` for integer `x >= 0` is false, with a warning.
  *cvt:pre:feastol=0* gives exact decisions; *cvt:pre:feastol=0
  cvt:pre:inteps=-1 cvt:pre:ulps=0* the previous results.
  *cvt:pre:all=0* also disables all rounding and tolerances.
- Fix: a non-strict comparison of a continuous expression close to its
  bound, for example `x <= 5e-5 or y >= 1` with `x` in `[0, 1e-4]`,
  was replaced by equality to the bound (`x == 0`), which cut off
  solutions or made the model infeasible; and `not (1e-5*x > 5e-5)`
  could make a model infeasible.
- An equality whose right-hand side differs from a bound of its
  expression only by rounding error, such as `0.3*b == 0.1*3` for a
  binary `b`, holds at that bound.
- Fix: `numberof x in (...)` with a variable `x` could give a wrong
  optimum when `x` came after the compared expressions in the NL file
  (for example, with `abs()` elsewhere in the model): preprocessing
  (*cvt:pre:sort*) reordered `x` with the compared values.
- Fix: `numberof` with compared expressions such as `max()`, `min()`,
  `abs()` or another `numberof` failed with "Part of expression
  redefinition could be lost" or, with a nested `numberof`, could give a
  wrong optimum.
- Starting points for continuous models (no unfixed integer variables,
  also after reformulation), and for solvers without MIP starts: incoming
  primal values (for example, from `let x := ...` or `var x := ...`)
  are now passed to the solver without dual values, and dual values
  without primal values. Previously a starting point was passed only
  when AMPL sent both, so some drivers ignored a primal start alone,
  for example for nonlinear barrier algorithms. Primal or dual values are passed only
  if AMPL sent them; a variable without an initial value gets 0, moved
  into its bounds. A warning reports incoming primal values that the
  driver does not pass for the model (for example, COPT takes an LP
  start only with both primal and dual values), and incoming dual values
  when no primal values are passed either. New option
  *alg:startpartial* (default 3); *alg:startpartial=0* restores the
  previous behavior. A solver or algorithm may still ignore an
  incomplete start. With *alg:start=2*, the incoming basis is not used
  also with an incomplete start.
- Clarified: *alg:start=0* does not disable an incoming LP basis; use
  *alg:basis=0* for that.
- Fix: the solver's name and version, for example "COPT 9.0.0: ", appeared
  twice at the start of the solve message wherever the output is not a
  terminal (notebooks, IDEs, logs, captured output). When AMPL runs the
  solver, drivers now ask AMPL to skip the repeated name in the solve
  message, as the ASL-based drivers do; the output in a terminal and AMPL's
  `solve_message` are unchanged. Run by another program with `-AMPL` (for
  example, Pyomo), drivers no longer print backspaces, and the message in
  the solution file has none.
- Fix: some errors, for example reading a corrupt NL file, were reported
  with solve result 1 (`solved`), so that AMPL kept the previous values of
  the variables as if they were a solution. Such errors are now reported as
  `failure` (500). The AMPLS C API functions now return an error's solve
  result code whenever it is above 99 (for example, 150 for a failed
  solution check with *chk:fail*), not only codes 200..599.
- AMPLS C API: new `AMPLSGetSolveResult()` gives the solve result code
  (AMPL's `solve_result_num`) of the last solve, or of the last failed
  call. The return value of `AMPLSSolve()` and the other calls is 0 on
  success, whatever the solve result.
- Fix: the solution check missed fractional values of integer variables
  unless they were below 0.5 in absolute value (for example, 10.5 passed).
  Now such values are reported ("variable integrality", or "aux var
  integrality" for variables added by reformulations). The default of
  *sol:chk:inttol* is raised from 1e-5 to 1e-4; the new *sol:chk:inttolrel*
  (default 0, off) adds a relative tolerance for very large values. If the
  solver runs with a larger integrality tolerance, raise *sol:chk:inttol*;
  *sol:chk:inttol=0.5* turns the integrality check off.
- Messages that refer to variables or constraints (solution check
  warnings, contradictory bounds, errors in constraints, SDP matrix
  errors) suggest `option <solver>_auxfiles rc;` when AMPL did not pass
  names; when it did, the AMPL name is shown next to `_svar[i]` or
  `_scon[i]`. Warnings about contradictory bounds of a variable are now
  titled "Contradictory variable bounds" and, for bounds deduced by the
  conversion, "Contradictory deduced bounds" (were "_svar bounds" and
  "_svar best-known bounds").
- Fix: when a model was solved repeatedly through the AMPLS API, the solve
  message of each report repeated the notes of the earlier solves.
- Option `cvt:pre:vbnd` is now called `cvt:pre:domsync` (synonym
  `domsync`); `cvt:pre:vbnd` and `vbnd` remain accepted.
- `1^x` in an NL file not written by AMPL (AMPL simplifies it) is taken as
  the constant 1, instead of an exponential with base 1, which the
  piecewise-linear approximation could not handle.
- Fixed false infeasibility when linearizing a product with a bounded integer
  factor whose original NL bounds are fractional.
- Preserve exact integer bounds during rounding, including large bounds that
  enforce one side of a complementarity constraint.
- Fixed `div`, `mod`, `floor`, `ceil`, `round` and `trunc`:
  - the bounds of their results that follow from the bounds of the
    arguments are now passed to the solver also when *cvt:pre:boundsbest*
    is off (MP2NL's default). Without them, a solver behind MP2NL that
    linearizes logical expressions stopped with "Set bounds on variables
    participating in logical expressions", and SCIP could report a wrong
    optimum;
  - `x div y` with a variable divisor whose `abs(y)` also appears in the
    model stopped with "Part of expression redefinition could be lost".
- Contradictory equalities of an integer variable (`x == 1 and x == -1`)
  stopped with the error "vector"; the model is now reported infeasible.
- Fixed *alg:feasrelax* with Gurobi, CPLEX, COPT and HiGHS: constraints and
  variables that the conversion adds (for example, for `abs()`, `max()` or
  logical expressions) could be relaxed at no cost, so a reported "solved"
  relaxation could violate the original constraints. They are no longer
  relaxed. Limitations: logical constraints are in general not relaxed
  (their linearization consists of added constraints, and solvers do not
  relax native ones); and the converter deduces variable bounds from the
  constraints, which are not relaxed: they can simplify other constraints
  away, and contradictory deduced bounds are passed to the solver
  (*cvt:pre:infeasvarbnd=0*), so a relaxation can be costlier than needed
  or still violate a constraint. AMPL's presolve also turns constraints of
  one variable into bounds; use `option presolve 0;` to relax them as
  constraints.
- New warning when the domain of a divisor of `/`, `div` or `mod`
  contains 0: the divisor is treated as nonzero in the whole model, also
  where the division is not evaluated (for example, in an if-then-else
  branch not taken).
- Warnings and error messages are wrapped at 60 characters (tables and
  explicit line breaks are kept). New option *tech:msgwidth* sets the
  width; 0 restores the previous, unwrapped output.

- Fixed wrong results of the conversion of logical constructs (all solvers):
  - a comparison shared by a disjunction and an if-then-else condition lost
    one of its two implications when the if-then-else result was bounded by a
    strict comparison (`<`, `>`);
  - an integer variable fixed at a fractional value (`2*x == 1`) was made
    continuous and the model reported optimal;
  - `abs()` of a variable fixed at a negative value made `abs(x) > c` infeasible;
  - comparisons `x == c` for constants outside the domain of a variable that was
    fixed later were left free in the unary encoding, so infeasible models
    were reported optimal;
  - contradictory bounds of an auxiliary variable stopped the conversion with
    an error instead of reporting infeasibility;
  - a domain hole of a variable (from `not (x != 0 and x != 3)`) made the
    rewriting of `x == ub` and `x >= ub` recurse endlessly (stack overflow).

- Option `cvt:pre:infeasvarbnd`: with 1, the converter declares the model
  infeasible when the bounds of a variable, also best-known bounds deduced
  by the converter, become contradictory. Default 0 (as before): a warning
  is given and the bounds are passed on to the solver; Xpress, SCIP and
  COPT then stop with an error if such a variable is an indicator binary
  ("not a binary variable"). MP2NL uses default 1, because several NL
  subsolvers stop with an error on such bounds.
- The bounds of integer and implied-integer variables are rounded inwards;
  a bound within the integrality tolerance *cvt:pre:inteps* (see above)
  of an integer is that integer.
- An indicator constraint whose flag became fixed only after the indicator
  was created (for example, `if x <= 1 then ...` with `x <= 1` implied by
  another constraint) made Xpress and SCIP stop with an error ("not a binary
  variable"). Such indicators are now replaced by their implied constraint,
  or dropped: option `cvt:pre:indfixedflag` (default 1; 0 restores the old
  behavior; off with `cvt:pre:all=0`).
- Linear (and quadratic) subexpressions are now also inlined into the
  expressions defining other variables, for example the linearization of
  `not` inside an equivalence. With NL output (MP2NL, `acc:_expr=1`), such
  linear terms were written as nonlinear expressions; CPLEX and Xpress
  through MP2NL then rejected the model ("not convex quadratic"), and
  nonlinear solvers received linear parts as nonlinear. New bit 32 of
  option `cvt:pre:unnest` (default now 63).
- The inlining of shared subexpressions is limited by the fill it can cause
  (upper bound (u-1)(k-1)-2 for k terms and u uses): options
  `cvt:pre:agglinfill` (default 10) and `cvt:pre:aggquadfill` (quadratic terms
  only, default -1: no limit). A linear subexpression not inlined stays a
  variable with its linear definition and the bounds implied by it.
  `cvt:pre:unnest=31 cvt:pre:agglinfill=-1` restores the old behavior.
- Fixed a wrong result with expressions passed to the solver (e.g., Gurobi 14,
  `acc:_expr=1`): a function result relaxed to one inequality by context
  (e.g., `v <= exp(.)` inside `log(1 + v^2)`) lost its implied bounds
  (`v >= 0`), on which the context relied. They are now passed to the solver.

- Fixed wrong results of the multi-objective emulator with nonlinear
  objectives (for example `exp(x)`) when the first objective has a negative
  `.objweight`: the pass objective's nonlinear part was relaxed in the wrong
  direction, so a solver could report an unbounded or a wrong optimal
  solution. Solvers that take nonlinear objectives natively were affected
  with `cvt:nlobj=0` or when the objective was reformulated.
- Fixed a wrong objective in models with a convex quadratic objective whose
  squared terms have a coefficient other than 1 (for example `maximize
  3*x - 0.8*x^2`) when the model also has quadratic constraints recognized
  as second-order cones: the objective's cone used the coefficient instead of
  its square root, so the objective term was off (0.8 became 0.64).
  Observed with BARON through MP.
- Fixed second-order cones whose first variable is fixed by its bounds,
  for solvers that receive the cone as a quadratic constraint (for example
  Gurobi, Xpress, CPLEX, BARON): `sqrt(x^2 + y^2) <= t` with `t` fixed at 3
  became `x^2 + y^2 <= 3` instead of `<= 9`. A negative fixed value keeps
  the constraint infeasible.
- Fixed a regression in 20261029: with variable names passed to the solver
  (for example `option highs_auxfiles rc;`), HiGHS failed some logical
  models with "Error -1 for call ... Highs_passColName". Columns that MP
  adds for constant values now get generated names (`_svar<n>_`).
- Result bounds of `exp`, `a^x`, `log` and `log10`/`log_a` now follow
  from the bounds of their argument (they were only `>= 0` for the
  exponentials and unbounded for the logarithms). Logical and indicator
  constraints that contain such a function with a bounded argument, for
  example `b == 0 ==> exp(x) <= 2` with `x <= 5`, no longer stop with
  "Set bounds on variables participating in logical expressions" when the
  solver needs big-M (`acc:ind...=0`, or no native indicators). New option
  *cvt:pre:funcbndtol* (default 1e-12) is the relative margin added to the
  computed bounds against rounding; a negative value restores the previous
  behavior.
- Fixed a regression in 20261029: a second `AMPLSSolve()` on the same model
  (AMPLS C API) returned without solving. It now solves again; with
  emulated multiple objectives it reports that a repeated solve is not
  supported.
- AMPLS C API: `AMPLSSolve()` now returns 0 on success; on failure, a
  nonzero solve result code or -1, as `AMPLSLoadNLModel()` does. An
  error during the solve (for example, a failing solver call or the
  repeated solve above) no longer escapes from `AMPLSSolve()` as a C++
  exception, which could terminate the calling program. Its message is
  among `AMPLSGetMessages()`. After a failed solve, report the error with
  `AMPLSReportError()` if needed, not with `AMPLSReportResults()`, which
  would report the previous solve. Existing callers keep working.

## 20261029

**Compatibility and correctness**

- **Original piecewise-linear expressions** (AMPL `pl_linearize=0`, recommended):
  fixed incorrect solutions for chained terms and breakpoint-related failures
  without native SOS2. Nonconvex MIP reformulations now require finite argument
  bounds instead of silently truncating at +-1e6; `cvt:pl:domain` explicitly
  permits truncation, with a warning (use cautiously). Native SOS2 now handles
  unbounded argument domains exactly. AMPL's default `pl_linearize=1` bypasses
  this MP reformulation path. See [MP Modeling Guide: Numerical accuracy — Piecewise-linear functions](https://mp.ampl.com/modeling-numerics.html#piecewise-linear-functions).
- Fixed underestimated objectives from indicator linearization and unrelated
  constraints marked as IIS members in models with logical/conditional expressions.
  Fixed `cvt:quadcon=0` for quadratic range constraints, including solvers
  with native quadratic support.
- Fixed **multiobjective emulation** (`obj:multi=2`, or `obj:multi=1` without
  native support) producing wrong solutions when opposite-sense objectives
  share a priority.
- **MP2NL:** `cvt:quadcon=0` now requests quadratic linearization in NL output
  (previously ignored). See [MP Modeling Guide: Tools & details — Meta-driver MP2NL](https://mp.ampl.com/modeling-tools.html#mp2nl)
  for native output, complementarity and subsolver settings.

**New modeling support**

- Added `floor`, `ceil`, `round(x[,n])` (ties to even), `trunc(x[,n])`, `div`
  and `mod`, including real-valued and variable divisors. MIP floor/ceil
  reformulations use `cvt:mip:eps`.
- Added exact MIP reformulations for products with a bounded integer factor
  and integer powers.
  `cvt:mip:intprodbits` defaults to 30; 0 restores piecewise-linear approximation.
  CP-SAT supports exact integer-variable products natively (`acc:mul`), with
  `cvt:multoutcard=1` retaining products of sums. See
  [MP Modeling Guide: Expressions supported — Quadratic and power operators](https://mp.ampl.com/modeling-expressions.html#quadratic-and-power-operators).

**Options and reporting**

- Added preprocessing controls for disjunctions (`cvt:pre:disj2ind`, default 11)
  and variable-bound logic (`cvt:pre:vbnd`, default 3). Fixed Xpress rejection
  of indicator binaries reused from conditional expressions.
- Extended `cvt:pre:boundsbest`: 2/3 select inferred bounds by expression type;
  4 uses inferred bounds for linear models and hard bounds otherwise. MP2NL
  uses 4 in its nonlinear-subsolver configurations. `cvt:pre:auxint` defaults
  to 0, allowing continuous auxiliaries where integrality is implied; 1 keeps
  integer auxiliaries (CP-SAT defaults to 1). See [MP Modeling Guide: Efficiency considerations — Automatic bounds and auxiliary integrality](https://mp.ampl.com/modeling-efficiency.html#automatic-bounds-and-auxiliary-integrality)
  and driver option help.
- Added MIP basis recovery for linear algebraic MIPs in COPT, MOSEK, HiGHS and
  SCIP, plus independent `mip:basis:timelimit` and `mip:basis:iterlimit` controls.
- Improved emulated multiobjective pass statistics: separate bounds from gaps
  and report stopped/unrun passes. Fixed Debug assertions when requesting
  ordinary gap/bound suffixes with several objectives; these suffixes do not
  provide independent certificates for every objective. Use `obj:multi:stats=1`
  and the [MP Modeling Guide: Multiple objectives — Multi-objective solve statistics (experimental)](https://mp.ampl.com/modeling-mo.html#multi-objective-solve-statistics-experimental).
- Solution checking now applies `sol:chk:feastol` to recomputed conditional
  comparisons and recomputes only eliminated auxiliaries.
- Fixed objective-conversion failures with `cvt:quadobj=0`, MOSEK dual-ray
  requests for constraint-free models, missing auxiliary names in
  `tech:writegraph`, and rejection of 64-bit option values above 2^31.
- **AMPLS:** fixed incomplete emulated multiobjective solves and crashes through
  `AMPLSSolve`, including ASL drivers with `mp2nl=1`. `AMPLSReportError` writes
  SOL error results after model-loading or result-reporting failures.


## 20260908
- Set membership operator discouraged
  during rework, see Modeling Guide.
- Don't fail when *mip:return_gap* or *mip:bestbound*
  ask for values but they are not available
  (warning issued).
- Option *obj:multi:stats=1* outputs multi-objective
  pass statistics in the *.objpass[_...]* suffixes.
- Fix a bug in complementarity reformulation
  with *cvt:compl=2,3*.
- Options *mip:plateau:...* control MIP
  tailing-off in supporting solvers.
  Can be used in multi-objective solves,
  see option *obj:multi:options*. 
- Option *cvt:nlobj* controls whether nonlinear
  objective terms are passed to the objective,
  vs moved to a constraint or piecewise-linearly
  approximated. Setting *cvt:nlobj=0*
  *cvt:quadobj=0* enables native
  multiobjective method (*obj:multi=1*)
  if only supported for linear objectives.
- Option *cvt:pre:unnest*: bit 4 (value 16) turns on
  inlining of algebraic subexpressions in indicator
  constraints. When it is off, the body of indicator
  constraints is always just a variable
  (possibly assigned to a larger expression).
  On by default.


## 20260624
- CMake option *MP_LINK_WITH_SHARED_LIB* links the
  driver as a shared library (Gurobi, IBM ILOG CPLEX,
  Cardinal Optimizer, FICO Xpress;
  for them, it's the default).
- Option *debug=1* outputs all options in problem
  suffixes *option___...*.
- Option *cvt:pow2_as_qp* prefers quadratization
  of the ^2 operator for affine arguments (default).
- Faster parsing of quadratic expressions with
  defined variables.


## 20260515
- Fixed a bug in solve result name in the Multi-Objective
  Emulator.
- Similar to inlining, nonlinear expressions with nontrivial
  external bounds are explicified into an auxiliary variable,
  depending on options *cvt:pre:boundsbest* and
  *cvt:pre:boundlogarg*.
- Fixed a debugging feature on Windows.
- Fixed a non-critical bug in the inlining of algebraic
  subexpressions resulting from redefinitions,
  option *cvt:pre:unnest*. In particular
  with *cvt:pre:boundsbest=1*, subexpressions
  with stronger result bounds are not inlined.  


## 20260414
- Fixed a bug in the multi-objective emulator
  (see option *obj:multi*) when nonlinear terms
  were not reset for consecutive objectives.
- Option *cvt:pre:continuous_fixed_vars* to declare
  fixed variables continuous. Default 1.
- Option *cvt:pre:boundsbest* defaults to making
  auxiliary subexpression result variables free,
  unless the result is constant.  This allows
  stronger nonlinear presolve.
- Additional case where the signpow function
  is recognized: abs(C*x)*x
  (C constant, x variable).
- Auxiliary variables and constraints receive names
  describing the redefinition path. Such names are
  only given if AMPL option (solver_)auxfiles rc is
  set, or solver option *cvt:names=2,3* provided.
- Added option *tech:stats* to return solution 
  statistics via problem suffixes.


## 20260109
- Remove excessive '_' symbols in output 
  model constraint and variable names.


## 20251213
- Recognize the logistic function 1/(1+e^(-x)),
  option *cvt:pre:logistic*.
- Recognize the signpow function sign(x)*x^k
  from several kinds of algebra, such as abs(x)*x^k.
  Option *cvt:pre:signpow*.
- *iis:find* now the main name for option
  *alg:iisfind* (which is now a synonym).
- Option *obj:multi:options*: fix integer-valued
  options.


## 20251121
- Meta-driver MP2NL, see
  https://mp.ampl.com/modeling-tools.html#mp2nl.
- New option *obj:multi:options* to control
  whether multi-objective option suffixes
  are used.
- Catch duplicate solver options.


## 20251021
- Option *cvt:expr:nlassign* adjusting expression
  inlining. Improves problem building performance
  for some nonlinear solvers on large models.


## 20251015
- Option *cvt:pre:ctx2count* is deprecated,
  with a new default value 7. Its former behavior
  is specialized and generalized by the new option
  *cvt:pre:ctx2bndeq*.
- Option *cvt:pre:ctx2bndeq* controls whether
  conditional comparisons with upper/lower
  bounds are reformulated in exact context,
  vs. always mixed. Default yes, value 1.
- Workaround for #270 (redefinition of
  an expression as an existing expression can
  lose connection to the original).
- Fixed *multiobj=2* for a single objective.
  Running multiobjective for a single-objective
  model makes sense with a negative `.objweight`
  suffix.
- Options *cvt:pre:ineq2bndeq*,
  *cvt:pre:ineq2related* control unification
  of conditional equlaities and inequalities.
- *cvt:compl:eps=1e-9* new default.
- *alg:start=3* default (new value),
  which passes both the
  initial guesses and basis for LP
  (previously basis only).
- Option *cvt:uenc:negctx:max* disregards lb(x), ub(x).
- Options *cvt:pre:ctx:abs* etc. control propagation
  of expression context.


## 20250814
- Fixed reformulation of numberof and alldiff.
  Previously could produce an inefficient
  reformulation.
- Option *cvt:pre:sort*:
  - Sorting and eliminating duplicates
    in arguments of AND, OR, MIN, MAX.
    This is necessary for some solvers.
  - Sorting arguments of COUNT, ATLEAST,
    EXACTLY, ATMOST, NUMBEROF, ALLDIFF.
- Option *cvt:pre:unnest* extended for MIN, MAX
  (bit 3).



## 20250811
- The legacy fmt::MemoryWriter has double output
  precision of 17 significant digits
  (previously 6).
- Option *cvt:pre:boundlogarg* (default 0) to bound
  arguments of logarithms nonnegative. Previously
  always done, sometimes deteriorating performance
  of nonlinear solvers.


## 20250801
- Tolerances set by options *cvt:pre:feastol*,
  *cvt:pre:feastolrel* both need to be violated
  to produce a warning on contradicting
  variable/constraint bounds. Previously
  the preprocessor failed on any violation,
  without letting the solver try.
- Options *cvt:compl*, *cvt:compl:eps* control
  complementarity reformulations.


## 20250617
- Option *cvt:unnest*: bits 1 and 2 switch on
  inlining of linear and quadratic subexpressions
  produced during reformulations (by default on).
- Options *cvt:pre:ctx2ineq*, *cvt:pre:ctx2count*
  to control context propagation into conditional
  comparisons #267.


## 20250616
- Multi-objective emulator: added support for 
  objective-specific options via objective suffixes
  beginning with *option_*.


## 20250429
- Fixed a bug in parsing of quadratic expressions,
  which could wrongly parse products of unequal
  linear expressions, such as (x-3)*(x-z-5).


## 20250424
- Option *alg:sens* returns synonym suffixes
  `.down/.up/.current` for objective coefficients
  and `.down/.up` for right-hand sides.
- Improved parsing speed of large sums of quadratic
  and polynomial expressions.
  - Setting option *cvt:qp2passes=0* switches to
    a simpler method, usually slightly slower.
  - Option *cvt:multoutcard* to limit the size of
    out-multiplied QP expressions. Can improve speed
    on large models, but prone to numerical issues.


## 20250308
- Fixed a bug in monotonicity presolve which could
  have led to wrong models using a^x with a<1.
- Reduced memory fragmentation for linear
  and quadratic expressions.
- Option *tech:outlev_mp=1* (implied by *tech:outlev=1*)
  prints initial and transformed model statistics.
- Improved presolve for disequality: more cases when
  no logical disjunction is needed.


## 20250205
- Fixed a bug in model reformulations which could
  cause some constraints to be lost, see issue 248.
- Option *alg:sens=1* now returns suffix `.sensobj` with
  current objective coefficients in the solver model
  (this is the model corresponding to AMPL command `solexpand`.)
- Option *acc:_all* is overridden by individual acceptance
  options, e.g., *acc:or*.
- Substitute AMPL defined variables
  into linear, quadratic, and polynomial expressions
  (option *cvt:dvelim*).
  This can simplify quadratic and polynomial models
  (linear substitutions are already performed by AMPL,
  see AMPL options *linelim* and *substout.)
- *cvt:prod=7* default for LP and convex solvers,
  logicalizing also products of just 2 binary variables.
- More presolve for logical expressions
  (Not, And, Or, Indicator), and new options to control
  some of them: *cvt:pre:ineqresult*, *cvt:pre:ineqrhs*.
- Fix lower bound calculation of the division result.
- [BREAKING] Option *acc:pow* now affects only expressions x^y
  with both x, y variable; previous meaning of *acc:pow*
  is now with *acc:powconstexp*.
- Option *tech:writemodel:index* to choose the iteration
  when solver model is exported
  in the multi-objective emulator.
- SCIP (and any solver with linear objective
	and non-linear constraints): improve reformulation
	of QP objectives.
- Fix reformulation of non-linear objective expressions
  for the multi-objective case (option *obj:multi*) when
  negative objective weights are used (*obj:multi:weight*).


## 20240724
- Option *acc:_all*
  - Useful to disable all reformulations (*acc:_all=2*),
    or force linearization (*acc:_all=0*).
- Faster input of quadratic expressions.


## 20240617
- *Multi-objective emulator*
  - *obj:multi=2* forces emulation, even if MO natively supported.
	- Fixed a bug in the objective degradation suffixes
    `.objabstol`, `.objreltol`.


## 20240604
- Presolve division by constant, resulting in fewer constraints
- Fix no-solution case in multi-objective emulator


## 20240529
- *Multi-objective emulator*
  - All flat MP solvers support multi-objective mode (*obj:multi=1*),
		either natively, or via emulation.
  - Suffixes `.objpriority`, `.objweight`, `.objabstol`, `.objreltol`.
  - [BREAKING] Default intuitive handling of `.objweight`,
    see option `obj:multi:weight`, even when natively supported.


## 20240429
- *Options `report_times` and `timing`*.
  - Options `timing` and `report_times` now have the same effect:
    setting their value to 1 returns basic timing information, 
    setting it to 2 returns more granular info.


## 20240320
- *SOS constraints*.
  - Fixed handling of SOS2 constraints created by AMPL
    as reformulations of PL expressions (`option
    pl_linearize 1`, default; set to 0 to use the solver's
    native PL functions if supported, or MP linearization.)
  - Disallow repeated weights for SOS constraints
    (suffixes `.sosno`/`.ref`).
- *Reformulation explorer*.
  - Upgraded option *writegraph* exports the reformulation
    graph which can be explored with the script in
    support/modelexplore (WIP).
- *Native handling of POW(x, INT)*.
  - Power expressions with positive integer exponent
    are passed natively to the solvers accepting them,
    vs previously quadratic or linear reformulation.
  - For best performance, global solving capability
    might be needed (e.g., Gurobi: `global=1`.)
- *Option `report_times`*.
- *Unused `acc:` options*.
  - The constraint acceptance options `acc:...`
    for non-handled constraints are ignored
    (previously triggered error.)
- *NLWPY*.
  - Python NL Writer API for MIQP models.
  - Available on PyPI.


## 20240115
- *Solve result codes*.
  - List codes by running (solver) -!
  - [BREAKING] Standardized codes. Major changes:
    - 100-199 (solved?) means solution candidate
      provided, but can be suboptimal/infeasible
    - 300-349 means unbounded problem but
      feasible solution returned
    - 400-449 means limit/interrupt but feasible
  - [BREAKING] *sol:chk:fail* returns code 150 (solved?)
- Improved translation of *SOCP constraints*.
  - Options *cvt:socp*, *cvt:socp2qc*.
- Compact solution check warnings
- Fixed presolve of the power function #226.
- Fixed graceful exit on Ctrl-C from AMPL in Linux.


## 20231103
- Improved translation of logical constraints:
  inlining of nested disjunctions and conjunctions;
  fewer auxiliary binary variables.


## 20231017
- Fixed a bug in NL reader on Windows.


## 20230919
- *mp_options*.
	Receive mp_options from AMPL (for all MP solvers).
  They are parsed before `(solvername)_options`.
- Solution checking: relative tolerance
  *sol:chk:feastolrel*; options *sol:chk:round*, *sol:chk:prec*.


## 20230831
- Solution checking, options *sol:chk:...* (experimental).
- Preprocess And/Or constraints.


## 20230817
- Alternative solutions: solve status equal to that
  of the final solution.
- Fixed a bug causing repeated names for
  auxiliary variables and constraints.
- Option values can be assigned without '='.
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
