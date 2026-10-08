# TASK-003: Sequential unit-test repairs

- Status: done (steps 1-4 repairs and steps 5-6 planning); maintainer: Codex; updated: 2026-10-08 (Australia/Sydney).
- Related: [TASK-001](../TASK-001-unit-test-investigation/task.md), [TASK-002](../TASK-002-format-library-investigation/task.md), [specification](../../../specifications.md).
- Baseline: a2a3b5ae; branch: codex/unit-test-repairs.

## Scope and acceptance
Perform roadmap steps 1-4 sequentially, preferring a Google Test update. Commit each distinct logical change with verification evidence. Flag and document any necessary major production API change without implementing that step. Plan steps 5 and 6 in detail without implementing them. Preserve unrelated artifacts. Work in this chat; the completed format investigation is the only other chat involved.

## Progress and decisions
Read guidance, specification, both investigations, and current working tree. TASK-002 corrected the initial suffix crash classification: the diagnostic global CXX_FLAGS override removed /EHsc. A new normal build must preserve compiler defaults. Inherited uncommitted TASK-002 and TASK-001 correction will be delivered as their own documentation checkpoint before repair code.

## Validation and next step
Steps 1-4 are delivered; the final full-suite run passes 16/21 CTest targets with default MSVC exception flags. Error constructor redesign is deferred because it would alter public call meaning; existing signatures are retained. Detailed proposed follow-up is in [plan.md](plan.md). Steps 5-6 are not implemented.

## Delivery and recovery
Delivered commits: 2ce534b8 (baseline/tracking correction), 7f3fb008 (framework), d529c8de (bool storage), a7a5ea7d (diagnostics), 4044c81d (solution counting), cec0c376 (NL fixtures), followed by the documentation commit containing final evidence and step 5/6 plans. Baseline logs and unrelated untracked artifacts remain in place. No push or PR requested. Descriptive commits and this record provide history; a duplicate history file is unnecessary.

## Step 1: Modern test framework (complete)
Vendored unchanged Google Test/Mock v1.18.0 headers/sources at upstream commit 063de7e9578f82b369302001269680b4b1553359 with license and provenance. Removed old amalgamations and obsolete compiler workarounds. Updated helper assertions and optional ILOG CP source compilation; no MP public API change. Added a local-object unwinding regression to detect missing exception settings. Documentation explains offline configuration and compiler-default preservation.

Fresh build/unit-repairs with default MSVC flags builds all Debug tests and examples. Across disjoint CTest selections: 14/21 targets pass; 7 fail (util, expr, expr-visitor, nl-reader, problem, solver, sp). Both converter targets, all suffix cases, OS tests, and the unwinding regression pass. Solver still aborts at NSolSuffix. Logs: build/unit-repairs-step1-ctest.log and build/unit-repairs-step1-extra-ctest.log. Optional solver/ASL configurations unavailable and not validated. Step 1 delivered in the framework-update commit containing this checkpoint.

Execution constraint updated: all subsequent builds use --parallel 1; tests and commands execute serially; no additional workers. Earlier --parallel 4 build has completed. Keep detailed changes in commit messages and verification/deferrals here; a separate history file is not needed yet.

Next: type-correct bool options, then explicit internal Error formatting without changing public constructor signatures.

## Step 2a: Type-correct boolean option storage (complete)
Step 1 commit: 7f3fb008. Replaced BoolOption's reinterpretation of bool as int with a TypedSolverOption<int> adapter that reads/writes the actual bool and accepts only 0/1. Public option signatures, names, values and storage layout are unchanged. Added regressions covering defaults, reset, invalid values, independent bool values, and adjacent debug state. Single-worker Debug solver-test build succeeds; CountSolutionsOption and BooleanOptionsPreserveIndependentValues both pass. Logs: build/unit-repairs-bool-build.log and build/unit-repairs-bool-tests.log. Delivered in the bool-storage commit containing this checkpoint.

Next: explicit formatting at internal integer Error call sites; public Error constructor redesign deferred as a consequential API change.

## Step 2b: Diagnostic formatting (complete)
Bool-storage fix commit: d529c8de. Preformatted integer diagnostics in function definition, undefined-function handling, option range checking, and test process helpers. Added explicit-message/code compatibility coverage and documented constructor selection in the header. A global Error constructor redesign is deliberately not implemented: changing the meaning of existing (message, int) calls requires a public API migration decision. Existing explicit exit-code calls retain their meaning. Build and test results will be recorded before committing this logical change.

Step 2b validation: single-worker Debug build succeeds for error, util, expr and NL-reader targets. Serial error-test and util-test pass; ExprFactoryTest.AddFunctions passes with the numeric function index in its diagnostic. Explicit exit-code regression passes. Logs: build/unit-repairs-error-build.log and build/unit-repairs-error-expr-tests.log. No public Error signature or exit-code interpretation was changed. Delivered in the internal-diagnostic formatting commit containing this checkpoint. Next: guard the missing-suffix test and reconcile its current registration boundary.

## Step 3: Multiple-solution tests (complete)
Step 2b commit: a7a5ea7d. Current SolutionWriter creates output nsol/npool only with primal values; BasicSolver capability registration exposes options, not entries in deprecated Solver's suffix list. Replaced the unsafe obsolete registration-list dereference with a capability/option test. Counting fixtures now contain a variable, an objective and primal values; the matcher checks both problem/objective nsol and npool safely before reading values. Added the no-primal-output case and documented the existing writer boundary in specifications. Production suffix behavior and interfaces are unchanged.

Step 3 validation: complete solver-test runs to completion without exclusions or aborts: 108 cases, 89 pass, 19 fail. All multiple-solution option/writer cases pass, including both suffix scopes and no-primal handling. Remaining failures concern historical option/reporting/signal and NL callback expectations; these are explicitly not claimed fixed. Logs: build/unit-repairs-step3-build.log and build/unit-repairs-step3-ctest.log. Delivered in the multiple-solution fixture commit containing this checkpoint. Next: repair explicit text/binary NL fixture headers and current builder callbacks/objective selection, preserving parser semantics.

## Step 4: NL fixture contracts (complete)
Step 3 commit: 4044c81d. Text fixture serialization now explicitly chooses TEXT, while binary header fixtures explicitly choose their flags. Optional-header parser tests explicitly initialize their option fields rather than relying on obsolete defaults. Strict builder fixtures expect current zero-size variable blocks and objective-choice notifications. Multi-objective tests enable multiobj deliberately; single-objective selection expects one resulting objective and its selection callback. Bounds-order mocks explicitly request objective input. Bound diagnostic expectations retain existing max_digits10 output precision. Added a public text/binary header-and-bounds round trip that checks actual numeric bounds. No reader behavior or public API is changed.

First serial NL-reader run after fixture repair: 107/110 cases pass; three remaining precision/objective-interest expectations corrected in the next iteration. Logs retained at build/unit-repairs-step4-nl-ctest.log.

Final step 4 validation: single-worker Debug reader/solver build succeeds. All 111 NL-reader cases pass, including text/binary bounds round trips. Solver completes all 108 cases: 92 pass, 16 fail; the three objective-selection/header callback failures are fixed. Logs: build/unit-repairs-step4-build.log and build/unit-repairs-step4-ctest.log. Next: refresh every target and run all CTest targets serially, then deliver detailed step 5/6 plans.

## Final acceptance evidence and follow-up

The complete Debug build succeeds with `cmake --build build/unit-repairs --config Debug --parallel 1`. All registered targets run without exclusions using `ctest --test-dir build/unit-repairs -C Debug -j 1 --output-on-failure --timeout 60`: 16/21 pass in 53.08 seconds. Both compulsory converter targets pass, as do NL-reader, suffix, OS, error, util and problem-builder targets. No target is missing or times out. Local logs: build/unit-repairs-final-build.log and build/unit-repairs-final-ctest.log (not committed).

Remaining failures: expr-test (2 allocator cases), expr-visitor-test (InvalidExpr), problem-test (LinearExpr capacity), solver-test (16 cases), sp-test (CommonExprInExpectation). These 21 individual failing cases across five targets are recorded and planned in [plan.md](plan.md); the full suite is not green. Plans cover allocator alignment/ownership/API review, meaningful fixture assertions, option/status/signal contracts, core CI, sanitizers, optional dependencies and separate end-to-end outputs.

Configuration: MSVC 14.44.35207, C++17, tests/examples enabled, optional BUILD empty, default /EHsc preserved. Release, Unix, sanitizers, optional ASL/solver SDK tests and licensed end-to-end validation were not run. NOSE2 is unavailable, so Python support-test is not registered. Existing unrelated generated files and the diagnostic build tree are preserved. All work after the worker-count correction used one build worker and serial commands/tests; no additional agents were launched.

Authorized scope is complete: safe steps 1-4 committed independently, consequential Error API change flagged and deferred, and steps 5-6 discussed/planned without implementation. Next authorized follow-up would begin with step 5a in plan.md, then allocation review before visitor repair. fmt modernization remains the separate TASK-002 proposal.
