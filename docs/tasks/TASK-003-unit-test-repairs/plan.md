# Remaining unit-test and coverage plan

Proposed follow-up to [TASK-003](task.md), 2026-10-08. Steps 5 and 6 are plans, not delivered changes. Use one build worker and serial tests; keep at most this chat plus one investigation chat. Complete each logical change and its validation before starting the next, and commit each separately.

## Step 5: Resolve remaining fixtures and allocation contracts

### 5a. Small, implementation-independent fixture repairs

For `problem-test.LinearExpr`, replace the reversed capacity comparisons (0/1 >= capacity) and exact reserved-capacity expectation (10, currently 12) with capacity-at-least-size/request and growth/content invariants. Current initial capacity is 6. Check retained coefficients and variable indices across growth and clearing; do not prescribe a particular small-buffer strategy.

For `sp-test.CommonExprInExpectation`, supply the common-expression position/known-state metadata required by the current builder. Keep the expected stochastic expression and indexing assertions. Cover genuinely missing definitions separately so the repair does not hide an invalid-model error. Validate the complete problem and SP targets before committing.

### 5b. Expression allocation and dependent visitor tests

The two allocator tests and `expr-visitor-test.InvalidExpr` require a production design review first. `BasicExprFactory<Alloc>` currently bypasses allocator hooks with `new char*[sizeof(Impl) + extra_bytes]` and releases through `delete[] (char*)`; the pointer-array allocation also scales the requested byte count by pointer size. InvalidExpr injects corruption through the unused mock allocator, so its dispatch assertions cannot test their intended condition.

Audit the public allocator template contract, historical issue-174 workaround, alignment of every expression/function implementation, variable-size overflow checks, ownership, and exception cleanup. Prefer restoring the existing allocator contract if it can be done without a public API change. Pair allocation/deallocation using the same storage mechanism, with sufficient alignment and checked sizes. Do not just remove mock allocation expectations or keep a mismatched delete.

Before implementation, document whether stateful/custom allocators and over-aligned storage are supported by the existing interface. If repairing these requires changing the public allocator API or supported semantics, flag it for a separate migration decision. Add regressions for hook use, allocation failure, cleanup, representative variable-size nodes/functions, and alignment. Validate full expr and visitor targets on Windows and GCC/Clang, then run available ASan/UBSan configurations. Only after allocation is correct should the visitor fixture be reconciled with actual invalid-kind dispatch and assertions.

### 5c. Solver option and reporting contracts

The current solver target has 16 failing cases, grouped below. Inspect the public option documentation and implementation together; classify each discrepancy as obsolete fixture or current defect before editing expectations.

| Group | Cases | Proposed check |
| --- | --- | --- |
| Option registration/parsing | SolverOption, UnknownOption, ParseOptionRecovery, FormatOption, ErrorOnKeywordOptionValue, ParseOptionsHandlesOptionErrorsInParse, NoEchoOnErrors, OptionEcho, EQOption | Use valid nonempty option names and current canonical names/aliases; verify unknown-option rejection, recovery to subsequent valid options, assignment syntax, and error/echo ordering. Preserve unsuccessful-option status and diagnostics. |
| Numeric rendering and timing | DoubleOptionHelper, TimingOption, ReportInputTime | Distinguish numeric-value correctness from display precision. Review supported timing modes (2 is now accepted) and current `tech:timing` echo before updating assertions. Validate finite parsed values and meaningful timing output. |
| Application output | ReportError, SolverOptionsEchoedByDefault | Separate throwing/reporting API tests from process exit behavior. Assert documented exit status, error text, default echo behavior, and explicit echo suppression at the relevant boundary. |
| Signal handling | SignalHandler, SignalHandlerExitOnTwoSIGINTs | Verify current solver display-name initialization and install/restore handler lifecycle. Reproduce the second-SIGINT case in a separate child process and check actual termination/status. Do not globally disable death tests or accept missing termination. |

Modern-framework/default-exception builds already pass the OS and suffix targets, so the earlier Windows SEH concern does not justify changes there. Run individual solver groups while diagnosing, then the complete solver target without exclusions and the full component suite. If option/status/signal behavior needs a consequential public change, document it and stop that implementation pending a migration decision. Commit independent fixture and production fixes separately.

### 5d. Deferred Error API decision

Internal integer diagnostics now preformat explicitly; the public `(message, int exit_code)` constructor still keeps its existing meaning. A tagged exit-code API, named factory, or changed formatting overload selection would require a documented public migration. Audit downstream-independent examples and public call patterns first, specify coexistence/deprecation, then add compatibility tests. This is separate from fixture repair and from updating fmt.

## Step 6: Make coverage reproducible

Start core CI after step 5 produces an understood passing baseline, or explicitly track unresolved failures rather than making them silent skips. The inspected wheel workflows do not establish core MP coverage.

1. Add a standalone core configure/build/CTest workflow using vendored test dependencies, C++17, `BUILD_TESTS=ON`, `BUILD_DOC=OFF`, and examples. Exercise Windows/MSVC and Linux/GCC or Clang in Debug and Release. Preserve compiler defaults; use `--parallel 1` and `ctest -j 1`. Require the expected test inventory for each configuration and run both compulsory converter targets without filters. Preserve configure/build/CTest logs and failure output as artifacts.
2. Add a separate Linux sanitizer configuration for expression allocation, readers, suffixes, and converters. Document compiler/linker settings and coverage; investigate failures rather than weakening assertions. Keep sanitizer and normal Release results distinct.
3. Inventory optional ASL, Python support-test, and solver-specific test prerequisites. Add an explicit ASL configuration when its public dependencies are available. Report missing optional dependencies and the resulting test inventory. Keep licensed SDK jobs on appropriate runners and avoid pretending an unregistered target passed.
4. Add representative driver/model checks: use public fixtures with the mock visitor and available converter driver for interface/model traversal, then available real solvers for numeric/status/suffix behavior. Select cases by actual claimed capability, including solution counting and multiple objectives where supported. Run the documented release solution-checker configurations (`chk:fail`, `chk:feastol=1e-2`) where applicable; record per-case numerical expectations and unsupported features.
5. Run end-to-end checks from a copied fixture/output tree so generated NL/SOL/log files do not modify the source checkout. Record driver/compiler versions and license-dependent skips without storing credentials or confidential models. Keep the existing NL Writer wheel job independently useful.

Suggested commits: core CI baseline; sanitizer job; optional ASL coverage; public driver/end-to-end harness; licensed runner configuration only when available. Acceptance is an independent checkout reproducing registered tests, required converters passing, meaningful failures failing CI, and optional coverage being explicitly visible. A green wheel build or core-only job must not imply all driver/release checks passed.

## fmt investigation relationship

[TASK-002](../TASK-002-format-library-investigation/task.md) treats the locally patched fmt 3.0.1 modernization separately. It does not resolve Error's constructor selection or allocation/fixture drift. A future fmt task should first inventory exposed fmt types and local patches, establish compatibility adapters, and run core tests plus public examples before proposing a migration. Do not bundle that API migration into steps 5/6.
