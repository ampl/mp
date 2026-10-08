# TASK-001: Unit test investigation

- Status: done (investigation and roadmap only); maintainer: Codex.
- Last meaningful update: 2026-10-08 (Australia/Sydney).
- Repository: MP; inspected HEAD: b92cc0e89aab2906f274cf0bd3c06b6e319d391e, detached.
- References: [task rules](../../../TASK_RULES.md), [specification](../../../specifications.md), [testing manual](../../../doc/source/testing.rst), [test configuration](../../../test/CMakeLists.txt).
- No issue, implementation commit, or PR is associated with this investigation. The investigation record and task index are delivered in the documentation commit containing this record; commit authorized on 2026-10-08.

## Goal, scope, and acceptance

Investigate the available Debug unit-test build, distinguish defects from environment/dependency limitations, and prepare a prioritized repair roadmap with validation steps. Completion requires configuration inspection, execution of available tests, and evidence with limitations. Behavior fixes, dependency replacement, publication, and release certification are excluded. Committing the investigation documentation was subsequently authorized. These investigation criteria are met; repairs remain follow-up work.

## Verified baseline and decisions

Read AGENTS.md, TASK_RULES.md, specifications.md, and test documentation/configuration. No existing task index or matching record was present. Preserve pre-existing untracked documentation and end-to-end artifacts.

The existing build uses Visual Studio 17 2022, MSVC 14.44.35207, CMake/CTest 3.31.2, BUILD_TESTS=ON, and empty BUILD (optional modules disabled). The assumed build/Debug directory does not exist: CMake expects executables under build/bin/Debug. All 20 C++ test executables were absent. Only the CMake test could run.

The existing Debug build fails in gtest-extra: bundled Google Test/Google Mock assumes std::tr1 is available on modern MSVC, producing C2039/C3083 and cascading errors. This is a test-framework/toolchain compatibility defect, not a solver license failure.

To investigate beyond that blocker without changing source or the original cache, created local build-unit-investigation with CMAKE_CXX_FLAGS=/DGTEST_USE_OWN_TR1_TUPLE=1. This is an experimental build setting, not a delivered fix. It successfully builds all 20 C++ test executables. No tracked production source or specification was changed.

## Validation actually executed

Commands below are reproducible from the repository root; generated logs are local and are not public deliverables.

```console
ctest --test-dir build -C Debug -N
ctest --test-dir build -C Debug --output-on-failure --timeout 60
cmake --build build --config Debug --parallel 4
cmake -S . -B build-unit-investigation -DBUILD_TESTS=ON -DBUILD_DOC=OFF -DBUILD_EXAMPLES=OFF -DCMAKE_CXX_FLAGS=/DGTEST_USE_OWN_TR1_TUPLE=1
cmake --build build-unit-investigation --config Debug --parallel 4
ctest --test-dir build-unit-investigation -C Debug -E converter --output-on-failure --timeout 60
ctest --test-dir build-unit-investigation -C Debug -R converter --output-on-failure --timeout 60
```

- Original CTest: 1 passed; 20 Not Run (missing executables). Original build: failed in test framework.
- Diagnostic build: successful. Across two disjoint CTest selections: 12 of 21 targets passed; 9 failed. Counts are targets, not individual Google Test cases.
- Passed: cmake-test, assert-test, clock-test, common-test, error-test, expr-writer-test, option-test, problem-builder-test, rstparser-test, safeint-test, converter-flat-test, converter-mip-test.
- Failed: util-test, expr-test, expr-visitor-test, nl-reader-test, os-test, problem-test, solver-test, sp-test, suffix-test.
- Both compulsory converter targets passed. They use test model APIs/mock backends; no licensed optimization engine was needed.
- A follow-up solver run from build-unit-investigation/test using ../bin/Debug/solver-test.exe --gtest_filter=-SolverTest.NSolSuffix completed: 105 cases, 83 passed, 22 failed. This exclusion diagnosed cases after the original abort; it is not an acceptable final validation filter.
- Isolated SuffixSetTest.AddDuplicateSuffix reproduced access violation 0xc0000005.

| Target | Observed evidence | Initial classification |
| --- | --- | --- |
| util-test | 1/5 failed: expected code 42, message contains literal {} | Confirmed diagnostic formatting defect; test helper calls Error with an integer second argument |
| expr-test | 3/88 failed: allocator hooks unused; duplicate-function message retains {} | Allocator expectations conflict with direct allocation in current implementation; Error overload defect also applies |
| expr-visitor-test | InvalidExpr fails on unhandled/unsupported dispatch expectations | Mock contract drift requires review |
| nl-reader-test | 78 cases fail, including offset 89 expected bound; header defaults and new builder callbacks differ; cascading mock failures/leak reports | Large fixture/mock drift plus any remaining reader defects; do not count each cascade as independent defect |
| os-test | 2/22 fail: unmapped-memory death tests report an exception; Unicode path test passes | Windows SEH/death-test handling incompatibility is the leading explanation |
| problem-test | 1/37 fails on exact LinearExpr capacity (actual initial 6, later 12) | Implementation-detail expectation drift, not evidence of incorrect algebra |
| solver-test | Original abort 0xc0000409 at NSolSuffix; countsolutions reads -858993664; filtered run reveals 22 failures | Confirmed unsafe option storage plus fixture/API drift and unsafe test iterator usage |
| sp-test | 1/39 fails: defined variable 0 not provided | Fixture omits current common-expression position/known-state metadata |
| suffix-test | 1/26 fails: duplicate suffix raises access violation | Reproducible crash; precise cause still needs debugger/sanitizer investigation |

## Prioritized repair roadmap (proposed, not implemented)

1. **Restore a reproducible default test build.** Scope the tuple workaround to the framework and its consumers, or modernize bundled Google Test/Mock with an explicit migration. The experiment proves the own-tuple switch works here. Validate a fresh default Windows C++17 build with no diagnostic global flags, check that all expected binaries exist, and run CTest. Review Unix compatibility and custom gtest-extra macros before framework replacement.
2. **Fix production safety and error-constructor ambiguity.** In src/solver.cc, BoolOption binds bool storage as int using *(int*)&value: reads beyond the bool and writes adjacent storage. Replace with type-correct access and regress default/set/reset for every affected bool option. In include/mp/error.h, Error(message, int exit_code) competes with formatted Error(format, integer): observed messages retain placeholders. Introduce an unambiguous distinction and audit integer-format/explicit-code call sites while preserving status mapping. Reproduce duplicate-suffix failure under a debugger/ASan and fix its established cause rather than weakening the test. Validate isolated regressions plus solver/suffix/expr/util targets, Debug and Release, and memory checks where supported.
3. **Make crash-prone tests diagnostic.** SolverTest.NSolSuffix dereferences find_if without checking for end; nsol is not registered by the inspected BasicSolver constructor. Determine the intended registration boundary and assert presence before dereference. Reconcile suffix registration and solution counting with current backend behavior. Run the complete solver target without exclusions after repair; retain regressions for the observed unsafe option storage.
4. **Repair NL fixtures and mock lifecycle/contracts.** include/mp/nl-header.h now defaults to binary and three AMPL options; text fixture builders must explicitly select TEXT and intended options rather than inherit defaults. Review ReadArithKind and missing-option contracts separately, add current AddVars/NotifyObjChoice callbacks, correct objective selection expectations, and isolate remaining parser failures after these changes. Validate text and binary fixtures, malformed headers, bounds, callbacks, reader/solver/problem-builder targets, and a public writer-to-reader round trip. Treat leak reports following mock exceptions as unresolved until isolated.
5. **Refresh implementation-dependent fixtures and platform tests.** Replace exact LinearExpr capacity assertions with meaningful growth/content invariants. Decide whether BasicExprFactory's documented custom allocator contract remains supported: current direct new char*[] plus delete[] through char* bypasses allocator hooks and deserves allocation/deallocation/alignment review, not just deletion of assertions. Set common-expression position/known metadata for SP fixtures. Update visitor expectations and timing/option parsing cases against current documentation. Investigate Windows SEH death-test capture using a targeted framework configuration without disabling safety assertions globally. Validate affected targets on Windows and a Unix toolchain, plus sanitizer checks for allocation changes.
6. **Close coverage gaps.** Add an MP core build-and-CTest CI job: inspected GitHub/Azure jobs build NL Writer wheels, not the full MP suite. Add explicit optional ASL and representative driver configurations when prerequisites are available; keep unavailable dependencies visible. Run representative end-to-end model tests and documented solution-checker configurations after production solver changes, using a separate fixture/output location to preserve local artifacts.

## Limitations and recovery

No ASL or optional solver-specific unit targets are registered in this configuration. NOSE2 is not found, so support-test is absent. CPLEX was reported unavailable during diagnostic configuration, but is not required by this core build. No SDK/license-dependent validation or end-to-end suite was run; their availability was not established. No runtime failure above was shown to originate from a missing solver license. Release, Unix, and sanitizer validation remain planned. No fixes are delivered and the original build remains blocked.

Local artifacts: build/unit-investigation-build.log, build/unit-investigation-configure.log, build/unit-investigation-tuple-build.log, build/unit-investigation-core-ctest.log, build/unit-investigation-converter-ctest.log, build/unit-investigation-suffix.log, build-unit-investigation/solver-investigation.log, and generated CTest logs. The diagnostic build tree is new and untracked; original generated artifacts were preserved. CTest temporary logs can be overwritten by later runs, so use the named logs for this baseline.

Next concrete operation for an authorized repair task: re-read this record, verify HEAD/working tree, implement a test-scoped framework compatibility fix, then address BoolOption storage and Error constructor ambiguity with regressions. Re-run the complete suite and investigate the suffix crash with debugger evidence. Update specifications only when an agreed repair changes documented public contracts.
