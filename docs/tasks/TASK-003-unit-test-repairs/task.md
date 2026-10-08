# TASK-003: Sequential unit-test repairs

- Status: active; maintainer: Codex; updated: 2026-10-08 (Australia/Sydney).
- Related: [TASK-001](../TASK-001-unit-test-investigation/task.md), [TASK-002](../TASK-002-format-library-investigation/task.md), [specification](../../../specifications.md).
- Baseline: a2a3b5ae; branch: codex/unit-test-repairs.

## Scope and acceptance
Perform roadmap steps 1-4 sequentially, preferring a Google Test update. Commit each distinct logical change with verification evidence. Flag and document any necessary major production API change without implementing that step. Plan steps 5 and 6 in detail without implementing them. Preserve unrelated artifacts. Work in this chat; the completed format investigation is the only other chat involved.

## Progress and decisions
Read guidance, specification, both investigations, and current working tree. TASK-002 corrected the initial suffix crash classification: the diagnostic global CXX_FLAGS override removed /EHsc. A new normal build must preserve compiler defaults. Inherited uncommitted TASK-002 and TASK-001 correction will be delivered as their own documentation checkpoint before repair code.

## Validation and next step
No repair tests run yet. Next: pin and vendor modern Google Test/Mock, adapt test-only compatibility helpers, and rebuild with default MSVC exception flags. Error constructor redesign would alter public call meaning; retain signatures and investigate explicit preformatting of ambiguous internal calls instead.

## Delivery and recovery
No source changes yet. Baseline build logs and untracked artifacts remain in place. No push or PR requested. Update this record at every repair/validation checkpoint.

## Step 1: Modern test framework (complete)
Vendored unchanged Google Test/Mock v1.18.0 headers/sources at upstream commit 063de7e9578f82b369302001269680b4b1553359 with license and provenance. Removed old amalgamations and obsolete compiler workarounds. Updated helper assertions and optional ILOG CP source compilation; no MP public API change. Added a local-object unwinding regression to detect missing exception settings. Documentation explains offline configuration and compiler-default preservation.

Fresh build/unit-repairs with default MSVC flags builds all Debug tests and examples. Across disjoint CTest selections: 14/21 targets pass; 7 fail (util, expr, expr-visitor, nl-reader, problem, solver, sp). Both converter targets, all suffix cases, OS tests, and the unwinding regression pass. Solver still aborts at NSolSuffix. Logs: build/unit-repairs-step1-ctest.log and build/unit-repairs-step1-extra-ctest.log. Optional solver/ASL configurations unavailable and not validated. Step 1 delivered in the framework-update commit containing this checkpoint.

Execution constraint updated: all subsequent builds use --parallel 1; tests and commands execute serially; no additional workers. Earlier --parallel 4 build has completed. Keep detailed changes in commit messages and verification/deferrals here; a separate history file is not needed yet.

Next: type-correct bool options, then explicit internal Error formatting without changing public constructor signatures.
