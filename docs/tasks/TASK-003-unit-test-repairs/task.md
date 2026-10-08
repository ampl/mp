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
