---
name: mp-release-notes
description: Consolidate MP and driver changelogs for a release candidate or release PR, or when explicitly asked to review release notes. Reconcile final shipped behavior, customer impact, migration guidance and public documentation; ordinary implementation commits need only a short outcome note.
---

# Consolidate MP release notes

Read AGENTS.md, specifications.md and the affected changelogs/reference pages
from this independent checkout. Use the authorized canonical record; keep private
investigation evidence out of public documents. This workflow does not authorize
publishing, tagging, packaging, or changing solver behavior.

## Establish the release boundary

Identify the previous shipped release and the target release/candidate from
release records, tags and the user's information; do not infer shipping from a
commit date or a reused version string. Compare final behavior to that baseline.
If shipping status is consequential and unresolved, flag it while reviewing
unaffected material. Preserve dates and published boundaries unless an editorial
correction is authorized. Highlight shipped regressions and affected versions;
omit intermediate breakages fixed within the same release. Retain intentional
changes to defaults and explicit controls even when default behavior is restored.

## Consolidate by customer outcome

Review shared MP and all affected driver entries together. Lead with customer
action/compatibility and correctness, then capabilities and reporting. Merge
successive changes to one feature. State the trigger, affected users and outcome;
keep option names/defaults/migration controls when they change customer action.
Move full value tables, bit masks and algorithm descriptions into maintained
public documentation and link to them. Internal-only changes need no release
bullet. Use a soft roughly two-minute reading target; never remove material
limitations or recovery guidance to meet a word count.

Check conditions that change applicability: AMPL versus driver options,
nondefault input paths, native versus emulated multiobjective solving, solver
capabilities/versions, local versus global results, and unavailable statistics.
Do not turn placeholders, assertions or build compatibility into feature claims.
For blended multiobjective passes, distinguish weighted pass certificates from
individual-objective bounds and gaps.

## Check coverage and hand off

Keep a concise map of customer-visible changes to the final bullets or their
public reference destinations in the authorized record. Verify defaults and
compatibility against source and meaningful test assertions; flag contract or
test gaps rather than hiding them through wording. Check links, formatting and
available documentation rendering. Report checks actually run and limitations;
solver reruns are needed only when behavior changes or unresolved claims require
them. Report entry size, omitted/moved detail and outstanding release questions.

Run this review at release-candidate preparation or release PR opening. Before
tagging/packaging, verify that it covers the final diff; repeat the affected
portion after late changes. AGENTS.md provides the explicit reminder even when
a client does not automatically discover this repository-owned skill.
