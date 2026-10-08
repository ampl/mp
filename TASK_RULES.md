# Task tracking rules

## Purpose and automatic use

Task records preserve enough context to continue work after an interruption, context loss, a change of developer/agent, or a later session. They are working documents, not replacements for specifications, source history, or issue review.

Agents must create or reuse a canonical record automatically when the task has any of these characteristics:

- It is likely to span multiple sessions or need a handoff.
- It has several substantial stages of investigation, design, implementation, and validation whose results would be expensive to reconstruct.
- It coordinates changes across repositories or independently maintained components.
- It involves significant unresolved decisions, experiments, blockers, or delivery dependencies.
- The user explicitly requests tracking.

Use judgment based on recovery needs, not merely elapsed time, file count, or the length of the prompt. A small, straightforward edit or short read-only answer normally needs no record. If the task grows beyond that scope, start a record then. When uncertain, prefer a short record over a large reporting framework.

Creating and updating a record within the authorized task scope requires no separate user request. Briefly identify the record and its location when introducing it. This does not authorize committing, publishing, messaging others, or expanding the task.

## Scope and clarification

Start from the user's request, the referenced issue, applicable specifications, and available evidence. Record the goal, scope, and acceptance criteria before substantial implementation when practical.

Ask for clarification when an unresolved choice materially affects the desired result, compatibility, visibility, or acceptance. Do not ask routine questions whose answers are already available, or require the user to fill out a form before useful work can begin. Record provisional assumptions explicitly and continue work that does not depend on the unanswered choice. Do not treat an unanswered question as approval.

## Canonical record and visibility

Use one canonical task record. Search for an existing matching record before creating another; start with the task index when present.

For repository-owned tracking, use:

```text
docs/tasks/
  README.md
  GH-123-short-goal/
    task.md
    design.md       # optional
  TASK-001-short-goal/
    task.md
```

Include the issue's repository and URL in the record; issue numbers are not globally unique. Use a local `TASK-NNN` identifier when there is no issue and select an unused ID. Keep names stable as implementation evolves. Do not create empty indexes or supporting documents before there is a task to track.

The repository's visibility determines what may be stored here. Public task records, examples, reports, and links must be safe for public distribution. Private tracking must live in an explicitly designated private repository or other authorized private destination, not in an ignored file inside a public checkout.

If the user or an authorized project workflow supplies a canonical record outside this repository, maintain it there instead of creating a duplicate public record. Independent development must not require that external destination. If private tracking is needed and no suitable destination is known, clarify the destination before writing sensitive details; continue unaffected work.

A separate public issue/PR summary may describe deliverables safely, but must not duplicate the private progress log or expose its contents or location.

## Record contents

Keep `task.md` concise and current. Use the following fields or equivalent headings:

- **Identity and links:** stable ID, issue/repository, related tasks, applicable specifications, and PRs when they exist.
- **Owner and status:** responsible developer when known, current maintainer, lifecycle state, and last meaningful update with date/time zone when needed.
- **Goal and scope:** problem, intended outcome, exclusions, and provisional assumptions.
- **Acceptance criteria:** observable conditions for completion, including whether merge, release, or downstream integration is required.
- **Design and decisions:** selected approach, rationale, rejected alternatives when useful, and unresolved questions.
- **Progress and next step:** completed work, remaining steps, concrete blockers, and the next actionable operation.
- **Validation:** checks actually run, commands/configurations, results, skipped checks, and limitations.
- **Delivery and recovery:** affected repositories, branches or detached revisions, commits/PRs, important uncommitted work, necessary prerequisites, and artifact/report locations.

Separate verified facts, proposed plans, and unanswered questions. Do not invent an owner, deadline, approval, support guarantee, or test result.

Start with one file. Add `design.md` for substantial design reasoning; add other documents only when they improve usability. Link large reports/artifacts to authorized durable locations, noting local-only availability. Do not store secrets, confidential models, raw chat transcripts, or complete command logs.

## Lifecycle

Use `proposed`, `active`, `blocked`, `review`, `done`, or `cancelled`.

- `proposed`: scope or work is being considered.
- `active`: work is underway; record outstanding questions and incomplete stages.
- `blocked`: a concrete condition prevents the next necessary step; identify what would unblock it.
- `review`: identify the review or acceptance still required.
- `done`: the stated acceptance criteria are met, with delivery and validation evidence.
- `cancelled`: work was stopped; record why and the disposition of partial changes.

An interruption is not completion. Preserve the last accurate state and a recovery checkpoint; verify reality on resumption. Implementation, component commit, integration, merge, and release are distinct milestones and only required when included in the agreed scope.

## Updates and disruption recovery

Update the record at meaningful checkpoints: initial scope, consequential findings, design/scope changes, completed stages, blockers, validation results, and delivery changes. Do not log every command or make speculative claims about work still running.

Before a planned handoff or session boundary, save:

1. The current objective, scope, assumptions, and outstanding user decisions.
2. Completed changes and the actual repository/commit or uncommitted-file state.
3. Executed checks and their results, separately from planned checks.
4. The next concrete step and prerequisites needed to perform it.
5. Relevant artifact locations and any local-only or temporary limitations.

Unexpected interruption may prevent a final update, so checkpoint after substantial findings or changes rather than only at the end. A record can be maintained while work continues; it does not justify ending authorized work early.

On resumption, read the applicable instructions, specification, and canonical task record. Verify branch/HEAD, working tree, outstanding operations, and delivery state before acting. Treat the record as a recovery aid, not as proof that files, tests, or external state still match its last update.

## Ownership, index, and closure

When tasks exist, `docs/tasks/README.md` explains organization and links to active records with their owners/status. Keep it small and synchronize it when lifecycle states change. Retain completed records at stable paths so historical links remain valid; avoid maintaining a second full progress summary in the index.

One maintainer should coordinate updates to a shared task record. Contributors record their own results without overwriting unrelated progress; resolve conflicting changes explicitly.

Keep lasting requirements, compatibility commitments, and significant architecture in `specifications.md`. Task records explain the work that led to a change. Update the owning specification and relevant tests when the authorized change affects those contracts.

Before closure, record delivered revisions, acceptance evidence, remaining limitations, and any follow-up tasks. Follow-up ideas are not automatically authorized implementation work. Code and documentation remain in their owning repositories, and task-record updates follow the existing commit/review workflow.
