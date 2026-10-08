# Agent instructions

## Specifications

- Before changing behavior, read this repository's `specifications.md` if it exists and the relevant documentation. If it is missing, use existing code and documentation and report consequential uncertainty; do not invent requirements.
- Keep requirements and architectural contracts in the specification, and working instructions in this file.
- Update the specification when an authorized change affects documented behavior, public interfaces, compatibility guarantees, or significant design decisions.
- Keep documentation and development instructions self-contained for an independent checkout. Do not reference private consuming repositories or require access to their specifications or infrastructure.

## Changes and validation

- Preserve unrelated local changes and generated artifacts. Follow existing conventions and keep changes scoped to the task.
- Consider model semantics, numerical behavior, solver capabilities, options, suffixes, and result mapping when changing solver interfaces or transformations.
- Use this repository's build and test configuration. Run relevant unit, regression, and end-to-end checks for the affected behavior; a consuming project's build does not establish that component tests passed.
- Add regression coverage for substantive correctness changes when practical. Do not assume that all solver SDKs or licenses are available.
- Report what changed, what was checked, and any validation limitations.
- Keep credentials, confidential models, and private consumer details out of documentation and test fixtures.
