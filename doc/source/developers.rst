Developer documentation
-----------------------

.. toctree::
   :maxdepth: 2
   
   Library introduction <intro>
   howto
   API manual <components>
   testing
   details

Customer-facing release documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Solver READMEs explain current use; changelogs explain changes for customers.
When preparing release documentation:

* Describe the user-visible trigger and outcome of each fix or feature.
* Retain version requirements, beta or experimental status, limitations,
  option defaults, and controls that restore previous behavior.
* Check claims against supported behavior. A clearer error for an unsupported
  feature does not mean that feature is now supported, and build compatibility
  does not establish a shipping or performance guarantee.
* Keep local test results, implementation mechanisms, and investigation history
  in task records. Never link private records from public documentation.
* Put supported build and extension instructions in the developer manual.
  Internal-only changes do not require a customer changelog entry.

During development, update one concise customer-outcome note and maintain the
owning reference documentation. At release-candidate preparation or release PR
opening, apply the repository's
:download:`release-notes workflow <../../skills/mp-release-notes/SKILL.md>`.
Before tagging or packaging, confirm that review still covers the final diff;
repeat affected checks after late changes. Consolidate related changes, target
roughly two minutes of main-entry reading, and retain compatibility and
correctness caveats even when they require more space.

Compare final behavior with the previous shipped release. Omit intermediate
breakages fixed before shipping, but retain genuine shipped regressions and
intentional changes in defaults or explicit controls. Validate customer claims
against test semantics: a test that asserts placeholder zeros does not establish
meaningful per-objective reporting. Record a coverage map in the task/release
record rather than adding implementation detail to the changelog.
