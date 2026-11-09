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
