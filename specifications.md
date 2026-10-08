# MP specifications

Documentation baseline: 2026-10-08. Initial specification of the checked-out public project. **Current** statements describe inspected source/documentation; **proposed** statements describe acceptance or maintenance policy for review. Existing documented release checks are identified explicitly. Open decisions are in section 8; no new blanket compatibility or performance guarantee is implied.

## 1. Project overview

MP is public infrastructure for implementing AMPL solver interfaces and exchanging optimization models/results. It provides model reading and representation, solver-driver abstractions, model transformations, result mapping, and NL writing/SOL reading facilities.

Consumers include solver-interface developers, maintainers of MP-based drivers, modeling applications using the writer APIs, and users of resulting solver executables. The value is reusable infrastructure that lets drivers accept AMPL models and communicate consistent results while adapting to different solver capabilities.

Scope includes the core library, transformations, drivers, examples, NL Writer components, documentation, and tests. Optimization engines remain external dependencies.

## 2. Requirements

### 2.1 Functional requirements

**Current contracts and capabilities:**

- Read NL models through reusable reader interfaces and supported problem/model-manager representations; report parsing failures through the corresponding error mechanisms.
- Provide backend and model-API abstractions for constructing solver drivers. The documented recommended setup uses the backend hierarchy with flat or expression model APIs; legacy setups also exist.
- Represent expressions and constraints, negotiate supported constructs, and transform models when needed by the backend and selected options.
- Expose driver options and suffix mechanisms and map solver results into AMPL solution/status output. Particular features depend on the driver and solver.
- Provide NL Writer facilities and SOL reading for applications that submit models to AMPL-compatible solvers; these form a separate usable API surface.
- Provide examples and a mock `visitor` driver for understanding model traversal and starting a new driver.

**Proposed acceptance contracts:**

- Supported transformations preserve the documented meaning of the original model, including declared approximations, restrictions, and numerical tolerances. Unsupported constructs must not be silently treated as equivalent supported constructs.
- Postsolve mappings must return applicable values and status information in the original model's coordinates. Optional outputs such as duals, bases, or certificates require feature-specific contracts; availability is not universal.
- Errors, infeasibility, unboundedness, limits, and successful solution outcomes must be distinguished consistently with the driver documentation.

Human setup, build, test, and runtime procedures are in section 2.8 and the linked developer manual.

### 2.2 Data requirements

Inputs include NL models, optional name/auxiliary files, suffixes, options, and programmatic model data. NL supports a wide range of continuous/discrete, linear/nonlinear, and logical constructs; this does not imply support by every backend.

Outputs include solver models, SOL results and suffixes, optional diagnostic exports, and test reports. NL Writer accepts data through its APIs without requiring callers to adopt the main MP intermediate representation. Use the headers and component documentation as the authority for format/API details.

End-to-end fixtures consist of models/data/scripts or NL files and `modellist.json` descriptions containing feature tags, options, expected objective/values, and optional output checks. Missing or invalid fixture metadata should be investigated separately from solver correctness. Numerical validation uses the case's bounds/tolerances, not exact floating-point equality by default.

### 2.3 Non-functional requirements and compatibility

**Current:** the core build requires C++17. Optional driver dependencies and build choices constrain portability. CMake supports standalone and embedded use; generated binaries/libraries are placed in configured `bin/` and `lib/` directories.

**Proposed:** maintain documented format and interface behavior and use reproducible cases for numerical and performance regressions. Record conversion/runtime/memory separately where useful. Do not imply exact repeatability across solvers, thread settings, or floating-point platforms. Stable API/ABI scope, supported compiler/platform versions, and deprecation timelines require explicit agreement; see section 8.

### 2.4 Operational and quality requirements

| Validation layer | Source | Responsibility |
| --- | --- | --- |
| Unit and converter tests | [test/CMakeLists.txt](test/CMakeLists.txt), [testing manual](doc/source/testing.rst) | Library behavior, model representation, conversions, readers/writers, and supporting utilities. |
| End-to-end model tests | [test/end2end/run.py](test/end2end/run.py), [cases](test/end2end/cases) | Driver features, reformulations, numerical output, statuses, suffixes, and options. Requires chosen binaries and applicable AMPL/solver licenses. |
| Documentation/examples | [doc](doc), [examples](examples), [NL Writer](nl-writer2) | Public API usage and independently usable entry points. |
| CI/package jobs | [azure-pipelines.yml](azure-pipelines.yml), [.github/workflows/build-and-test.yaml](.github/workflows/build-and-test.yaml) | Inspect selected jobs for actual coverage. The inspected GitHub workflow builds NL Writer Python wheels; it is not evidence of full MP unit/driver testing. |

The existing [testing manual](doc/source/testing.rst) calls `converter-flat-test` and `converter-mip-test` compulsory and requires representative solvers and solution-checker configurations for release. It supplies `chk:fail` and `chk:feastol=1e-2` as a release-check configuration; this setting is not a universal accuracy guarantee for every model.

**Proposed:** correctness fixes include a representative regression case, important new features receive unit and/or end-to-end coverage, and reports distinguish passes from unsupported cases and skips. Use feature tags to select applicable tests without concealing regressions in features a driver claims to support. Sanitizer and performance checks should match the affected code and available toolchain.

### 2.5 External integrations

| Integration | Purpose and constraints |
| --- | --- |
| ASL | Optional adapters and model evaluation support; build selection and available ASL targets determine use. Independent public dependency. |
| Solver SDKs/runtimes | Backend optimization and feature-specific APIs; versions, headers, licensing, and supported platforms vary by driver. |
| Third-party source | Supporting libraries/test infrastructure; inspect build definitions and license notices before replacing or redistributing. |
| AMPL, Python tooling | End-to-end model execution and reports; not required for every core-library consumer. |
| Documentation and wheel tools | Optional documentation and NL Writer Python packaging; dependencies are maintained with those workflows. |

Submodule identities are in [.gitmodules](.gitmodules). A public standalone checkout must not require private consuming repositories or their infrastructure.

### 2.6 Technology constraints

[CMakeLists.txt](CMakeLists.txt) defines C++17, module selection through `BUILD`, examples, documentation, unit tests, library linkage, and optional sanitizer/profiler settings. [solvers/CMakeLists.txt](solvers/CMakeLists.txt) and driver sources define SDK discovery and driver-specific constraints. Version requirements should be maintained there rather than copied as a second changing inventory here. The minimum CMake declaration alone does not establish compatibility for every optional module.

### 2.7 Explicit non-goals

- Implement vendor optimization algorithms or guarantee the same feature set in every backend.
- Guarantee that every NL construct can be transformed for every solver.
- Define downstream licensing, deployment, or distribution policies.
- Require confidential models, user accounts, or a hosted service for independent use.

### 2.8 Human manual: standalone build, tests, and usage

These are source-derived reference commands, not execution results from preparing this document. Run from an independent MP checkout with Git, CMake, and a C++17 toolchain available. Initialize required submodules in a fresh checkout; preserve existing local changes.

```console
git submodule update --init --recursive
cmake -S . -B build-spec -DBUILD_TESTS=ON -DBUILD_DOC=OFF -DBUILD_EXAMPLES=ON
cmake --build build-spec --config Release
ctest --test-dir build-spec -C Release -N
ctest --test-dir build-spec -C Release --output-on-failure
```

Choose `-DCMAKE_BUILD_TYPE=Release` for a single-configuration generator when desired. `-N` shows what tests are registered. Optional module selection and ASL availability can change that list.

To build the documented mock driver, use a separate build tree:

```console
cmake -S . -B build-visitor-spec -DBUILD=visitor -DBUILD_TESTS=ON -DBUILD_DOC=OFF
cmake --build build-visitor-spec --config Release
```

Run the generated `visitor` executable on a model NL file as described in [the driver HOWTO](doc/source/howto.rst). It traverses the model; it is not an optimization engine. Binary locations can include the configuration name on multi-configuration generators.

For a real driver, select its documented build configuration and install the corresponding SDK/runtime and license. With AMPL and the chosen driver on the search path, a typical AMPL session is:

```ampl
model example.mod;
data example.dat;
option solver gurobi;
solve;
```

The model, data, and solver name are examples; use installed drivers and inputs compatible with their capabilities. See the selected driver's README under `solvers/` for options.

End-to-end tests are separate from CTest:

```console
python test/end2end/run.py gurobi --bin_path build-spec/bin --dir test/end2end/cases/categorized/fast
```

Install the needed Python dependencies from [test/end2end/requirements.txt](test/end2end/requirements.txt); `multiprocessing` is standard-library functionality rather than a required PyPI package. Adjust the binary path to the actual build. The [testing manual](doc/source/testing.rst) explains `--ampl`, reports, test subsets, feature tags, and release configurations. [NL Writer README](nl-writer2/README.md) covers its independent library/examples; [Python README](nl-writer2/nlwpy/README.md) covers its wrapper.

## 3. Design choices

### 3.1 Project structure and ownership

| Location | Responsibility |
| --- | --- |
| `include/mp/`, `src/` | Public headers, expressions, model representations, solver abstractions, I/O, and utilities |
| `include/mp/flat/`, related converter/backend headers | Flat representations, capability-dependent transformations, and backend integration |
| `src/asl/` | Optional ASL adapters |
| `solvers/` | Concrete backend/driver implementations, usage READMEs, and driver changelogs |
| `nl-writer2/` | NL Writer library, C/C++ examples, and Python wrapper |
| `test/`, `test/end2end/` | Component tests and model-based driver tests |
| `doc/source/`, `examples/` | Public developer/user documentation and examples |
| `support/`, root CMake/CI files | Build configuration and maintenance tooling |
| `thirdparty/` | Bundled dependencies and submodules, with their own licenses |
| Build trees and generated reports | Generated artifacts; not authoritative API/requirement sources |

### 3.2 Interfaces and dependency direction

Driver applications use model managers, converters, backend classes, and model APIs to translate models and communicate with solver SDKs. Capability-specific backend code depends on reusable abstractions; shared transformations should not acquire an unrelated vendor dependency. The NL Writer API can be used separately from the driver stack. ASL integration is optional and belongs to its adapter boundary.

Public specifications, documentation, and build instructions remain self-contained. Consumers reference these contracts and own their integration requirements; consumer-specific configuration does not redefine MP's public guarantees.

### 3.3 Important decisions

- **Backend/model-API separation (documented):** recommended driver architecture separates common driver services from solver model APIs. Consequence: shared behavior can be reused, while backend capabilities remain explicit.
- **Capability-dependent transformations (current):** MP can reformulate unsupported native constructs where a supported conversion exists. Consequence: conversion correctness and numerical assumptions need tests independently of solver convergence.
- **Result conversion (current):** model changes require mappings back to original entities. Consequence: objective values alone do not validate all postsolve behavior.
- **Separate reader/writer APIs (documented):** model interchange can be reused without adopting the entire driver architecture. Consequence: examples and tests should cover these entry points independently.
- **Legacy driver setups (documented):** older setups coexist with the recommended architecture. The HOWTO notes that support may be discontinued; no date is established here.

### 3.4 Sources of truth and related documents

| Information | Source |
| --- | --- |
| Public API/contracts | [include/mp](include/mp), [component manual](doc/source/components.rst), [NL Writer](nl-writer2/README.md) |
| Driver architecture and extension procedure | [howto.rst](doc/source/howto.rst), [developers.rst](doc/source/developers.rst) |
| Features and modeling guidance | [features-guide.rst](doc/source/features-guide.rst), [model-guide.rst](doc/source/model-guide.rst), driver READMEs |
| Build options/library version | [CMakeLists.txt](CMakeLists.txt), [support/cmake](support/cmake) |
| Driver SDK/capability/version definitions | [solvers](solvers) and [solvers/CMakeLists.txt](solvers/CMakeLists.txt) |
| Tests and release checks | [testing.rst](doc/source/testing.rst), [test](test) |
| Dependencies and licenses | [.gitmodules](.gitmodules), [thirdparty](thirdparty), [LICENSE.rst](LICENSE.rst) |
| Change history/agent workflow | [CHANGES.mp.md](CHANGES.mp.md), [AGENTS.md](AGENTS.md) |
| Task lifecycle and recovery records | [TASK_RULES.md](TASK_RULES.md); repository-owned records under `docs/tasks/` when needed |

### 3.5 Change coordination

Update affected public documentation, this specification, and tests alongside changes to contracts or transformations. Check every affected backend when changing a shared abstraction, recording unavailable SDKs/licenses. Coordinate real ASL adapter changes with that dependency's interfaces. Consumers decide how to pin and validate new revisions; independent MP validation must remain possible without their infrastructure.

## 4. Security and privacy

MP is public source with license terms and third-party notices maintained in the repository. Vendor SDKs may impose separate conditions. Runtime model input and diagnostic exports can contain confidential data even though the library is public.

Parsers, model conversions, user-defined function integration where enabled, and native solver calls are security-relevant input boundaries. **Proposed:** validate input and failures appropriately, exercise malformed-input and memory-safety cases where practical, use public/reproducible fixtures, and exclude credentials or confidential consumer data from reports and source. No application-level identity service or user-model persistence is provided by the core project. A full security audit is outside this documentation baseline.

## 5. Analytics and success criteria

Existing release expectations are in [testing.rst](doc/source/testing.rst). **Proposed acceptance:** applicable component/converter checks and representative driver tests pass; transformed solutions satisfy the specified checks in original-model terms; options, statuses, and suffixes match documented behavior; and public examples remain usable independently.

Record test counts, feature coverage, skips, numerical mismatches, and revision/configuration information. For performance work, record models, versions, seeds, threads, conversion/runtime/memory baselines, and target thresholds before comparing results. Global API/ABI guarantees, coverage percentages, and performance budgets are open decisions. User telemetry and adoption tracking are not required for core functionality.

## 6. Milestones

Proposed ongoing stages; no dates or completion claims are implied.

| Stage | Deliverable | Acceptance/dependency |
| --- | --- | --- |
| Review baseline | Confirm contract boundaries and open policies | Maintainer review |
| Implement capability or fix | Library/backend change, docs, regression cases | Necessary SDKs and reproducible fixtures |
| Validate public release | Converter/unit checks, representative end-to-end checks, examples | Existing release test requirements and recorded configuration results |
| Maintain compatibility | Changelog and documented migration where needed | Reviewed impact on affected APIs/backends/consumers |

## 7. Use cases

- **Create a driver:** start from the documented setup or `visitor`, declare backend capabilities, implement model/result interfaces, and register representative model tests.
- **Add a reformulation:** specify valid domains and numerical assumptions, implement the conversion and mappings, and test original-model correctness across relevant backends.
- **Use model interchange APIs:** supply a model through NL Writer or consume NL through reader APIs without requiring a vendor backend.
- **Investigate a regression:** reduce the model, capture options and versions, isolate conversion versus solver/result behavior, and add a durable fixture.

## 8. Open decisions and known gaps

- Define the intended stability scope for public C++/C APIs and ABI, supported toolchains/platforms, and deprecation policy.
- Confirm the release representative-solver set, ownership of skipped configurations, and report retention requirements.
- Establish numerical/performance baselines beyond case-specific expectations; avoid turning example tolerances into universal guarantees.
- Review Python test dependency metadata, including the standard-library `multiprocessing` entry.
- Document actual CI coverage by job; NL Writer wheel builds do not demonstrate full library/driver correctness.
- Confirm legacy-driver maintenance and migration plans; the current HOWTO does not give a retirement schedule.
