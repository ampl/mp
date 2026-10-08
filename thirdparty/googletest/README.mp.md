# Vendored Google Test and Google Mock

Version: v1.18.0
Upstream: https://github.com/google/googletest
Commit: 063de7e9578f82b369302001269680b4b1553359
License: BSD-3-Clause; see LICENSE and upstream file notices.

The googletest/ and googlemock/ include and src directories are copied unchanged
from the pinned tag. MP compiles their amalgamated sources in test/CMakeLists.txt
so its tests use the same compiler and runtime as MP. No network access or
installed Google Test package is required when configuring an independent checkout.
The upstream tests, examples, and build scripts are not needed by this integration.
To update, replace these source directories from a reviewed upstream release,
retain the license, update this provenance, and validate all MP test configurations.
