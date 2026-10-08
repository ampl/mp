MP component tests use Google Test and Google Mock v1.18.0, vendored in
``thirdparty/googletest``. See that directory's ``README.mp.md`` for the exact
upstream revision and license. The framework is compiled with MP's toolchain
and runtime; configuration does not download dependencies.

From an independent checkout, build and run the registered tests with::

    cmake -S . -B build-tests -DBUILD_TESTS=ON -DBUILD_DOC=OFF
    cmake --build build-tests --config Debug
    ctest --test-dir build-tests -C Debug --output-on-failure

Use ``ctest --test-dir build-tests -C Debug -N`` to inspect the test inventory.
Optional modules change which tests are registered. End-to-end driver tests
are separate; see ``doc/source/testing.rst``. Compiler exception flags must
remain enabled: do not replace CMake's default ``CMAKE_CXX_FLAGS`` with a
single diagnostic macro on MSVC.
