# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import build_ext
from setuptools import setup, Extension
import pybind11
import glob

__version__ = "0.0.1b0"

# The main interface is through Pybind11Extension.
# * You can add cxx_std=11/14/17, and then build_ext can be removed.
# * You can set include_pybind11=false to add the include directory yourself,
#   say from a submodule.
#
# Note:
#   Sort input source files if you glob sources to ensure bit-for-bit
#   reproducible builds (https://github.com/pybind/python_example/pull/53)


def compile_args():
    from platform import system

    if system() == "Windows":
        return ["/std:c++17"]
    elif system() == "Linux":
        ignore_warnings = [
            "-Wno-stringop-truncation",
            "-Wno-catch-value",
            "-Wno-unused-variable",
        ]
        return ["-std=c++17"] + ignore_warnings
    elif system() == "Darwin":
        ignore_warnings = [
            "-Wno-unused-variable",
        ]
        return [
            "-std=c++17",
            "-mmacosx-version-min=10.15",
        ] + ignore_warnings
    else:
        return []


ext_modules = [
    Extension(
        "nlwpy",
        ["nlwpy/src/nlw_bindings.cc"] + glob.glob("./src/" + "*.cc"),
        extra_compile_args=compile_args(),
        include_dirs=["include", pybind11.get_include()],
        # Example: passing in the version to the compiled code
        define_macros=[("VERSION_INFO", __version__)],
    ),
]

setup(
    name="nlwpy",
    version=__version__,
    author="Gleb Belov",
    author_email="gleb@ampl.com",
    url="https://github.com/ampl/mp",
    description="Python API for the AMPL NL Writer library",
    long_description="",
    ext_modules=ext_modules,
    extras_require={"test": "pytest"},
    # Currently, build_ext only provides an optional "highest supported C++
    # level" feature, but in the future it may provide more features.
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.7",
)
