# Available at setup time due to pyproject.toml
from setuptools import setup, Extension
import pybind11
import glob
import os

__version__ = "0.0.1b0"


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


def ls_dir(base_dir):
    return [
        os.path.join(dirpath, fname)
        for (dirpath, dirnames, files) in os.walk(base_dir)
        for fname in files
    ]


def package_content():
    files = ls_dir("nlwpy/")
    return [
        fpath.replace("nlwpy/", "", 1)
        for fpath in files
        if "__pycache__" not in fpath and fpath.endswith(".py")
    ]


setup(
    name="nlwpy",
    license="BSD-3",
    version=__version__,
    author="Gleb Belov",
    author_email="gleb@ampl.com",
    url="https://github.com/ampl/mp",
    description="Python API for the AMPL NL Writer library",
    long_description="",
    packages=["nlwpy"],
    ext_modules=[
        Extension(
            "_nlwpy",
            ["nlwpy/src/nlw_bindings.cc"] + sorted(glob.glob("./src/" + "*.cc")),
            extra_compile_args=compile_args(),
            include_dirs=["include", pybind11.get_include()],
            define_macros=[("VERSION_INFO", __version__)],
        ),
    ],
    package_data={"": package_content()},
    extras_require={"test": "pytest"},
    zip_safe=False,
    python_requires=">=3.7",
)
