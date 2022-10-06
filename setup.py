import sys

try:
    from skbuild import setup
except ImportError:
    print(
        "Please update pip, you need pip 10 or greater,\n"
        " or you need to install the PEP 518 requirements in pyproject.toml yourself",
        file=sys.stderr,
    )
    raise

from setuptools import find_packages

setup(
    name="PyMatching",
    version="2.0.0",
    description="A package for decoding quantum error correcting codes using minimum-weight perfect matching.",
    author="Oscar Higgott and Craig Gidney",
    license="Apache",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    cmake_install_dir="src/pymatching",
    include_package_data=True,
    extras_require={"test": ["pytest"]},
    python_requires=">=3.6",
    install_requires=["networkx", "retworkx", "stim"]
)