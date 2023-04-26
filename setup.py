#!/usr/bin/env python

from setuptools import setup
import os

version_file = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    'src',
    'qmflows',
    '_version.py',
))
version: "dict[str, str]" = {}
with open(version_file, 'r', encoding='utf8') as f:
    exec(f.read(), version)


def readme() -> str:
    """Load the readme file."""
    with open('README.rst', 'r', encoding='utf8') as f:
        return f.read()


docs_require = [
    'sphinx>=2.1,!=3.1.1',
    'sphinx-autodoc-typehints',
    'sphinx_rtd_theme',
    'nbsphinx',
    'jupyter',
    'pandoc',
]

tests_no_optional_require = [
    'assertionlib>=3.1.0',
    'pytest>=6.0',
    'pytest-cov',
    'pytest-mock',
    'typing_extensions'
]

tests_require = tests_no_optional_require.copy()
tests_require += docs_require
tests_require.append("rdkit>=2018.03.1")

setup(
    name='qmflows',
    version=version['__version__'],
    description='Automation of computations in quantum chemistry',
    license='LGPLv3',
    url='https://github.com/SCM-NV/qmflows',
    author='Felipe Zapata',
    author_email='f.zapata@esciencecenter.nl',
    keywords='chemistry workflows simulation materials',
    long_description=readme(),
    package_dir={'': 'src'},
    packages=["qmflows",
              "qmflows.components",
              "qmflows.data",
              "qmflows.data.dictionaries",
              "qmflows.examples",
              "qmflows.examples.Conditional_workflows",
              "qmflows.examples.Constrained_and_TS_optimizations",
              "qmflows.packages",
              "qmflows.parsers",
              "qmflows.templates"],
    package_data={
        "qmflows": ['data/dictionaries/*yaml',
                    'py.typed']
    },
    python_requires='>=3.8',
    classifiers=[
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Typing :: Typed',
    ],
    install_requires=[
        'more-itertools',
        'h5py',
        'numpy',
        'pandas',
        'noodles>=0.3.3',
        'plams>=1.5.1',
        'pyparsing!=3.0.0,<3.1.0',
        'pyyaml>=5.1',
        'filelock',
        'packaging>=1.16.8',
    ],
    extras_require={
        'test': tests_require,
        'test_no_optional': tests_no_optional_require,
        'doc': docs_require,
    }
)
