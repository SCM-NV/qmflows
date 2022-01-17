#!/usr/bin/env python

from setuptools import setup
import os

here = os.path.abspath(os.path.dirname(__file__))
version = {}
with open(os.path.join(here, 'src', 'qmflows', '__version__.py')) as f:
    exec(f.read(), version)


def readme():
    with open('README.rst') as f:
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
    'pytest>=5.4',
    'pytest-cov',
    'pytest-mock',
    'pytest-pycodestyle',
    'pytest-pydocstyle>=2.1',
    'typing_extensions'
]

tests_require =  [
    'mypy',
    'types-PyYAML',
    'types-setuptools',
]
tests_require += tests_no_optional_require
tests_require += docs_require

setup(
    name='qmflows',
    version=version['__version__'],
    description='Automation of computations in quantum chemistry',
    license='LGPLv3',
    url='https://github.com/SCM-NV/qmflows',
    author=['Felipe Zapata'],
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
    python_requires='>=3.6',
    classifiers=[
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
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
        'noodles==0.3.3',
        'plams>=1.5.1',
        'pyparsing<3.0',
        'pyyaml>=5.1',
        'filelock',
    ],
    extras_require={
        'test': tests_require,
        'test_no_optional': tests_no_optional_require,
        'doc': docs_require,
    }
)
