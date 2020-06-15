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
    'sphinx>=2.1',
    'sphinx-autodoc-typehints',
    'sphinx_rtd_theme',
    'nbsphinx',
    'jupyter'
]

tests_require =  [
    'assertionlib>=2.2.0',
    'mypy',
    'pytest>=5.4',
    'pytest-cov',
    'pytest-mock',
    'pytest-pycodestyle',
    'pytest-pydocstyle>=2.1',
    'typing_extensions'
]
tests_require += docs_require

setup(
    name='qmflows',
    version=version['__version__'],
    description='Automation of computations in quantum chemistry',
    license='Apache 2.0',
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
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.7',
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry'
    ],
    install_requires=['more-itertools', 'h5py', 'numpy', 'pandas', 'noodles==0.3.3',
                      'plams@git+https://github.com/SCM-NV/PLAMS@a5696ce62c09153a9fa67b2b03a750913e1d0924',
                      'pyparsing', 'pyyaml>=5.1', 'filelock'],
    extras_require={
        'test': tests_require,
        'doc': tests_require
    }
)
