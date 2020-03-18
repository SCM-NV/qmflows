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
        "qmflows": ['data/dictionaries/*yaml']
    },
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.7',
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry'
    ],
    install_requires=['more-itertools', 'h5py', 'numpy', 'pandas', 'noodles==0.3.3',
                      'plams@git+https://github.com/SCM-NV/PLAMS@master',
                      'pyparsing', 'pyyaml>=5.1', 'filelock'],
    extras_require={
        'test': ['pytest', 'pytest-cov', 'pytest-mock', 'nbsphinx', 'assertionlib'],
        'doc': ['sphinx>=2.1', 'sphinx-autodoc-typehints', 'sphinx_rtd_theme', 'nbsphinx']
    }
)
