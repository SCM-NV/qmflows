#!/usr/bin/env python

from setuptools import setup


def readme():
    with open('README.rst') as f:
        return f.read()


setup(
    name='qmflows',
    version='0.2.1',
    description='Automation of computations in quantum chemistry',
    license='Apache 2.0',
    url='https://github.com/SCM-NV/qmflows',
    author='felipe zapata',
    author_email='f.zapata@esciencecenter.nl',
    keywords='chemistry workflows simulation materials',
    long_description=readme(),
    package_dir={'': 'src'},
    packages=["qmflows", "qmflows.components",
              "qmflows.data",
              "qmflows.data.dictionaries",
              "qmflows.examples",
              "qmflows.hdf5",
              "qmflows.packages", "qmflows.parsers",
              "qmflows.templates"],
    package_data={
        "qmflows": ['data/templates/*json', 'data/dictionaries/*json']
    },
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.5',
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry'
    ],
    install_requires=['h5py', 'numpy', 'noodles', 'plams>=1.1.1',
                      'pymonad', 'pyparsing', 'six', 'tinydb',
                      'noodles[prov, xenon, numpy]', 'pyxenon',
                      'filelock', 'msgpack-python',
                      'sphinx_rtd_theme'],
    extras_require={
        'test': ['nose', 'coverage', 'nbsphinx'],
        'doc': ['sphinx', 'sphinx_rtd_theme', 'nbsphinx']
    },
    dependency_links=[
        "https://github.com/SCM-NV/plams/tarball/master#egg=plams"]
)
