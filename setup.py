#!/usr/bin/env python

from setuptools import setup


setup(
    name='qmworks',
    version='0.1.4',
    description='Automation of computations in quantum chemistry',
    license='',
    url='https://github.com/SCM-NV/qmworks',
    author_email='',
    keywords='chemistry workflows simulation materials',
    package_dir={'': 'src'},
    packages=["qmworks", "qmworks.components",
              "qmworks.data",
              "qmworks.data.dictionaries",
              "qmworks.hdf5",
              "qmworks.packages", "qmworks.parsers",
              "qmworks.templates"],
    package_data={
        "qmworks": ['data/templates/*json', 'data/dictionaries/*json']
    },
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.5',
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry'
    ],
    install_requires=['h5py', 'numpy', 'noodles', 'plams',
                      'pymonad', 'pyparsing', 'six', 'tinydb',
                      'noodles[prov, xenon, numpy]', 'pyxenon',
                      'filelock', 'msgpack-python'],
    extras_require={'test': ['nose', 'coverage']},
    dependency_links=[
        "https://github.com/SCM-NV/plams/tarball/master#egg=plams"
    ]

)


