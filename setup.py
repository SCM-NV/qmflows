#!/usr/bin/env python

from setuptools import setup


setup(
    name='qmworks',
    version='0.1.2',
    description='Automation of computations in quantum chemistry',
    license='',
    url='https://github.com/SCM-NV/qmworks',
    author_email='',
    keywords='chemistry workflows simulation materials',
    packages=["qmworks", "qmworks.components",
              "qmworks.data",
              "qmworks.hdf5",
              "qmworks.packages", "qmworks.parsers",
              "qmworks.templates"],
    package_data={
        "qmworks": ['data/templates/*json', 'data/dictionaries/*json']
    },
    classifiers=[
        'Intended Audience :: Science/Research',
        'programming language :: python :: 3.5',
        'development status :: 4 - beta',
        'intended audience :: science/research',
        'topic :: scientific/engineering :: chemistry'
    ],
    install_requires=['cython', 'h5py', 'numpy', 'noodles', 'plams',
                      'pymonad', 'pyparsing', 'six', 'tinydb'],
    extras_require={'test': ['nose', 'coverage']}

)


