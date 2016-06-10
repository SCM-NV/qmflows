#!/usr/bin/env python

from setuptools import setup, find_packages


setup(
    name='QuantumWorkflows',
    version='0.1.2',
    description='Automation of computations in quantum chemistry',
    license='',
    url='',
    author_email='',
    keywords='chemistry workflows simulation materials',
    packages=["qmworks", "qmworks.components",
              "qmworks.data",
              "qmworks.hdf5",
              "qmworks.packages", "qmworks.parsers",
              "qmworks.templates", "qmworks.sshConfig"],
    package_data={
        "qmworks": ['data/templates/*json', 'data/dictionaries/*json']
    },
    classifiers=[
        'Intended Audience :: Science/Research',
        'programming language :: python :: 3.5',
        'development status :: 3 - alpha',
        'intended audience :: science/research',
        'topic :: scientific/engineering :: chemistry'
    ],
    install_requires=['cython', 'h5py', 'numpy', 'plams',
                      'pymonad', 'pyparsing', 'six'],
    dependency_links = ['https://github.com/SCM-NV/qmworks/tarball/master#egg=package-0.1.2']
)
