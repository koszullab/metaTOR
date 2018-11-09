#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""A pipeline for binning metagenomic datasets from 3C data.
"""

from setuptools import setup, find_packages
import setuptools.command.install as install

CLASSIFIERS = [
    "Development Status :: 3 - Alpha",
    "Intended Audience :: Science/Research",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3.6",
    "Topic :: Scientific/Engineering",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Topic :: Scientific/Engineering :: Visualization",
    "Operating System :: Microsoft :: Windows",
    "Operating System :: POSIX",
    "Operating System :: Unix",
    "Operating System :: MacOS",
]

name = "metator"

MAJOR = 0
MINOR = 1
MAINTENANCE = 0
VERSION = f"{MAJOR}.{MINOR}.{MAINTENANCE}"

LICENSE = "GPLv3"
URL = "https://github.com/koszullab/metator"

with open("requirements.txt", "r") as f:
    REQUIREMENTS = f.read().splitlines()

with open("metator/version.py", "w") as f:
    f.write("__version__ = '{}'\n".format(VERSION))


setup(
    name=name,
    author="lyam.baudry@pasteur.fr",
    description=__doc__,
    version=VERSION,
    license=LICENSE,
    classifiers=CLASSIFIERS,
    url=URL,
    packages=find_packages(),
    package_data={"metator": ("bin/*.sh", "share/*")},
    include_package_data=True,
    install_requires=REQUIREMENTS,
    entry_points={"console_scripts": ["metator=metator.metator:main"]},
)
