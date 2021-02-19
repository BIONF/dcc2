#!/bin/env python

#######################################################################
#  Copyright (C) 2020 Vinh Tran
#
#  This file is part of dcc2.
#
#  dcc2 is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  dcc2 is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with greedyFAS.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

from setuptools import setup, find_packages

with open("README.md", "r") as input:
    long_description = input.read()

setup(
    name="dcc2",
    version="0.1.0",
    python_requires='>=3.7.0',
    description="Dynamic core ortholog compilation tool",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Vinh Tran",
    author_email="tran@bio.uni-frankfurt.de",
    url="https://github.com/BIONF/dcc2",
    packages=find_packages(),
    package_data={'': ['*']},
    install_requires=[
        'biopython',
        'bs4',
        'omadb',
        'pyfaidx',
        'tqdm'
    ],
    entry_points={
        'console_scripts': ["dcc.prepare = dcc2.prepareDcc:main",
                            "dcc.parseOrthoxml = dcc2.parseOrthoxml:main",
                            "dcc.parseOmaById = dcc2.parseOmaById:main",
                            "dcc.parseOmaBySpec = dcc2.parseOmaBySpec:main"],
    },
    license="GPL-3.0",
    classifiers=[
        "Environment :: Console",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: End Users/Desktop",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
    ],
)
