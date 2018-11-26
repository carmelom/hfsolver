from os.path import join as pjoin

# Format expected by setup.py and doc/source/conf.py: string of form "X.Y.Z"
_version_major = 0
_version_minor = 9
_version_micro = ''  # use '' for first of series, number for 1 and above
#_version_extra = 'dev'
_version_extra = ''  # Uncomment this for full releases

# Construct full version string from these.
_ver = [_version_major, _version_minor]
if _version_micro:
    _ver.append(_version_micro)
if _version_extra:
    _ver.append(_version_extra)

__version__ = '.'.join(map(str, _ver))

CLASSIFIERS = ['Development Status :: 3 - Alpha',
               'Environment :: Console',
               'Intended Audience :: Science/Research',
               'License :: OSI Approved :: GNU GPL v3',
               'Operating System :: OS Independent',
               'Programming Language :: Python',
               'Topic :: Scientific/Engineering']

# Description should be a one-liner:
description = 'hfsolver: a numerical solver for the Hartree-Fock equation of state of a homogeneous Bose gas'
# Long description will go up on the pypi page
long_description = '''
hfsolver
========
``hfsolver`` is a numerical solver for the Hartree-Fock equation of state of a homogeneous Bose gas

License
=======
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Copyright (c) 2018 Carmelo Mordini
'''

NAME = 'hfsolver'
MAINTAINER = 'Carmelo Mordini'
MAINTAINER_EMAIL = 'carmelo.mordini@unitn.com'
DESCRIPTION = description
LONG_DESCRIPTION = long_description
URL = 'https://github.com/carmelom/hfsolver'
DOWNLOAD_URL = ''
LICENSE = 'GNU GPL v3'
AUTHOR = 'Carmelo Mordini'
AUTHOR_EMAIL = 'carmelo.mordini@unitn.com'
PLATFORMS = 'OS Independent'
MAJOR = _version_major
MINOR = _version_minor
MICRO = _version_micro
VERSION = __version__
PACKAGES = ['hfsolver']
PACKAGE_DATA = {'hfsolver': [pjoin('data', '*')]}
REQUIRES = []
