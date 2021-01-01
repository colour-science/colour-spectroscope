# -*- coding: utf-8 -*-
"""
Colour - Spectroscope
=====================

Analysis of the *Fraunhofer* lines in images captured with the homemade
spectroscope.

Subpackages
-----------
-   fraunhofer : Analysis of the *Fraunhofer* lines.
"""

from __future__ import absolute_import

import numpy as np
import os
import subprocess

import colour

from .fraunhofer import calibrated_RGB_spectrum, fraunhofer_lines_plot

__author__ = 'Colour Developers'
__copyright__ = 'Copyright (C) 2015-2021 - Colour Developers'
__license__ = 'New BSD License - https://opensource.org/licenses/BSD-3-Clause'
__maintainer__ = 'Colour Developers'
__email__ = 'colour-developers@colour-science.org'
__status__ = 'Production'

__all__ = [
    'calibrated_RGB_spectrum',
    'fraunhofer_lines_plot',
]

RESOURCES_DIRECTORY = os.path.join(os.path.dirname(__file__), 'resources')

__application_name__ = 'Colour - Spectroscope'

__major_version__ = '0'
__minor_version__ = '1'
__change_version__ = '0'
__version__ = '.'.join(
    (__major_version__,
     __minor_version__,
     __change_version__))  # yapf: disable

try:
    version = subprocess.check_output(
        ['git', 'describe'], cwd=os.path.dirname(__file__)).strip()
    version = version.decode('utf-8')
except Exception:
    version = __version__

colour.utilities.ANCILLARY_COLOUR_SCIENCE_PACKAGES[
    'colour-spectroscope'] = version

# TODO: Remove legacy printing support when deemed appropriate.
try:
    np.set_printoptions(legacy='1.13')
except TypeError:
    pass
