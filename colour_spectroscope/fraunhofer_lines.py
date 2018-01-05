#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Fraunhofer Lines
================

Defines objects for the homemade spectroscope spectrum images of Fraunhofer
lines analysis.

References
==========
.. [1]  http://en.wikipedia.org/wiki/Fraunhofer_lines
"""

from __future__ import division, unicode_literals

import bisect
import matplotlib.pyplot
import numpy as np
import os
import pylab
import re

from colour import sRGB_COLOURSPACE, read_image
from colour.plotting import boundaries, decorate, display
from analysis import (
    RGB_spectrum,
    luminance_spd)

__author__ = 'Colour Developers'
__copyright__ = 'Copyright (C) 2013-2018 - Colour Developers'
__license__ = 'New BSD License - http://opensource.org/licenses/BSD-3-Clause'
__maintainer__ = 'Colour Developers'
__email__ = 'colour-science@googlegroups.com'
__status__ = 'Production'

__all__ = ['RESOURCES_DIRECTORY',
           'SUN_SPECTRUM_IMAGE',
           'FRAUNHOFER_LINES_PUBLISHED',
           'FRAUNHOFER_LINES_ELEMENTS_MAPPING',
           'FRAUNHOFER_LINES_NOTABLE',
           'FRAUNHOFER_LINES_CLUSTERED',
           'FRAUNHOFER_LINES_MEASURED',
           'fraunhofer_lines_plot']

RESOURCES_DIRECTORY = os.path.join(os.path.dirname(__file__), 'resources')
SUN_SPECTRUM_IMAGE = str(os.path.join(RESOURCES_DIRECTORY,
                                  'Fraunhofer_Lines_001.png'))

FRAUNHOFER_LINES_PUBLISHED = {
    'y': 898.765,
    'Z': 822.696,
    'A': 759.370,
    'B': 686.719,
    'C': 656.281,
    'a': 627.661,
    'D1': 589.592,
    'D2': 588.995,
    'D3': 587.5618,
    'e': 546.073,
    'E2': 527.039,
    'b1': 518.362,
    'b2': 517.270,
    'b3': 516.891,
    'b4': 516.733,
    'c': 495.761,
    'F': 486.134,
    'd': 466.814,
    'e': 438.355,
    'G': 430.790,
    'h': 410.175,
    'H': 396.847,
    'K': 393.368,
    'L': 382.044,
    'N': 358.121,
    'P': 336.112,
    'T': 302.108,
    't': 299.444}

FRAUNHOFER_LINES_ELEMENTS_MAPPING = {
    'y': 'O2',
    'Z': 'O2',
    'A': 'O2',
    'B': 'O2',
    'C': 'H Alpha',
    'a': 'O2',
    'D1': 'Na',
    'D2': 'Na',
    'D3': 'He',
    'e': 'Hg',
    'E2': 'Fe',
    'b1': 'Mg',
    'b2': 'Mg',
    'b3': 'Fe',
    'b4': 'Mg',
    'c': 'Fe',
    'F': 'H Beta',
    'd': 'Fe',
    'e': 'Fe',
    'G"': 'H Gamma',
    'G': 'Fe',
    'G': 'Ca',
    'h': 'H Delta',
    'H': 'Ca+',
    'K': 'Ca+',
    'L': 'Fe',
    'N': 'Fe',
    'P': 'Ti+',
    'T': 'Fe',
    't': 'Ni'}

FRAUNHOFER_LINES_NOTABLE = (
    'A',
    'B',
    'C',
    'D1',
    'D2',
    'D3',
    'E2',
    'F',
    'G',
    'H',
    'K')

FRAUNHOFER_LINES_CLUSTERED = {'b[1-4]': ('b2', ('b4', 'b3', 'b1',), 'b\n4-1'),
                              'D[1-3]': ('D3', ('D2', 'D1'), 'D\n3-1')}

FRAUNHOFER_LINES_MEASURED = {
    'G': 134,
    'F': 371,
    'b4': 502,
    'E2': 545,
    'D1': 810,
    'a': 974,
    'C': 1095}


def fraunhofer_lines_plot(image=SUN_SPECTRUM_IMAGE):
    """
    Plots the Fraunhofer lines of given image.

    Parameters
    ----------
    image : unicode
        Path to read the image from.

    Returns
    -------
    bool
        Definition success.
    """

    spectrum = RGB_spectrum(read_image(image),
                            FRAUNHOFER_LINES_PUBLISHED,
                            FRAUNHOFER_LINES_MEASURED)

    height = len(spectrum.R.values) / 8
    spd = luminance_spd(spectrum).normalise(height)

    pylab.title('The Solar Spectrum - Fraunhofer Lines')

    wavelengths = spectrum.wavelengths
    input, output = min(wavelengths), max(wavelengths)
    pylab.imshow(sRGB_COLOURSPACE.encoding_cctf(
        np.dstack([spectrum.R.values,
                   spectrum.G.values,
                   spectrum.B.values])),
        extent=[input, output, 0, height])

    pylab.plot(spd.wavelengths,
               spd.values, color='black',
               linewidth=1)

    fraunhofer_wavelengths = np.array(
        sorted(FRAUNHOFER_LINES_PUBLISHED.values()))
    fraunhofer_wavelengths = fraunhofer_wavelengths[
        np.where(np.logical_and(fraunhofer_wavelengths >= input,
                                fraunhofer_wavelengths <= output))]
    fraunhofer_lines_labels = [
        tuple(FRAUNHOFER_LINES_PUBLISHED.keys())[
            tuple(FRAUNHOFER_LINES_PUBLISHED.values()).index(i)]
        for i in fraunhofer_wavelengths]

    y0, y1 = 0, height * .5
    for i, label in enumerate(fraunhofer_lines_labels):

        # Trick to cluster siblings fraunhofer lines.
        from_siblings = False
        for pattern, (first, siblings,
                      specific_label) in FRAUNHOFER_LINES_CLUSTERED.items():
            if re.match(pattern, label):
                if label in siblings:
                    from_siblings = True

                label = specific_label
                break

        power = bisect.bisect_left(wavelengths, fraunhofer_wavelengths[i])
        scale = (spd.get(wavelengths[power]) / height)

        is_large_line = label in FRAUNHOFER_LINES_NOTABLE

        pylab.vlines(fraunhofer_wavelengths[i], y0, y1 * scale,
                     linewidth=1 if is_large_line else 1)

        pylab.vlines(fraunhofer_wavelengths[i], y0, height,
                     linewidth=1 if is_large_line else 1, alpha=0.075)

        if not from_siblings:
            pylab.text(fraunhofer_wavelengths[i],
                       y1 * scale + (y1 * 0.025),
                       label,
                       clip_on=True,
                       ha='center',
                       va='bottom',
                       fontdict={'size': 'large' if is_large_line else'small'})

    r = lambda x: int(x / 100) * 100
    matplotlib.pyplot.xticks(np.arange(r(input), r(output * 1.5), 20))

    settings = {'x_tighten': True,
                'y_tighten': True,
                'x_label': u'Wavelength λ (nm)',
                'y_label': False,
                'legend': False,
                'limits': [input, output, 0, height],
                'y_ticker': True,
                'grid': True}

    boundaries(**settings)

    decorate(**settings)

    return display(**settings)


if __name__ == '__main__':
    fraunhofer_lines_plot()
