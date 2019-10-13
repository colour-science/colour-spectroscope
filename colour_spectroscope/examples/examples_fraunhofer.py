# -*- coding: utf-8 -*-
"""
Showcases examples of analysis of the *Fraunhofer* lines in images captured
with the homemade spectroscope.
"""

import os

from colour import read_image
from colour.utilities import message_box
from colour.plotting import COLOUR_STYLE_CONSTANTS
from colour_spectroscope import RESOURCES_DIRECTORY, fraunhofer_lines_plot
from colour_spectroscope.fraunhofer.plotting import FRAUNHOFER_LINES_MEASURED

message_box(('The Solar Spectrum - "Fraunhofer" Lines'))

SUN_SPECTRUM_IMAGE = str(
    os.path.join(RESOURCES_DIRECTORY, 'Fraunhofer_Lines_001.png'))

fraunhofer_lines_plot(
    COLOUR_STYLE_CONSTANTS.colour.colourspace.cctf_decoding(
        read_image(SUN_SPECTRUM_IMAGE)), FRAUNHOFER_LINES_MEASURED)
