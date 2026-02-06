#!/usr/bin/env python3

"""

    Filename :   header.py
    Notes :      Header file for importing packages/modules and declaring global variables required for working with MSA-3D data.
    Author :    Ayan
    Created: 06-02-26
"""

import numpy as np
import argparse
import os
import copy
import importlib
import glob
import sys
import re
import json
import shutil
import itertools
import io

from datetime import datetime, timedelta
from pathlib import Path
import pandas as pd
from importlib import reload
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.odr import ODR, Model, RealData


import requests
from urllib.parse import quote as urlencode
from urllib.request import urlretrieve
import gzip

from astropy.table import Table, join
from astropy import units as u
from astropy.coordinates import SkyCoord, angular_separation
from astropy.nddata.utils import Cutout2D
from astropy import wcs as pywcs
from astropy.io import fits
from astropy.cosmology import Planck18
from astropy.cosmology import FlatLambdaCDM

from uncertainties import unumpy as unp
from uncertainties import ufloat

import warnings
warnings.filterwarnings("ignore")

import pprint
pp = pprint.PrettyPrinter(indent=4)

from os import PathLike

from NebulaBayes import NB_Model
import NebulaBayes

import matplotlib
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib import colors as mplcolors
from matplotlib import colormaps as mplcolormaps
from matplotlib import cm as mpl_cm
from matplotlib.backends.backend_pdf import PdfPages
import mplcyberpunk
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection

from matplotlib import pyplot as plt
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['ytick.right'] = True
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['xtick.top'] = True

HOME = Path.home()

rest_wave_dict = {'Lya': 121.6, 'OII': 372.7, 'NeIII-3867':386.9, 'Hd': 434.0, 'OIII-4363': 436.3, 'Hb': 486.1, 'OIII': 500.7,
                  'OI-6302':630., 'Ha+NII': 655.5, 'Ha': 656.3, 'SII': 671.7, 'ArIII-7138':713.6,
                  'SIII': 953.0, 'PaD': 1004.6, 'PaG': 1093.5, 'PaB': 1281.4,
                  'PaA': 1874.5, 'OII-7325':732.0, 'HeI-5877':587.6, 'Hg':434.0, 'Hd':410.2}  # approximate wavelengths in nm

# Set up necessary variables for cosmological calculations.
cosmo = FlatLambdaCDM(H0=69.5, Om0=0.285, Ob0=0.0461)

# --- Enable autoreload when running in IPython ---
try:
    from IPython import get_ipython
    ipython = get_ipython()
    if ipython is not None:
        ipython.run_line_magic('load_ext', 'autoreload')
        ipython.run_line_magic('reload_ext', 'autoreload')
        ipython.run_line_magic('autoreload', '2')
except Exception:
    pass
