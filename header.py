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
import requests
import gzip

from datetime import datetime, timedelta
from pathlib import Path
import pandas as pd
from importlib import reload
from scipy.interpolate import interp1d, UnivariateSpline
from scipy.ndimage import gaussian_filter1d, uniform_filter1d
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.odr import ODR, Model, RealData

from os import PathLike
from multiprocessing import Pool, cpu_count
from functools import partial

from astropy.table import Table, join
from astropy import units as u
from astropy.coordinates import SkyCoord, angular_separation
from astropy.nddata.utils import Cutout2D
from astropy import wcs as pywcs
from astropy.io import fits
from astropy.cosmology import Planck18
from astropy.cosmology import FlatLambdaCDM
from astropy.constants import c

from uncertainties import unumpy as unp
from uncertainties import ufloat

import warnings
warnings.filterwarnings("ignore")

import pprint
pp = pprint.PrettyPrinter(indent=4)

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
import matplotlib.ticker as ticker
from matplotlib.collections import PatchCollection
from matplotlib import pyplot as plt
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['ytick.right'] = True
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['xtick.top'] = True

HOME = Path.home()

rest_wave_dict = {'OIII-4363': 4363.209,
                  'H-beta': 4861.333,
                  'OIII-5007': 5008.239,
                  'OI-6302': 6302.047,
                  'NII-6548': 6549.86,
                  'H-alpha': 6564.632,
                  'NII-6584': 6585.273,
                  'SII-6717': 6718.30,
                  'SII-6730': 6732.674,
                  }  # approximate wavelengths in Angstroms

# Set up necessary variables for cosmological calculations.
cosmo = FlatLambdaCDM(H0=69.5, Om0=0.285, Ob0=0.0461)
c_km_s = c.to(u.km / u.s).value

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
