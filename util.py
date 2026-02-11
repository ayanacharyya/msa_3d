##!/usr/bin/env python3

"""

    Filename :   util.py
    Notes :      Contains various generic utility functions and classes used by the other scripts in MSA-3D, including a function to parse args
    Author :    Ayan
    Created: 06-02-26
"""

from header import *

# -------------------------------------------------------------------------------------------------------
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

# --------------------------------------------------------------------------------------------------------------------
def parse_args():
    '''
    Parse command line arguments. Returns args object
    '''

    parser = argparse.ArgumentParser(description='Produces emission line maps for JWST-PASSAGE data.')

    # ---- common args used widely over the full codebase ------------
    parser.add_argument('--input_dir', metavar='input_dir', type=str, action='store', default=None, help='Where do the input files reside?')
    parser.add_argument('--output_dir', metavar='output_dir', type=str, action='store', default=None, help='Where do you want to store the outputs?')
    parser.add_argument('--system', metavar='system', type=str, action='store', default='hd', help='Which file system is the code being run on?')
    parser.add_argument('--code_dir', metavar='code_dir', type=str, action='store', default='/Users/acharyya/Work/astro/ayan_codes/passage/', help='Where is the source code?')
    parser.add_argument('--mappings_dir', metavar='mappings_dir', type=str, action='store', default='/Users/acharyya/Work/astro/Mappings', help='Where do MAPPINGS files reside?')
    parser.add_argument('--clobber', dest='clobber', action='store_true', default=False, help='Over-write existing plots? Default is no.')
    parser.add_argument('--silent', dest='silent', action='store_true', default=False, help='Suppress some generic print statements? Default is no.')
    parser.add_argument('--keep', dest='keep', action='store_true', default=False, help='Keep existing plot windows open? Default is no.')
    parser.add_argument('--forpaper', dest='forpaper', action='store_true', default=False, help='Format plots to paper quality? Default is no.')
    parser.add_argument('--fortalk', dest='fortalk', action='store_true', default=False, help='Format plots suitable for putting in talks? Default is no.')
    parser.add_argument('--fontsize', metavar='fontsize', type=int, action='store', default=15, help='fontsize of plot labels, etc.; default is 15')

    parser.add_argument('--id', metavar='id', type=str, action='store', default=None, help='Object ID. Default is None')
    parser.add_argument('--do_all_obj', dest='do_all_obj', action='store_true', default=False, help='Reduce spectra and make beam files for ALL detected objects? Default is no.')
    parser.add_argument('--plot_circle_at_arcsec', metavar='plot_circle_at_arcsec', type=float, action='store', default=None, help='Radius in arcseconds of a circle to be plotted on every 2D map; default is None')
    parser.add_argument('--no_text_on_plot', dest='no_text_on_plot', action='store_true', default=False, help='Skip putting text annotations on plot2D? Default is no.')

    # ------- args added for make_msa3d_line_maps.py ------------------------------
    parser.add_argument('--line_list', metavar='line_list', type=str, action='store', default='all', help='Which emission lines to look for? Default is all') # OR set default to 'Lya,OII,Hb,OIII,Ha,Ha+NII,SII,SIII,PaB,He-1083,PaA'
    parser.add_argument('--snr_cut', metavar='snr_cut', type=float, action='store', default=None, help='Impose an SNR cut on the emission line maps to; default is 0')
    parser.add_argument('--flam_max', metavar='flam_max', type=float, action='store', default=10, help='Maximum y-axis limit for f_lambda (in units of 1e-19 ergs/s/cm^2/A); default is None')
    parser.add_argument('--mask_window', metavar='mask_window', type=float, action='store', default=30, help='Wavelength window around expected emission lines to mask out, before fitting continuum, in Angstrom; default is 30')
    parser.add_argument('--group_gap', metavar='group_gap', type=float, action='store', default=150, help='Wavelength window to consider for make friends-of-friends neighbouring line list, in Angstrom; default is 150')
    parser.add_argument('--fit_padding', metavar='fit_padding', type=float, action='store', default=20, help='Wavelength window to pad on either side of a neighbouring line list, in Angstrom; default is 10')
    parser.add_argument('--tie_vdisp', dest='tie_vdisp', action='store_true', default=False, help='Tie the velocity dispersion of all lines to be the same? Default is no.')
    parser.add_argument('--n_cores', metavar='n_cores', type=int, action='store', default=None, help='Number of cores to use in parallel for line fitting; default is None, i.e. all available cores')

    parser.add_argument('--debug_linefit', dest='debug_linefit', type=str, action='store', default=None, help='Debug the line fitting at a specific pixel (comma separated)? Default is no.')
    parser.add_argument('--save_linefit_plot', dest='save_linefit_plot', action='store_true', default=False, help='Save the plot for emission line fitting in each spaxel? Default is no.')
    parser.add_argument('--show_log_flux', dest='show_log_flux', action='store_true', default=False, help='Display line fitting spectrum plots in log-scale in y-axis? Default is no.')

    parser.add_argument('--plot_flux_maps', dest='plot_flux_maps', action='store_true', default=False, help='Plot the line flux maps for a given galaxy? Default is no.')
    parser.add_argument('--plot_line_maps', dest='plot_line_maps', action='store_true', default=False, help='Plot the line flux and kinematic maps for a given galaxy and given line? Default is no.')
    parser.add_argument('--plot_ratio_maps', dest='plot_ratio_maps', action='store_true', default=False, help='Plot the line ratio maps for a given galaxy? Default is no.')
    parser.add_argument('--plot_snr', dest='plot_snr', action='store_true', default=False, help='Plot the SNR map for a given galaxy? Default is no.')

    # ------- args added for make_metallicity_sfr_maps.py ------------------------------
    parser.add_argument('--upto_kpc', metavar='upto_kpc', type=float, action='store', default=5, help='Radius in kpc within which integrated quantitites would be measured, and radial fits would be performed; default is 5')
    parser.add_argument('--plot_ionisation_parameter', dest='plot_ionisation_parameter', action='store_true', default=False, help='Plot the plot_ionisation_parameter map along with metallicity? Default is no.')
    parser.add_argument('--plot_DIG', dest='plot_DIG', action='store_true', default=False, help='Plot DIG diagnostics? Default is no.')
    parser.add_argument('--plot_radial_profiles', dest='plot_radial_profiles', action='store_true', default=False, help='Plot radial profiles corresponding to the 2D maps? Default is no.')

    parser.add_argument('--use_original_NB_grid', dest='use_original_NB_grid', action='store_true', default=False, help='Use the original, unmodified NebulaBayes grid? Default is no.')
    parser.add_argument('--plot_metallicity', dest='plot_metallicity', action='store_true', default=False, help='Plot the metallicity map? Default is no.')
    parser.add_argument('--plot_met_sfr', dest='plot_met_sfr', action='store_true', default=False, help='Plot the metallicity and SFR maps? Default is no.')
    parser.add_argument('--Zdiag', metavar='Zdiag', type=str, action='store', default='KD02_R23', help='Which metallicity diagnostic to use (choose between KD02_R23,R23,R3,O3S2,O3O2,S2,R2,RS32,Te,P25,NB? Default is KD02_R23')
    parser.add_argument('--Zbranch', metavar='Zbranch', type=str, action='store', default='low', help='Which R23 branch to be used (choose between high/low)? Default is low')
    parser.add_argument('--use_C25', dest='use_C25', action='store_true', default=False, help='Use the Cataldi+2025 calibrations rather than Curti+2020? Default is no.')
    parser.add_argument('--debug_Zdiag', dest='debug_Zdiag', action='store_true', default=False, help='Make additional plots to debug the metallicity diagnostic implementation? Default is no.')
    parser.add_argument('--ignore_combined_method', dest='ignore_combined_method', action='store_true', default=False, help='Ignore the combined method (S6 of KD02) while computing R23 metallicity and rely solely on R23? Default is no.')
    parser.add_argument('--exclude_lines', metavar='exclude_lines', type=str, action='store', default='', help='Which lines to be excluded for metallicity measurement with NB? Default is empty string, i.e., use all available lines')
    parser.add_argument('--dered_in_NB', dest='dered_in_NB', action='store_true', default=False, help='Make NebulaBayes de-redden the lines? Default is no (i.e., do de-reddening separately before calling NB)')

    parser.add_argument('--plot_BPT', dest='plot_BPT', action='store_true', default=False, help='Plot BPT? Default is no.')
    parser.add_argument('--plot_AGN_frac', dest='plot_AGN_frac', action='store_true', default=False, help='Plot AGN fraction 2D map (based on BPT diagram)? Default is no.')
    parser.add_argument('--AGN_diag', metavar='AGN_diag', type=str, action='store', default='None', help='Which AGN-SF BPT-like diagnostic to use (choose between VO87,H21,O2O3,O2Hb,Ne3O2? Default is None')
    parser.add_argument('--mask_agn', dest='mask_agn', action='store_true', default=False, help='Mask out the AGN-dominated pixels from all metallicity estimates? Default is no.')

    parser.add_argument('--debug_Zsfr', dest='debug_Zsfr', action='store_true', default=False, help='Debug the metallicity-sfr plots? Default is no.')
    parser.add_argument('--fit_correlation', dest='fit_correlation', action='store_true', default=False, help='Fit a slope between x and y? Default is no.')

    # ------- wrap up and processing args ------------------------------
    args = parser.parse_args()
    if args.line_list != 'all': args.line_list = [item for item in args.line_list.split(',')]

    if args.id is not None: args.id = [int(item) for item in args.id.split(',')]
    if args.debug_linefit is not None: args.debug_linefit = [int(item) for item in args.debug_linefit.split(',')]

    if args.system == 'hd' and not os.path.exists('/Volumes/Elements/'): args.system = 'local'
    if args.line_list == 'all': args.line_list = list(rest_wave_dict.keys())

    survey_name = 'msa_3d'
    root_dir = f'/Users/acharyya/Work/astro/{survey_name}' if 'local' in args.system else f'/Volumes/Elements/acharyya_backup/Work/astro/{survey_name}' if 'hd' in args.system else f'/Users/acharyya/Library/CloudStorage/GoogleDrive-ayan.acharyya@inaf.it/My Drive/{survey_name}' if 'gdrive' in args.system else ''
    args.root_dir = Path(root_dir)

    if args.input_dir is None:
        args.input_dir = args.root_dir / f'{survey_name}_data/'
    if args.output_dir is None:
        args.output_dir = args.root_dir / f'{survey_name}_output/'

    (args.output_dir / 'catalogs').mkdir(exist_ok=True, parents=True)
    (args.output_dir / 'plots').mkdir(exist_ok=True, parents=True)

    args.code_dir = Path(args.code_dir)
    args.mappings_dir = Path(args.mappings_dir)

    args.exclude_lines = args.exclude_lines.split(',') if len(args.exclude_lines) > 0 else []

    if args.fortalk:
        print(f'Setting up plots for talks..')
        setup_plots_for_talks()

    args.Zdiag_arr = args.Zdiag.split(',')
    if args.n_cores is None:
        if args.debug_linefit is not None: args.n_cores = 1
        else: args.n_cores = cpu_count() - 1

    args.extent = (-msa_arcsec_span_x/2, msa_arcsec_span_x/2, -msa_arcsec_span_y/2, msa_arcsec_span_y/2)

    return args

# ------------------------------------------------------------------------------------
def get_zrange_for_line(line, obs_wave_range=[800, 2200]):
    '''
    Computes the redshift range in which a certain emission line is available, for a given observed wavelength window (by default corresponds to NIRISS WFSS)
    Input wavelengths must be in nm
    Returns min and max redshift
    '''
    if type(line) in [float, int, np.float64, np.int64]: rest_wave = line
    elif type(line) in [str, np.str_]: rest_wave = rest_wave_dict[line] / 10 # to convert from Angstrom to nm

    z_min = (obs_wave_range[0] / rest_wave) - 1
    z_max = (obs_wave_range[1] / rest_wave) - 1

    return max(0, z_min), z_max

# ------------------------------------------------------------------------------------
def print_zrange_for_lines(lines, obs_wave_range=[800, 2200]):
    '''
    Prints the redshift range in which a certain emission line is available, for a given observed wavelength window (by default corresponds to NIRISS WFSS)
    Input wavelengths must be in nm
    '''
    lines = np.atleast_1d(lines)

    for line in lines:
        z_min, z_max = get_zrange_for_line(line, obs_wave_range=obs_wave_range)
        if type(line) is float or type(line) is int: line = '%.2f nm' %line

        print(f'Line: {line}, z=[{z_min:.2f}, {z_max:.2f}]')

# ---------------------------------------------------------------------------
def get_kpc_from_arc_at_redshift(arcseconds, redshift):
    '''
    Function to convert arcseconds on sky to physical kpc, at a given redshift
    '''
    d_A = cosmo.angular_diameter_distance(z=redshift)
    kpc = (d_A * arcseconds * u.arcsec).to(u.kpc, u.dimensionless_angles()).value # in kpc
    print('%.2f arcseconds corresponds to %.2F kpc at target redshift of %.2f' %(arcseconds, kpc, redshift))
    return kpc

# ----------------------------------------------------------------------------------------------------------
def get_crossmatch(df1, df2, sep_threshold=0.1, df1_idcol='id', df2_idcol='id'):
    '''
    Determines crossmatch between two dataframes df1 and df2
    df1 and df2 should have df1_idcol, and df2_idcol respectively, and they each have columns: ra, dec
    sep_threshold is in arcseconds
    Returns cross matched dataframe with IDs from both dataframes
    '''
    df1_coords = SkyCoord(df1['ra'], df1['dec'], unit='deg')
    df2_coords = SkyCoord(df2['ra'], df2['dec'], unit='deg')
    nearest_id_in_df2, sep_from_nearest_id_in_df2, _ = df1_coords.match_to_catalog_sky(df2_coords)

    df_crossmatch = pd.DataFrame({'df1_id': df1[df1_idcol].values, 'df2_id': df2[df2_idcol].iloc[nearest_id_in_df2].values, 'sep': sep_from_nearest_id_in_df2.arcsec})
    df_crossmatch = df_crossmatch[df_crossmatch['sep'] < sep_threshold]  # separation within XX arcsecond
    df_crossmatch = df_crossmatch.sort_values('sep').drop_duplicates(subset='df2_id', keep='first').reset_index(drop=True)  # to avoid multiple df1 objects being linked to the same df2 object

    return df_crossmatch

# ------------------------------------------------------------------------
def setup_plots_for_talks():
    '''
    Function to setup plto themes etc for talks
    '''
    plt.style.use('cyberpunk')
    background_for_talks = 'cyberpunk'  # 'dark_background' #'Solarize_Light2' #
    plt.style.use(background_for_talks)
    new_foreground_color = '#FFF1D0'
    #plt.rcParams['grid.color'] = new_foreground_color
    plt.rcParams['text.color'] = new_foreground_color
    plt.rcParams['xtick.color'] = new_foreground_color
    plt.rcParams['ytick.color'] = new_foreground_color
    plt.rcParams['xtick.color'] = new_foreground_color
    plt.rcParams['axes.titlecolor'] = new_foreground_color
    plt.rcParams['axes.labelcolor'] = new_foreground_color
    plt.rcParams['axes.edgecolor'] = new_foreground_color
    plt.rcParams['figure.edgecolor'] = new_foreground_color
    plt.rcParams['savefig.edgecolor'] = new_foreground_color
    plt.rcParams['axes.linewidth'] = 2

    new_background_color = '#120000'
    plt.rcParams['axes.facecolor'] = new_background_color
    plt.rcParams['figure.facecolor'] = new_background_color
    plt.rcParams['savefig.facecolor'] = new_background_color
    plt.rcParams['grid.alpha'] = 0.5
    plt.rcParams['grid.linewidth'] = 0.3

# ------------------------------------------------------------------------------------------------------
def get_sky_region_from_fits_header(header, CDELT1='CD1_1', CDELT2='CD2_2', ORIENT='ORIENTAT'):
    '''
    Function to make an astro RectanguleSkyRegion from a given fits header
    Returns SkyRegion
    '''
    center_ra = header['CRVAL1'] + (header['NAXIS1']/2 - header['CRPIX1']) * header[CDELT1]
    center_dec = header['CRVAL2'] + (header['NAXIS2']/2 - header['CRPIX2']) * header[CDELT2]
    width = np.abs(header[CDELT1]) * header['NAXIS1']
    height = np.abs(header[CDELT2]) * header['NAXIS2']
    angle = header[ORIENT] if ORIENT in header else 0.

    sky_center = SkyCoord(center_ra, center_dec, unit='deg')
    sky_region = RectangleSkyRegion(center=sky_center, width=width * u.deg, height=height * u.deg, angle=angle * u.deg)

    return sky_region

# ------------------------------------------------------------------------------------------------------
def is_point_in_region(sky_coord, data, CDELT1='CD1_1', CDELT2='CD2_2', ORIENT='ORIENTAT'):
    '''
    Function to check if an input sky coordinate lies within the footprint of a given fits file
    Returns True/False
    '''
    if type(data) == str: # in case the input is the filename
        data = fits.open(data)
    header = data[0].header
    sky_region = get_sky_region_from_fits_header(header, CDELT1=CDELT1, CDELT2=CDELT2, ORIENT=ORIENT)

    contains = sky_region.contains(sky_coord, pywcs.WCS(header))

    return contains

# --------------------------------------------------------------------------------------------------------------------
def distance(x, y, x0, y0):
    """
    Return distance between point
    P[x0,y0] and a curve (x,y)
    """
    d_x = x - x0
    d_y = y - y0
    dis = np.sqrt(d_x ** 2 + d_y ** 2)
    return dis

def min_distance(x, y, P, precision=5):
    """
    Compute minimum/a distance/s between
    a point P[x0,y0] and a curve (x,y)
    rounded at `precision`.

    ARGS:
        x, y      (array)
        P         (tuple)
        precision (int)

    Returns min indexes and distances array.
    """
    # compute distance
    d = distance(x, y, P[0], P[1])
    d = np.round(d, precision)
    # find the minima
    glob_min_idxs = np.argwhere(d == np.min(d)).ravel()
    return glob_min_idxs, d

def get_distance_from_line(xdata, ydata, func, method='K01'):
    '''
    Computes distance of each object in the given xdata and ydata (line ratios) arrays, from a given line func(x)
    Returns the distance as an array
    '''
    print(f'Computing distance from {method} line on the BPT diagram..')
    x = np.linspace(-2, 1, 100)
    y = func(x, method)

    min_dist_arr = []
    for P in zip(xdata.flatten(), ydata.flatten()):
        min_idxs, distances = min_distance(x, y, P)
        if len(min_idxs) > 0: min_dist = distances[min_idxs[0]]
        else: min_dist = np.nan
        min_dist_arr.append(min_dist)

    return np.reshape(min_dist_arr, np.shape(xdata))

# -----------------------------------------------------------------
def rebin(array, dimensions=None, scale=None):
    """ Return the array ``array`` to the new ``dimensions`` conserving flux the flux in the bins
    The sum of the array will remain the same

    >>> ar = numpy.array([
        [0,1,2],
        [1,2,3],
        [2,3,4]
        ])
    >>> rebin(ar, (2,2))
    array([
        [1.5, 4.5]
        [4.5, 7.5]
        ])
    Raises
    ------

    AssertionError
        If the totals of the input and result array don't agree, raise an error because computation may have gone wrong

    Reference
    =========
    +-+-+-+
    |1|2|3|
    +-+-+-+
    |4|5|6|
    +-+-+-+
    |7|8|9|
    +-+-+-+
    """
    if dimensions is not None:
        if isinstance(dimensions, float):
            dimensions = [int(dimensions)] * len(array.shape)
        elif isinstance(dimensions, int):
            dimensions = [dimensions] * len(array.shape)
        elif len(dimensions) != len(array.shape):
            raise RuntimeError('')
    elif scale is not None:
        if isinstance(scale, float) or isinstance(scale, int):
            dimensions = map(int, map(round, map(lambda x: x * scale, array.shape)))
        elif len(scale) != len(array.shape):
            raise RuntimeError('')
    else:
        raise RuntimeError('Incorrect parameters to rebin.\n\trebin(array, dimensions=(x,y))\n\trebin(array, scale=a')
    if np.shape(array) == dimensions: return array  # no rebinning actually needed
    import itertools
    # dY, dX = map(divmod, map(float, array.shape), dimensions)

    result = np.zeros(dimensions)
    for j, i in itertools.product(*map(range, array.shape)):
        (J, dj), (I, di) = divmod(j * dimensions[0], array.shape[0]), divmod(i * dimensions[1], array.shape[1])
        (J1, dj1), (I1, di1) = divmod(j + 1, array.shape[0] / float(dimensions[0])), divmod(i + 1,
                                                                                            array.shape[1] / float(
                                                                                                dimensions[1]))

        # Moving to new bin
        # Is this a discrete bin?
        dx, dy = 0, 0
        if (I1 - I == 0) | ((I1 - I == 1) & (di1 == 0)):
            dx = 1
        else:
            dx = 1 - di1
        if (J1 - J == 0) | ((J1 - J == 1) & (dj1 == 0)):
            dy = 1
        else:
            dy = 1 - dj1
        # Prevent it from allocating outide the array
        I_ = np.min([dimensions[1] - 1, I + 1])
        J_ = np.min([dimensions[0] - 1, J + 1])
        result[J, I] += array[j, i] * dx * dy
        result[J_, I] += array[j, i] * (1 - dy) * dx
        result[J, I_] += array[j, i] * dy * (1 - dx)
        result[J_, I_] += array[j, i] * (1 - dx) * (1 - dy)
    allowError = 0.1
    if array.sum() > 0: assert (array.sum() < result.sum() * (1 + allowError)) & (array.sum() > result.sum() * (1 - allowError))
    return result

# --------------------------------------------------------------------------------------------------------------------
def get_custom_cmap(cmap_name, cmap_path=None):
    '''
    Computes custom colormaps, code borrowed from Eduardo Vitral

    SET COLORMAPS FROM:
    Crameri, F., G.E. Shephard, and P.J. Heron (2020)
    The misuse of colour in science communication,
    Nature Communications, 11, 5444.
    https://doi.org/10.1038/s41467-020-19160-7
    ---
    Crameri, F.: Geodynamic diagnostics, scientific visualisation
    and StagLab 3.0, Geosci. Model Dev. Discuss.,
    https://doi.org/10.5194/gmd-2017-328, 2018.

    Returns the colormap
    '''

    if cmap_path is None: cmap_path = HOME / 'Work/astro/ayan_codes/ScientificColourMaps8'
    cpal_quant = np.loadtxt(cmap_path / cmap_name / f'{cmap_name}.txt')
    mycmap_quant = mplcolors.ListedColormap(cpal_quant, name=cmap_name)

    return mycmap_quant

def get_combined_cmap(breaks, cmaps, new_name='my_colormap'):
    '''
    Combines an arbitrary number of matplotlib colormaps at given break points
    Returns a new colormap
    Adapted from https://stackoverflow.com/questions/31051488/combining-two-matplotlib-colormaps
    '''
    colors = []
    breaks[0] = min(breaks[0], breaks[1])
    breaks[-1] = max(breaks[-1], breaks[-2])

    for index, cmap in enumerate(cmaps):
        ncolors = int(256 * (breaks[index + 1] - breaks[index]) / (breaks[-1] - breaks[0]))
        this_colors = mplcolormaps[cmap](np.linspace((breaks[index] - breaks[0]) / (breaks[-1] - breaks[0]), (breaks[index + 1] - breaks[0]) /(breaks[-1] - breaks[0]), ncolors))
        colors.append(this_colors)
    colors = np.vstack(colors)
    new_cmap = mplcolors.LinearSegmentedColormap.from_list(new_name, colors)

    return new_cmap

# --------------------------------------------------------------------------------------------------------------------
def unzip_and_delete(zip_file, destination_path):
    '''
    Unzips given zip file to given destination path and removes the zip file
    '''
    zip_file = Path(zip_file)
    print(f'Unpacking {zip_file} in to {destination_path}..')
    if zip_file.suffix == '.gz':
        with gzip.open(zip_file, 'rb') as gz_file:
            with open(zip_file.parent / zip_file.stem, 'wb') as out_file:
                shutil.copyfileobj(gz_file, out_file)
    else:
        shutil.unpack_archive(zip_file, destination_path)
    os.remove(zip_file)  # remove zipped files after unzipping

# ------------------------------------------------------------------------------
def insert_line_in_file(line, pos, filename, output=None):
    '''
    Function to inserting a line in a file, given the filename and what and where to insert
    '''
    f = open(filename, 'r')
    contents = f.readlines()
    f.close()

    if pos == -1: pos = len(contents)  # to append to end of file
    contents.insert(pos, line)

    if output is None: output = filename
    f = open(output, 'w')
    contents = ''.join(contents)
    f.write(contents)
    f.close()
    return

# --------------------------------------------------------------------------------------------------
class smart_dict(dict):
    '''
    In order to be able to return the key itself froma dict for missing keys
    From https://stackoverflow.com/questions/6229073/how-to-make-a-dictionary-that-returns-key-for-keys-missing-from-the-dictionary-i
    '''
    def __missing__(self, key):
        return key

# --------------------------------------------------------------------------------------------------------------------
def annotate_axes(ax, xlabel, ylabel, args=None, fontsize=10, fontfactor=1, label='', clabel='', hide_xaxis=False, hide_yaxis=False, hide_cbar=True, p=None, hide_cbar_ticks=False, cticks_integer=True):
    '''
    Annotates the axis of a given 2D image
    Returns the axis handle
    '''
    if args is not None: fontsize, fontfactor = args.fontsize, args.fontfactor
    ax.text(0.05, 0.9, label, c='k', fontsize=fontsize/fontfactor, ha='left', va='top', bbox=dict(facecolor='white', edgecolor='black', alpha=0.9), transform=ax.transAxes)

    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=3, prune='both'))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(4))
    if hide_xaxis:
        ax.tick_params(axis='x', which='major', labelsize=fontsize, labelbottom=False)
    else:
        ax.set_xlabel(xlabel, fontsize=fontsize)
        ax.tick_params(axis='x', which='major', labelsize=fontsize, labelbottom=True)

    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3, prune='both'))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(4))
    if hide_yaxis:
        ax.tick_params(axis='y', which='major', labelsize=fontsize, labelleft=False)
    else:
        ax.set_ylabel(ylabel, fontsize=fontsize)
        ax.tick_params(axis='y', which='major', labelsize=fontsize, labelleft=True)

    if not hide_cbar and p is not None:
        cax = inset_axes(ax, width="5%", height="100%", loc='right', bbox_to_anchor=(0.05, 0, 1, 1), bbox_transform=ax.transAxes, borderpad=0)
        cbar = plt.colorbar(p, cax=cax, orientation='vertical')
        cbar.set_label(clabel, fontsize=fontsize)

        cbar.locator = ticker.MaxNLocator(integer=cticks_integer, nbins=4)#, prune='both')
        cbar.update_ticks()
        if hide_cbar_ticks:
            cbar.ax.set_yticklabels([])
        else:
            cbar.ax.tick_params(labelsize=fontsize)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def save_fig(fig, fig_dir, figname, args, silent=False):
    '''
    Saves a given figure handle as a given output filename
    '''

    if args.fortalk:
        #mplcyberpunk.add_glow_effects()
        #try: mplcyberpunk.make_lines_glow()
        #except: pass
        try: mplcyberpunk.make_scatter_glow()
        except: pass

    fig_dir.mkdir(exist_ok=True, parents=True)
    figname = fig_dir / figname
    fig.savefig(figname, transparent=args.fortalk)
    if not silent: print(f'\nSaved figure as {figname}')
    plt.show(block=False)

    return

# --------------------------------------------------------------------------------------------------------------------
def make_colorbar_top(fig, axes, clabel, cmap, cmin, cmax, ncbins, fontsize, aspect=60):
    '''
    Creates a shared colorbar for the whole figure, at the top of the figure
    Returns figure handle
    '''
    norm = mplcolors.Normalize(vmin=cmin, vmax=cmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

    cbar = fig.colorbar(sm, ax=axes, location='top', shrink=0.95, pad=0.01, aspect=aspect)
    cbar.set_label(clabel, fontsize=fontsize, labelpad=5)    
    cbar.ax.tick_params(labelsize=fontsize)
    cbar.locator = ticker.MaxNLocator(integer=False, nbins=ncbins)#, prune='both')
    cbar.update_ticks()

    return fig

# --------------------------------------------------------------------------------------------------------------------
def take_safe_log_ratio(num_map, den_map, skip_log=False):
    '''
    Takes the log of ratio of two 2D masked arrays by properly accounting for bad values so as to avoid math errors
    Returns 2D masked array
    '''
    if np.ma.isMaskedArray(num_map):
        net_mask = num_map.mask | den_map.mask
        num_map = num_map.data
        den_map = den_map.data
    else:
        net_mask = False

    bad_mask = (unp.nominal_values(num_map) <= 0) | (unp.nominal_values(den_map) <= 0) | (~np.isfinite(unp.nominal_values(num_map))) | (~np.isfinite(unp.nominal_values(den_map))) | (~np.isfinite(unp.std_devs(num_map))) | (~np.isfinite(unp.std_devs(den_map)))
    num_map[bad_mask] = 1e-9  # arbitrary fill value to bypass unumpy's inability to handle math domain errors
    den_map[bad_mask] = 1e-9  # arbitrary fill value to bypass unumpy's inability to handle math domain errors
    ratio_map = num_map / den_map
    if not skip_log: ratio_map = unp.log10(ratio_map)
    ratio_map[bad_mask | net_mask] = -99.
    ratio_map = np.ma.masked_where(bad_mask | net_mask, ratio_map)

    return ratio_map

# --------------------------------------------------------------------------------------------------------------------
def take_safe_log_sum(map1, map2, skip_log=False):
    '''
    Takes the log of the sum of two 2D masked arrays by properly accounting for bad values so as to avoid math errors
    Returns 2D masked array
    '''
    if np.ma.isMaskedArray(map1):
        net_mask = map1.mask | map2.mask
        map1 = map1.data
        map2 = map2.data
    else:
        net_mask = False

    bad_mask = (unp.nominal_values(map1) <= 0) | (unp.nominal_values(map2) <= 0) | (~np.isfinite(unp.nominal_values(map1))) | (~np.isfinite(unp.nominal_values(map2))) | (~np.isfinite(unp.std_devs(map1))) | (~np.isfinite(unp.std_devs(map2)))
    map1[bad_mask] = 1e-9  # arbitrary fill value to bypass unumpy's inability to handle math domain errors
    map2[bad_mask] = 1e-9  # arbitrary fill value to bypass unumpy's inability to handle math domain errors
    sum_map = map1 + map2
    if not skip_log: sum_map = unp.log10(sum_map)
    sum_map[bad_mask | net_mask] = -99.
    sum_map = np.ma.masked_where(bad_mask | net_mask, sum_map)

    return sum_map

# --------------------------------------------------------------------------------------------------------------------
def get_kappa(x, i):
    '''
    To calculate kappa according to Clayton, Cardelli & Mathis 1989 dust law
    From ayan_codes/mage_project/ayan/mage.py
    '''
    Rv = 3.1  # Clayton Cardelli Mathis 1989
    x = np.array(x)
    if i == 1:
        a = 0.574 * x ** 1.61
        b = -0.527 * x ** 1.61
    elif i == 2:
        y = x - 1.82
        a = 1 + 0.17699 * y - 0.50447 * y ** 2 - 0.02427 * y ** 3 + 0.72085 * y ** 4 + 0.01979 * y ** 5 - 0.77530 * y ** 6 + 0.32999 * y ** 7
        b = 1.41338 * y + 2.28305 * y ** 2 + 1.07233 * y ** 3 - 5.38434 * y ** 4 - 0.62251 * y ** 5 + 5.30260 * y ** 6 - 2.09002 * y ** 7
    elif i == 3:
        a = 1.752 - 0.316 * x - 0.104 / ((x - 4.67) ** 2 + 0.341)
        b = -3.090 + 1.825 * x + 1.206 / ((x - 4.62) ** 2 + 0.263)
    elif i == 4:
        a = 1.752 - 0.316 * x - 0.104 / ((x - 4.67) ** 2 + 0.341) - 0.04473 * (x - 5.9) ** 2 - 0.009779 * (x - 5.9) ** 3
        b = -3.090 + 1.825 * x + 1.206 / ((x - 4.62) ** 2 + 0.263) + 0.2130 * (x - 5.9) ** 2 - 0.1207 * (x - 5.9) ** 3
    elif i == 5:
        a = -1.073 - 0.628 * (x - 8) + 0.137 * (x - 8) ** 2 - 0.070 * (x - 8) ** 3
        b = 13.670 + 4.257 * (x - 8) - 0.420 * (x - 8) ** 2 + 0.374 * (x - 8) ** 3
    return a * Rv + b

def get_full_kappa(wave, inAngstrom=True):
    '''
    To calculate kappa for a rang eof wavelengths, according to CCM89 dust law
    From ayan_codes/mage_project/ayan/mage.py
    '''
    flag = 0
    if type(wave) in [float, int, np.float64]:
        wave = [float(wave)]
        flag = 1
    wave = np.array(wave)
    if inAngstrom: wave /= 1e4  # to convert to micron
    x = 1. / wave
    k = np.zeros(len(x))
    k += get_kappa(x, 1) * ((x >= 0.3) & (x <= 1.1))
    k += get_kappa(x, 2) * ((x > 1.1) & (x <= 3.3))
    k += get_kappa(x, 3) * ((x > 3.3) & (x < 5.9))
    k += get_kappa(x, 4) * ((x >= 5.9) & (x <= 8.))
    k += get_kappa(x, 5) * ((x >= 8.) & (x <= 10.))
    if flag: k = k[0]  # if a single wavelength was input as a float, output a float (not array)
    return k

def get_dereddened_flux(obs_flux_map, wavelength, EB_V, inAngstrom=True):
    '''
    Calculates and returns dereddened fluxes
    Does not deal with uncertainties, for now
    From ayan_codes/mage_project/ayan/mage.py
    '''
    if EB_V is None: EB_V = np.zeros(np.shape(obs_flux_map))
    kappa = get_full_kappa(wavelength, inAngstrom=inAngstrom)
    A_map = kappa * EB_V
    flux_corrected_map = obs_flux_map * 10 ** (0.4 * A_map)

    return flux_corrected_map
# --------------------------------------------------------------------------------------------------------------------

