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
    parser.add_argument('--clobber', dest='clobber', action='store_true', default=False, help='Over-write existing plots? Default is no.')
    parser.add_argument('--silent', dest='silent', action='store_true', default=False, help='Suppress some generic print statements? Default is no.')
    parser.add_argument('--keep', dest='keep', action='store_true', default=False, help='Keep existing plot windows open? Default is no.')
    parser.add_argument('--forpaper', dest='forpaper', action='store_true', default=False, help='Format plots to paper quality? Default is no.')
    parser.add_argument('--fortalk', dest='fortalk', action='store_true', default=False, help='Format plots suitable for putting in talks? Default is no.')
    parser.add_argument('--fontsize', metavar='fontsize', type=int, action='store', default=15, help='fontsize of plot labels, etc.; default is 15')

    parser.add_argument('--id', metavar='id', type=str, action='store', default=None, help='Object ID. Default is None')

    # ------- args added for make_msa3d_line_maps.py ------------------------------
    parser.add_argument('--do_all_obj', dest='do_all_obj', action='store_true', default=False, help='Reduce spectra and make beam files for ALL detected objects? Default is no.')
    parser.add_argument('--zmin', metavar='zmin', type=float, action='store', default=0.5, help='minimum of redshift range within which to search for lines; default is 0.5')
    parser.add_argument('--zmax', metavar='zmax', type=float, action='store', default=1.0, help='maximum of redshift range within which to search for lines; default is None')
    parser.add_argument('--line_list', metavar='line_list', type=str, action='store', default='all', help='Which emission lines to look for? Default is all') # OR set default to 'Lya,OII,Hb,OIII,Ha,Ha+NII,SII,SIII,PaB,He-1083,PaA'
    parser.add_argument('--do_download', dest='do_download', action='store_true', default=False, help='Download fresh cube from Google Drive? Default is no.')
    parser.add_argument('--snr_cut', metavar='snr_cut', type=float, action='store', default=None, help='Impose an SNR cut on the emission line maps to; default is 0')
    parser.add_argument('--flam_max', metavar='flam_max', type=float, action='store', default=None, help='Maximum y-axis limit for f_lambda (in units of 1e-19 ergs/s/cm^2/A); default is None')
    
    parser.add_argument('--debug_offset', dest='debug_offset', action='store_true', default=False, help='Do extra plots and prints for debugging offset calculation from center? Default is no.')

    parser.add_argument('--plot_radial_profiles', dest='plot_radial_profiles', action='store_true', default=False, help='Plot radial profiles corresponding to the 2D maps? Default is no.')
    parser.add_argument('--plot_ratio_maps', dest='plot_ratio_maps', action='store_true', default=False, help='Plot the line ratio maps for a given 2D plot? Default is no.')

    parser.add_argument('--plot_ionisation_parameter', dest='plot_ionisation_parameter', action='store_true', default=False, help='Plot the plot_ionisation_parameter map along with metallicity? Default is no.')
    parser.add_argument('--plot_DIG', dest='plot_DIG', action='store_true', default=False, help='Plot DIG diagnostics? Default is no.')
    parser.add_argument('--plot_BPT', dest='plot_BPT', action='store_true', default=False, help='Plot BPT? Default is no.')

    parser.add_argument('--plot_metallicity', dest='plot_metallicity', action='store_true', default=False, help='Plot the metallicity map instead of the full diagnostic figure? Default is no.')
    parser.add_argument('--Zbranch', metavar='Zbranch', type=str, action='store', default='low', help='Which R23 branch to be used (choose between high/low)? Default is low')
    parser.add_argument('--Zdiag', metavar='Zdiag', type=str, action='store', default='KD02_R23', help='Which metallicity diagnostic to use (choose between KD02_R23,R23,R3,O3S2,O3O2,S2,R2,RS32,Te,P25,NB? Default is KD02_R23')
    parser.add_argument('--debug_Zdiag', dest='debug_Zdiag', action='store_true', default=False, help='Make additional plots to debug the metallicity diagnostic implementation? Default is no.')
    parser.add_argument('--ignore_combined_method', dest='ignore_combined_method', action='store_true', default=False, help='Ignore the combined method (S6 of KD02) while computing R23 metallicity and rely solely on R23? Default is no.')
    parser.add_argument('--exclude_lines', metavar='exclude_lines', type=str, action='store', default='', help='Which lines to be excluded for metallicity measurement with NB? Default is empty string, i.e., use all available lines')
    parser.add_argument('--dered_in_NB', dest='dered_in_NB', action='store_true', default=False, help='Make NebulaBayes de-redden the lines? Default is no (i.e., do de-reddening separately before calling NB)')
    parser.add_argument('--use_C25', dest='use_C25', action='store_true', default=False, help='Use the Cataldi+2025 calibrations rather than Curti+2020? Default is no.')

    parser.add_argument('--plot_AGN_frac', dest='plot_AGN_frac', action='store_true', default=False, help='Plot AGN fraction 2D map (based on BPT diagram)? Default is no.')
    parser.add_argument('--AGN_diag', metavar='AGN_diag', type=str, action='store', default='None', help='Which AGN-SF BPT-like diagnostic to use (choose between VO87,H21,O2O3,O2Hb,Ne3O2? Default is None')
    parser.add_argument('--mask_agn', dest='mask_agn', action='store_true', default=False, help='Mask out the AGN-dominated pixels from all metallicity estimates? Default is no.')
    
    parser.add_argument('--plot_circle_at_arcsec', metavar='plot_circle_at_arcsec', type=float, action='store', default=None, help='Radius in arcseconds of a circle to be plotted on every 2D map; default is None')
    parser.add_argument('--no_text_on_plot', dest='no_text_on_plot', action='store_true', default=False, help='Skip putting text annotations on plot2D? Default is no.')
    parser.add_argument('--plot_snr', dest='plot_snr', action='store_true', default=False, help='Plot the SNR map for a given 2D plot? Default is no.')

    parser.add_argument('--debug_Zsfr', dest='debug_Zsfr', action='store_true', default=False, help='Debug the metallicity-sfr plots? Default is no.')
    parser.add_argument('--histbycol', metavar='histbycol', type=str, action='store', default=None, help='Column name whose average per bin would be shown on line ratio histogram y-axis; Default is None, i.e., the usual counts will be on histogram y-axis')
    parser.add_argument('--clobber_mcmc', dest='clobber_mcmc', action='store_true', default=False, help='Over-write existing MCMC file? Default is no.')
    parser.add_argument('--fit_correlation', dest='fit_correlation', action='store_true', default=False, help='Fit a slope between x and y? Default is no.')

    # ------- wrap up and processing args ------------------------------
    args = parser.parse_args()
    if args.line_list != 'all': args.line_list = [item for item in args.line_list.split(',')]

    if args.id is not None: args.id = [int(item) for item in args.id.split(',')]

    if args.system == 'hd' and not os.path.exists('/Volumes/Elements/'): args.system = 'local'
    if args.line_list == 'all': args.line_list = ['Lya', 'OII', 'NeIII-3867', 'Hb', 'OIII-4363', 'OIII', 'Ha', 'NII','Ha+NII', 'SII', 'ArIII-7138', 'SIII', 'PaD','PaG','PaB','HeI-1083','PaA']

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

    return args

# ------------------------------------------------------------------------------------
def get_zrange_for_line(line, obs_wave_range=[800, 2200]):
    '''
    Computes the redshift range in which a certain emission line is available, for a given observed wavelength window (by default corresponds to NIRISS WFSS)
    Input wavelengths must be in nm
    Returns min and max redshift
    '''
    if type(line) in [float, int, np.float64, np.int64]: rest_wave = line
    elif type(line) in [str, np.str_]: rest_wave = rest_wave_dict[line]

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

# -------------------------------------------------------------------------------------------------------
def copy_from_hd_to_local(files_to_move=['*.txt', '*.png']):
    '''
    To copy all less heavy files (all txt and png) of each existing field in the HD to corresponding location locally
    '''
    hd_path = Path('/Volumes/Elements/acharyya_backup/Work/astro/passage/passage_output')
    local_path = Path('/Users/acharyya/Work/astro/passage/passage_output')

    available_fields = [os.path.split(item)[-1] for item in glob.glob(str(hd_path / 'Par*'))]

    for index, field in enumerate(available_fields):
        print(f'Doing field {field} which is {index+1} out of {len(available_fields)}..')
        dest_dir =local_path / field
        dest_dir.mkdir(parents=True, exist_ok=True)
        for files in files_to_move:
            target_files = str(hd_path / field / files)
            command = f'cp {target_files} {str(dest_dir)}/.'
            print(command)
            try: dummy = subprocess.check_output([command], shell=True)
            except: print('No such files. Skipping..')

    print('All done')

# -------------------------------------------------------------------------------------------------------
def get_fluxcols(args):
    '''
    Function to load or generate the list of filters and correpsonding flux columns in COSMOS2020 catalog
    Returns list of columns, and optionally the full cosmos2020 dataframe
    '''
    filepath = args.input_dir / 'COSMOS' / 'cosmos_fluxcols.npy'

    if os.path.exists(filepath):
        print(f'Reading flux columns from existing {filepath}')
        fluxcols = np.load(filepath)
        df_cosmos = None
    else:
        print(f'{filepath} does not exist, so preparing the list..')
        df_cosmos = read_COSMOS2020_catalog(args=args, filename=args.input_dir / 'COSMOS' / 'COSMOS2020_CLASSIC_R1_v2.2_p3.fits')

        all_flux_cols = [item for item in df_cosmos.columns if 'FLUX' in item and item != 'FLUX_RADIUS' and 'FLUXERR' not in item]
        filters = [item[:item.find('FLUX')] for item in all_flux_cols]
        fluxcols = [item + 'FLUX_AUTO' if item + 'FLUX_AUTO' in df_cosmos.columns else item + 'FLUX' for item in filters]
        fluxcols = list(dict.fromkeys(fluxcols)) # to remove duplicates
        np.save(filepath, fluxcols)

    return fluxcols, df_cosmos

# -------------------------------------------------------------------------------------------------------
def read_COSMOS2025_catalog(args=None, filename=None, aperture=1.0):
    '''
    Reads in the zCOSMOS galaxy catalog
    Returns as pandas dataframe
    '''
    aperture_dict = {0.2:0, 0.3:1, 0.5:2, 0.75:3, 1.0:4}
    if filename is None:
        if args is None: input_dir = '/Users/acharyya/Work/astro/passage/passage_data'
        else: input_dir = args.input_dir
        filename = Path(input_dir) / 'COSMOS' / 'COSMOSWeb_mastercatalog_v1_photom_primary.fits'

    print(f'Reading in {filename}, might take a while..')
    start_time2 = datetime.now()

    data = fits.open(filename)
    table = Table(data[1].data)

    multi_index_columns = [item for item in table.colnames if len(table[item].shape) > 1]
    single_index_columns = list(set(table.columns) - set(multi_index_columns))

    df  = table[single_index_columns].to_pandas()
    for thiscol in multi_index_columns: df[thiscol] = table[thiscol][:, aperture_dict[aperture]]

    print(f'Completed reading COSMOS2025 catalog in {timedelta(seconds=(datetime.now() - start_time2).seconds)}')

    return df

# -------------------------------------------------------------------------------------------------------
def read_COSMOSWebb_catalog(args=None, filename=None, aperture=1.0):
    '''
    Reads in the zCOSMOS galaxy catalog
    Returns as pandas dataframe
    '''
    aperture_dict = {0.1:0, 0.25:1, 0.5:2, 1.0:3, 1.5:4}
    if filename is None:
        if args is None: input_dir = '/Users/acharyya/Work/astro/passage/passage_data'
        else: input_dir = args.input_dir
        #filename = Path(input_dir) / 'COSMOS' / 'COSMOS_Web_for_Ayan_Dec24.fits'
        filename = Path(input_dir) / 'COSMOS' / 'COSMOSWeb_mastercatalog_v1_lephare_plus.fits'

    print(f'Reading in {filename}, might take a while..')
    start_time2 = datetime.now()

    table = Table.read(filename)

    multi_index_columns = [item for item in table.colnames if len(table[item].shape) > 1]
    single_index_columns = list(set(table.columns) - set(multi_index_columns))

    df  = table[single_index_columns].to_pandas()
    for thiscol in multi_index_columns: df[thiscol] = table[thiscol][:, aperture_dict[aperture]]

    df = df.rename(columns={'ID_SE++':'id', 'RA_DETEC':'ra', 'DEC_DETEC':'dec'})
    print(f'Completed reading COSMOSWebb catalog in {timedelta(seconds=(datetime.now() - start_time2).seconds)}')

    return df

# -------------------------------------------------------------------------------------------------------
def read_COSMOS2020_catalog(args=None, filename=None):
    '''
    Reads in the zCOSMOS galaxy catalog
    Returns as pandas dataframe
    '''
    if filename is None:
        if args is None: input_dir = '/Users/acharyya/Work/astro/passage/passage_data'
        else: input_dir = args.input_dir
        filename = Path(input_dir) / 'COSMOS' / 'COSMOS2020_CLASSIC_R1_v2.2_p3_subsetcolumns.fits'

    if not os.path.exists(filename) or (args is not None and args.clobber_cosmos): make_COSMOS_subset_table(filename, args)

    print(f'Reading in {filename}, might take a while..')
    start_time2 = datetime.now()

    data = fits.open(filename)
    table = Table(data[1].data)
    df  = table.to_pandas()
    df = df.rename(columns={'ID':'id', 'ALPHA_J2000':'ra', 'DELTA_J2000':'dec'})
    print(f'Completed reading COSMOS2020 catalog in {timedelta(seconds=(datetime.now() - start_time2).seconds)}')

    return df

# -------------------------------------------------------------------------------------------------------
def read_zCOSMOS_catalog(args=None, filename=None):
    '''
    Reads in the zCOSMOS galaxy catalog
    Returns as pandas dataframe
    '''
    if filename is None:
        if args is None: input_dir = '/Users/acharyya/Work/astro/passage/passage_data'
        else: input_dir = args.input_dir
        filename = input_dir / 'COSMOS/zCOSMOS-DR3' / 'zCOSMOS_VIMOS_BRIGHT_DR3_CATALOGUE.fits'

    data = fits.open(filename)
    table = Table(data[1].data)
    df  = table.to_pandas()
    df = df.rename(columns={'OBJECT_ID':'id', 'RAJ2000':'ra', 'DEJ2000':'dec'})

    return df

# -------------------------------------------------------------------------------------------------------
def make_COSMOS_subset_table(filename, args):
    '''
    Reads in the massive COSMOS2020 catalog and makes a smaller table with subset of columns and saves it
    '''
    suffix = '_subsetcolumns'
    filename = str(filename)
    if suffix in filename: filename = filename[:filename.find(suffix)] + '.fits'

    # -------determining flux column other columns to extract from df_cosmos-------
    fluxcols, _ = get_fluxcols(args)
    lp_cols_suffix = ['med', 'med_min68', 'med_max68', 'best']
    lp_cols = np.ravel([f'lp_{item}_{suffix}' for item in ['mass', 'SFR', 'sSFR'] for suffix in lp_cols_suffix])
    ez_cols_suffix = ['', '_p160', '_p500', '_p840']
    ez_cols = np.ravel([f'ez_{item}{suffix}' for item in ['mass', 'sfr', 'ssfr'] for suffix in ez_cols_suffix])
    flux_and_err_cols = np.ravel([[item, item.replace('FLUX', 'FLUXERR')] for item in fluxcols])
    cols_to_extract = np.hstack((['ID', 'ALPHA_J2000', 'DELTA_J2000', 'ID_COSMOS2015', 'ez_z_phot', 'lp_MK', 'lp_zBEST'], lp_cols, ez_cols, flux_and_err_cols)).tolist()

    print(f'Trying to read in  {filename}; can take a while..')
    data = fits.open(filename)
    table = Table(data[1].data)
    table_sub = table[cols_to_extract]

    outfilename = str(filename).split('.fits')[0] + suffix + '.fits'
    table_sub.write(outfilename, overwrite=True)
    print(f'Saved subset table as {outfilename}')

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

# -------------------------------------------------------------------------------------------------------
def split_COSMOS_subset_table_by_par(args, old_dir_format=True, fields=None):
    '''
    Reads in the subset of columns of COSMOS2020 catalog and splits it into smaller tables with only objects that are overlapping with individual PASSAGE fields
    '''
    # -------reading in the COSMOS2020 (sub)catalog------
    filename = Path(args.input_dir) / 'COSMOS' / 'COSMOS2020_CLASSIC_R1_v2.2_p3_subsetcolumns.fits'
    data = fits.open(filename)
    table_cosmos = Table(data[1].data)

    df_cosmos = table_cosmos.to_pandas()
    df_cosmos = df_cosmos.rename(columns={'ALPHA_J2000':'ra', 'DELTA_J2000':'dec'})

    if fields is None:
        field_list = [os.path.split(item[:-1])[1] for item in glob.glob(str(args.input_dir / args.drv / 'Par*') + '/')]
        field_list += [f'Par{item:03d}' for item in passage_fields_in_cosmos]
        field_list = list(np.unique(field_list))
        field_list.sort(key=natural_keys)
    else:
        field_list = fields

    for index, thisfield in enumerate(field_list):
        print(f'Starting {index+1} of {len(field_list)} fields..')
        # -------determining path to photometric catalog------
        if old_dir_format: product_dir = args.input_dir / thisfield / 'Products'
        else: product_dir = args.input_dir / 'data' / thisfield / 'DATA' / 'DIRECT_GRISM'
        catalog_file = product_dir / f'{thisfield}_photcat.fits'

        if os.path.exists(catalog_file):
            # -------reading in photometric catalog------
            catalog = GTable.read(catalog_file)
            df = catalog['id', 'ra', 'dec'].to_pandas()
            df['passage_id'] = thisfield + '-' + df['id'].astype(str)  # making a unique combination of field and object id

            # -------cross-matching RA/DEC of both catalogs------
            df_crossmatch = get_crossmatch(df, df_cosmos, sep_threshold=0.1, df1_idcol='passage_id', df2_idcol='ID')
            df_crossmatch = df_crossmatch.rename(columns={'df1_id': 'passage_id', 'df2_id': 'ID'})

            if len(df_crossmatch) > 0:
                table_crossmatch = Table.from_pandas(df_crossmatch)
                table_cosmos_thisfield = join(table_cosmos, table_crossmatch, keys='ID')

                if old_dir_format: outfilepath = args.input_dir / 'COSMOS'
                else: outfilepath = product_dir
                outfilename = outfilepath / f'cosmos2020_objects_in_{thisfield}.fits'
                
                table_cosmos_thisfield.write(outfilename, overwrite=True)
                print(f'Saved subset table as {outfilename}')
            else:
                print(f'No overlapping objects found in field {thisfield}')
        else:
            print(f'{catalog_file} does not exist, so skipping {thisfield}.')

# -------------------------------------------------------------------------------------------------------
def split_COSMOSWebb_table_by_par(args, filename=None, old_dir_format=True, fields=None):
    '''
    Reads in the COSMOSWebb catalog and splits it into smaller tables with only objects that are overlapping with individual PASSAGE fields
    '''
    # -------reading in the COSMOS2020 (sub)catalog------
    #if filename is None: filename = Path(args.input_dir) / 'COSMOS' / 'COSMOS_Web_for_Ayan_Dec24.fits'
    if filename is None: filename = Path(args.input_dir) / 'COSMOS' / 'COSMOSWeb_mastercatalog_v1_lephare_plus.fits'
    data = fits.open(filename)
    table_cosmos = Table(data[1].data)

    cosmos_idcol = 'id'
    df_cosmos = table_cosmos[[cosmos_idcol, 'ra', 'dec']].to_pandas()

    if fields is None:
        field_list = [os.path.split(item[:-1])[1] for item in glob.glob(str(args.input_dir / args.drv / 'Par*') + '/')]
        field_list += [f'Par{item:03d}' for item in passage_fields_in_cosmos]
        field_list = list(np.unique(field_list))
        field_list.sort(key=natural_keys)
    else:
        field_list = fields

    # The following method is slower because it involves cross matching all COSMOS objects with..
    # ..objects in ParXXX, but it leads to an accurate list of objects within ParXXX.
    for index, thisfield in enumerate(field_list):
        print(f'Starting {index+1} of {len(field_list)} fields..')
        # -------determining path to photometric catalog------
        if old_dir_format: product_dir = args.input_dir / thisfield / 'Products'
        else: product_dir = args.input_dir / 'data' / thisfield / 'DATA' / 'DIRECT_GRISM'
        catalog_file = product_dir / f'{thisfield}_photcat.fits'

        if os.path.exists(catalog_file):
            # -------reading in photometric catalog------
            catalog = GTable.read(catalog_file)
            df = catalog['id', 'ra', 'dec'].to_pandas()
            df['passage_id'] = thisfield + '-' + df['id'].astype(str)  # making a unique combination of field and object id

            # -------cross-matching RA/DEC of both catalogs------
            df_crossmatch = get_crossmatch(df, df_cosmos, sep_threshold=0.1, df1_idcol='passage_id', df2_idcol=cosmos_idcol)
            df_crossmatch = df_crossmatch.rename(columns={'df1_id': 'passage_id', 'df2_id': cosmos_idcol})

            if len(df_crossmatch) > 0:
                table_crossmatch = Table.from_pandas(df_crossmatch)
                table_cosmos_thisfield = join(table_cosmos, table_crossmatch, keys=cosmos_idcol)
                
                if old_dir_format: outfilepath = args.input_dir / 'COSMOS'
                else: outfilepath = product_dir
                outfilename = outfilepath / f'cosmoswebb_objects_in_{thisfield}.fits'
                
                table_cosmos_thisfield.write(outfilename, overwrite=True)
                print(f'Saved subset table as {outfilename}')
            else:
                print(f'No overlapping objects found in field {thisfield}')
        else:
            print(f'{catalog_file} does not exist, so skipping {thisfield}.')

    '''
    # The following method is faster because it does not do cross matching, just takes all objects contained..
    # ..within the region file corresponding to ParXXX, but since the region files are approximate..
    # ..the resulting list of COSMOS objects "contained" within ParXXX is not accurate either.
    for index, thisfield in enumerate(field_list):
        print(f'Starting {index+1} of {len(field_list)} fields..')
        # -------determining path to direct images------
        product_dir = args.input_dir / args.drv / thisfield / 'Products'
        image_files = glob.glob(str(product_dir) + f'/{thisfield}*clear*sci*.fits')

        if len(image_files) > 0:
            image_file = Path(image_files[0])
            hdul = fits.open(image_file)
            source_wcs = pywcs.WCS(hdul[0].header)

            region_dir = image_file.parent.parent / 'Regions'
            region_dir.mkdir(parents=True, exist_ok=True)
            region_file = region_dir / Path(filename.stem + '.reg')
            source_wcs.footprint_to_file(region_file, color='cyan', width=1)

            sky_region = Regions.read(region_file, format='ds9')[0]
            contained_ids = sky_region.contains(SkyCoord(df_cosmos['ra'], df_cosmos['dec'], unit='deg'), pywcs.WCS(hdul[0].header))
            table_contained = table_cosmos[contained_ids]

            if len(table_contained) > 0:
                outfilename = args.input_dir / 'COSMOS' / args.drv / f'cosmoswebb_objects_in_{thisfield}.fits'
                table_contained.write(outfilename, overwrite=True)
                print(f'Saved subset table as {outfilename}')
            else:
                print(f'No overlapping objects found in field {thisfield}')
        else:
                print(f'No image file exists, so skipping {thisfield}.')
    '''

# -------------------------------------------------------------------------------------------------------
def get_passage_filter_dict(args=None, filename=None):
    '''
    Reads PASSAGE spreadsheet and returns a dictionary with PASSAGE field names and which filters are present in them
    '''
    if filename is None:
        filename =  args.root_dir / 'passage_data' / 'v0.5' / 'JWST PASSAGE Cycle 1 - Cy1 Executed.csv'

    df = pd.read_csv(filename)
    df = df[['Par#', 'Obs Date', 'Filter']]
    df = df.ffill()
    df = df[~ df['Obs Date'].str.contains('SKIPPED')]

    dictionary = {item: np.unique(df[df['Par#'] == item]['Filter']).tolist() for item in np.unique(df['Par#'])}
    dictionary.update({'glass-a2744':['F115W', 'F150W', 'F120W'],
                        'Par999':['F115W', 'F150W', 'F120W']}) # for testing with mock data (Par999)

    return dictionary

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
def get_files_in_url(url, ext='fits', auth=None):
    '''
    Function to get list of fits files in a given URL
    Returns list of urls of the files
    From https://stackoverflow.com/questions/11023530/python-to-list-http-files-and-directories
    '''
    if auth is None: page = requests.get(url).text
    else: page = requests.get(url, auth=auth).text
    soup = BeautifulSoup(page, 'html.parser')
    url_list = [url + '/' + node.get('href') for node in soup.find_all('a') if node.get('href').endswith(ext)]

    return url_list

# ------------------------------------------------------------------------------------------------------
def download_files_from_url(url, outdir, ext='fits', match_strings=[''], auth=None):
    '''
    Function to download all fits files that contain <amch_string> in the filename from a given URL
    Saves the downloaded files in outdir
    '''
    start_time = datetime.now()
    outdir = Path(outdir)
    url_list = get_files_in_url(url, ext=ext, auth=auth)
    len_orig = len(url_list)
    for match_string in match_strings: url_list = [item for item in url_list if match_string in item]
    print(f'Found {len(url_list)} matching files out of {len_orig} at URL {url}')

    for index, this_url in enumerate(url_list):
        target_file = os.path.split(this_url)[-1]
        if os.path.isfile(outdir / target_file) or os.path.isfile(outdir / Path(target_file).stem):
            print(f'Skipping file {target_file} ({index + 1} of {len(url_list)}) because it already exists at target location')
        else:
            start_time2 = datetime.now()
            print(f'Downloading file {target_file} which is {index + 1} of {len(url_list)}..')
            if auth is None: response = requests.get(this_url)
            else: response = requests.get(this_url, auth=auth)
            with open(outdir / target_file, 'wb') as file: file.write(response.content)
            if target_file.endswith('.gz'):
                print(f'Unzipping downloaded file..')
                unzip_and_delete(outdir / target_file, outdir)
            print(f'Completed this download in {timedelta(seconds=(datetime.now() - start_time2).seconds)}')


    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')

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

# --------------------------------------------------------------------------------------------------------------------
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

# --------------------------------------------------------------------------------------------------------------------
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

# --------------------------------------------------------------------------------------------------------------------
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

# --------------------------------------------------------------------------------------------------------------------
def move_field_after_download(field_arr, args=None):
    '''
    Unzips and moves the reduced data from Box from Downloads folder to appropriate location by making the right dorectories
    '''
    origin_dir = HOME / 'Downloads'
    field_arr = np.atleast_1d(field_arr)

    for index, field in enumerate(field_arr):
        start_time = datetime.now()
        if 'Par' in str(field): field = f'Par{int(field.split("Par")[1]):03d}'
        else: field = f'Par{field:03d}'
        print(f'Starting field {field} which is {index + 1} out of {len(field_arr)}..')

        try:
            if args is None: destination_dir = Path('/Volumes/Elements/acharyya_backup/Work/astro/passage/passage_data/v0.5/data')
            else: destination_dir = args.input_dir / 'v0.5/data'
            Path(destination_dir / 'linelist').mkdir(exist_ok=True, parents=True)

            unzip_and_delete(origin_dir / f'{field}.zip', destination_dir)

            Path(destination_dir / field / 'DATA/DIRECT_GRISM').mkdir(exist_ok=True, parents=True)
            unzip_and_delete(destination_dir / field / f'{field}_spec1D.tar.gz', destination_dir / field)
            unzip_and_delete(destination_dir / field / f'{field}_spec2D.tar.gz', destination_dir / field)

            cat_files = glob.glob(str(destination_dir / field / '*cat.fits'))
            for this_cat_file in cat_files:
                print(f'Moving {os.path.split(this_cat_file)[-1]}..')
                shutil.move(this_cat_file, destination_dir / field / 'DATA/DIRECT_GRISM')

            fits_files = glob.glob(str(destination_dir / field / '*.fits'))
            for this_fits_file in fits_files:
                print(f'Moving {os.path.split(this_fits_file)[-1]}..')
                shutil.move(this_fits_file, destination_dir / field / 'DATA')

            filename = f'{field}lines.dat'
            print(f'Moving {filename}..')
            shutil.move(destination_dir / field / filename, destination_dir / 'linelist')
        except Exception as e:
            print(f'Probably skipping {field} due to following error: {e}')
            continue

        print(f'Completed moving {field} in {timedelta(seconds=(datetime.now() - start_time).seconds)}')

# --------------------------------------------------------------------------------------------------
def adjust_lightness(color, amount=0.5):
    '''
    Lightens (if amount > 1) or darkens (if amount < 1) matplitlib colors
    Adapted from https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib
    '''
    c = colorsys.rgb_to_hls(*mplcolors.to_rgb(color))
    color = colorsys.hls_to_rgb(c[0], max(0.1, min(0.9, amount * c[1])), c[2])
    return color

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

# --------------------------------------------------------------------------------------------------
def compare_SNR(filename, line_label, radius_arcsec = 0.25):
    '''
    Test function to compare various measures of SNRs of a given line in a given full.fits filename
    Returns integrated line flux, line flux map within segmentation map, and line flux map within a certain distance form the center, all with accompanying uncertainties
    '''
    full_hdu = fits.open(filename)
    header = full_hdu[0].header

    available_lines = np.array(header['HASLINES'].split(' '))
    line_index = np.where(available_lines == line_label)[0][0]

    line_flux_int = header[f'FLUX{line_index + 1:03d}']
    line_err_int = header[f'ERR{line_index + 1:03d}']
    line_int = ufloat(line_flux_int, line_err_int)

    line_flux_map = full_hdu['LINE', line_label].data * 1e-17
    line_err_map = 1e-17 / (full_hdu['LINEWHT', line_label].data ** 0.5)
    line_map = unp.uarray(line_flux_map, line_err_map)

    bad_mask = unp.nominal_values(line_map) <= 0 | ~np.isfinite(unp.std_devs(line_map))
    seg_map = full_hdu['SEG'].data
    seg_mask = seg_map != header['ID']

    line_map_seg = np.ma.masked_where(bad_mask | seg_mask, line_map)
    line_map_seg_sum = np.sum(line_map_seg)

    line_wcs = pywcs.WCS(full_hdu['DSCI'].header)
    pix_size_arcsec = utils.get_wcs_pscale(line_wcs)

    pixscale_kpc = pix_size_arcsec/ cosmo.arcsec_per_kpc_proper(header['REDSHIFT']).value # kpc
    center_pix = np.shape(line_map)[0] / 2.
    distance_map = np.array([[np.sqrt((i - center_pix)**2 + (j - center_pix)**2) for j in range(np.shape(line_map)[1])] for i in range(np.shape(line_map)[0])]) * pixscale_kpc # kpc
    radius_kpc = radius_arcsec / cosmo.arcsec_per_kpc_proper(header['REDSHIFT']).value  # converting arcsec to kpc
    distance_mask = distance_map > radius_kpc

    line_map_center = np.ma.masked_where(bad_mask | distance_mask, line_map)
    line_map_center_sum = np.sum(line_map_center)

    print(f'Obj: {header["ID"]}, line: {line_label}\n'
          f'integrated flux = {line_int}, SNR={line_int.n / line_int.s: .1f}\n'
          f'total sum within seg map = {line_map_seg_sum: .2e}, SNR={line_map_seg_sum.n / line_map_seg_sum.s: .1f}\n'
          f'total sum within {radius_arcsec: .2f} arcsec (={radius_kpc: .1f} kpc) radius = {line_map_center_sum: .2e}, SNR={line_map_center_sum.n / line_map_center_sum.s: .1f}\n')

    return line_int, line_map_seg, line_map_center

# --------------------------------------------------------------------------------------------------------------------
def make_rgb(im_r, im_g, im_b, interval=None, stretch=None):
    """
    This function was borrowed from ChatGPT
    Make an RGB image from three input images with customizable interval and stretch.
    
    Parameters
    ----------
    im_r, im_g, im_b : array-like
        Input images for red, green, blue channels.
    interval : instance of astropy.visualization.Interval, optional
        Interval object to define min/max scaling (e.g., ManualInterval, MinMaxInterval).
    stretch : instance of astropy.visualization.Stretch, optional
        Stretch object to define stretching function (e.g., LinearStretch, LogStretch, SqrtStretch).
    
    Returns
    -------
    rgb : ndarray
        (M, N, 3) RGB image with values scaled between 0 and 1.
    """
    # Default behavior if interval or stretch not provided
    if interval is None:
        interval = MinMaxInterval()
    if stretch is None:
        stretch = LinearStretch()

    def normalize(im):
        vmin, vmax = interval.get_limits(im)
        im_scaled = (im - vmin) / (vmax - vmin)
        im_scaled = np.clip(im_scaled, 0, 1)
        return stretch(im_scaled)

    r = normalize(im_r)
    g = normalize(im_g)
    b = normalize(im_b)

    rgb = np.stack([r, g, b], axis=-1)
    rgb = np.clip(rgb, 0, 1)
    return rgb

# --------------------------------------------------------------------------------------------------
def parse_latex_value(entry):
    '''
    This function was borrowed from ChatGPT
    Function to parse LaTeX string
    '''
    match = re.match(r"\$(?P<val>[\d\.]+)_{-(?P<minus>[\d\.]+)}\^{\+(?P<plus>[\d\.]+)}\$", entry)
    if match:
        value = float(match.group('val'))
        minus = float(match.group('minus'))
        plus = float(match.group('plus'))
        symmetric_error = (minus + plus) / 2
        return value, symmetric_error
    else:
        return None, None  # or raise an error if you prefer

# --------------------------------------------------------------------------------------------------------------
def trim_image(image, args=None, arcsec_limit=None, pix_size_arcsec=None, re_limit=None, re_arcsec=None, skip_re_trim=False):
    '''
    Trim a given 2D image to a given arcsecond dimension
    Returns 2D map
    '''
    if args is not None:
        pix_size_arcsec = args.pix_size_arcsec
        if args.re_limit is None or skip_re_trim:
            arcsec_limit = args.arcsec_limit
        else:
            arcsec_limit = args.re_limit * args.re_arcsec
    elif re_limit is not None and re_arcsec is not None:
        arcsec_limit = re_limit * re_arcsec

    image_shape = np.shape(image)
    center_pix = int(image_shape[0] / 2.)
    farthest_pix = int(arcsec_limit / pix_size_arcsec) # both quantities in arcsec

    image = image[center_pix - farthest_pix : center_pix + farthest_pix, center_pix - farthest_pix : center_pix + farthest_pix]
    # print(f'Trimming image of original shape {image_shape} to {args.arcsec_limit} arcseconds, which is from pixels {center_pix - farthest_pix} to {center_pix + farthest_pix}, so new shape is {np.shape(image)}')

    return image

# --------------------------------------------------------------------------------------------------
def extract_emission_line_map_from_full(line, args=None, field=None, id=None, arcsec_limit=1.):
    '''
    Reads in the fill.fits or maps.fits, as available, then reads in the emission line flux and err map
    Writes out the two maps as two separate fits files
    Reads the newly written files and returns as two 2D numpy arrays
    '''
    # -------reading full.fits---------------
    if id is not None:
        if field is None:
            if args is not None: field = args.field
            else: field = 'Par028'
        if args is not None:
            full_filename = args.input_dir / field / 'Products' / 'full' / f'{field}_{id:05d}.full.fits'
        else:
            full_filename = HOME / 'Work/astro/passage/passage_data/v0.5' / field / 'Products' / 'full' / f'{field}_{id:05d}.full.fits'
    else:
        full_filename = args.input_dir / args.field / 'Products' / 'full' / f'{args.field}_{args.id:05d}.full.fits'
    if not os.path.exists(full_filename): full_filename = Path(str(full_filename).replace('full', 'maps'))
    print(f'Reading in {full_filename}..')
    full_hdu = fits.open(full_filename)

    # -----reading in required emission line-------
    try:
        line_hdu = full_hdu['LINE', line]
        line_map_err = 1e-17 / (full_hdu['LINEWHT', line].data ** 0.5)  # ERR = 1/sqrt(LINEWHT) = flux uncertainty; in units of ergs/s/cm^2
    except KeyError:
        if line == 'OIII':
            line_hdu = full_hdu['LINE', 'OIII-5007']
            line_map_err = 1e-17 / (full_hdu['LINEWHT', 'OIII-5007'].data ** 0.5)  # ERR = 1/sqrt(LINEWHT) = flux uncertainty; in units of ergs/s/cm^2
 
    line_map = line_hdu.data * 1e-17 # in units of ergs/s/cm^2

    # --------trimming the images------------
    line_wcs = pywcs.WCS(full_hdu['DSCI'].header)
    pix_size_arcsec = utils.get_wcs_pscale(line_wcs)
    line_map = trim_image(line_map, pix_size_arcsec=pix_size_arcsec, arcsec_limit=args.arcsec_limit if args is not None else arcsec_limit)
    line_map_err = trim_image(line_map_err, pix_size_arcsec=pix_size_arcsec, arcsec_limit=args.arcsec_limit if args is not None else arcsec_limit)

    # -----writing out required emission line-------
    map_hdul = fits.HDUList(fits.PrimaryHDU(data=line_map))
    err_hdul = fits.HDUList(fits.PrimaryHDU(data=line_map_err))

    if args is not None: outfilename = args.input_dir / args.field / 'emission_maps' / 'unbinned' / f'{args.field}_{args.id:05d}_{line}_map.fits'
    else: outfilename = HOME / 'Work/astro/passage/passage_data/v0.5' / field / 'emission_maps' / 'unbinned' / f'{field}_{id:05d}_{line}_map.fits'
    
    map_hdul.writeto(outfilename, overwrite=True)
    err_hdul.writeto(str(outfilename).replace('map.fits', 'map_u.fits'), overwrite=True)
    print(f'Saved as {outfilename}')

    data_map = fits.open(outfilename)[0].data
    data_err = fits.open(str(outfilename).replace('map.fits', 'map_u.fits'))[0].data

    return data_map, data_err

# ---------------------------------------------------------------------------------------------------------------------------------------------
def plotline(id, path='/Users/acharyya/Work/astro/passage/passage_data/v0.5/Par028/Products/maps/', line='Ha', arc=1, cmin=None, cmax=None, savefig=False, hide=False):
     hdul = fits.open(path + f'Par028_{id:05d}.maps.fits')
     header = hdul[0].header
     z = header['REDSHIFT']
     line_hdu = hdul['LINE', line]
     pix_size_arcsec = line_hdu.header['PIXASEC']
     line_map = line_hdu.data * 1e-17
     line_map = trim_image(line_map, arcsec_limit=arc, pix_size_arcsec=pix_size_arcsec)
     fig, ax = plt.subplots(1)
     im = ax.imshow(np.log10(line_map), cmap='plasma', extent = [-1, 1, -1, 1], vmin=cmin, vmax=cmax)
     ax.set_title(f'#{id}; z={z:.1f}')
     ax.set_xlabel('arcsec')
     ax.set_ylabel('arcsec')
     cbar = plt.colorbar(im)
     cbar.set_label(f'Log {line} flux')
     if savefig:
         figdir = Path('/Users/acharyya/Downloads/linemaps')
         figdir.mkdir(exist_ok=True, parents=True)
         figname = figdir / f'Par028_{id:05d}_{line}_map.png'
         fig.savefig(figname)
         print(f'Saved as {figname}')
     if hide: plt.close()
     else: plt.show(block=False)
     return fig

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
def annotate_PAs(pa_arr, ax, fontsize=10, color='k'):
    '''
    Annotates a given plot with the PAs of all available filters in the given axis
    Returns the axis handle
    '''
    x_cen = ax.get_xlim()[1] - 0.15 * np.diff(ax.get_xlim())[0]
    y_cen = ax.get_ylim()[1] - 0.15 * np.diff(ax.get_ylim())[0]
    len = np.diff(ax.get_xlim())[0] * 0.1

    for pa in pa_arr:
        x_comp = len * np.sin(pa * np.pi / 180)
        y_comp = len * np.cos(pa * np.pi / 180)
        ax.plot([x_cen, x_cen - x_comp], [y_cen, y_cen + y_comp], lw=1, c=color)
        ax.text(x_cen - 1.2 * x_comp, y_cen + 1.2 * y_comp, r'%d$^\circ$' % pa, color=color, fontsize=fontsize, ha='center', va='center', rotation=pa)

    return ax

# --------------------------------------------------------------------------------------------------------------
def plot_rgb_image(objid, images_dir, filters=['F200W', 'F150W', 'F115W'], arcsec_limit=1, fontsize=15, hide_filter_names=False, hide_pa=False, offset=0, skip_annotate=False, make_transparent=False):
    '''
    Plots the RGB cutout of a given object ID (objid) by cutting out available filter images in the provided path (images_dir)
    Returns figure
    '''
    PSF_FWHM_arcsec_dict = {'F115W': 0.094, 'F150W': 0.092, 'F200W': 0.106} # value of PSF FWHM derived from Table 1 of https://jwst-docs.stsci.edu/jwst-near-infrared-imager-and-slitless-spectrograph/niriss-performance/niriss-point-spread-functions#gsc.tab=0
    col_arr = ['r', 'lightgreen', 'cornflowerblue']
    hdu_filename = Path(glob.glob(f'{str(Path(images_dir))}/spec2D/*_{objid:05d}.spec2D.fits')[0])
    print(f'Reading in {hdu_filename}..')

    # ------object specific quantitites-------
    field = hdu_filename.stem.split('_')[0]
    hdul = fits.open(hdu_filename)
    pos = SkyCoord(hdul[0].header['RA'], hdul[0].header['DEC'], unit = 'deg')
    size = 2 * arcsec_limit * u.arcsec
    filter_ind = 0
    pa_arr = [hdul[0].header[f'{filters[filter_ind]}{item + 1:02d}'] for item in range(hdul[0].header[f'N{filters[filter_ind]}'])]
    
    print(f'Plotting the RGB images with filters {filters} and PAs {pa_arr}..')
    fig, ax = plt.subplots(1, figsize=(6.5, 6))

    # -------get image for each filter---------
    image_arr, exptime_arr = [], []
    for index, filter in enumerate(filters):
       #filename = Path(images_dir) / f'{field}_{filter.lower()}-clear_drz_sci.fits'
        filename = Path(images_dir) / f'{field}_{filter.lower()}_drz_sci.fits'
        data = fits.open(filename)
        exptime = data[0].header['EXPTIME']
        exptime_arr.append(exptime)
        
        print(f'Reading {filter} image from {filename}, which has exptime={exptime} s..')

        image = data[0].data
        source_header = data[0].header
        wcs = pywcs.WCS(source_header)

        cutout = Cutout2D(image, pos, size, wcs=wcs)
        cutout_data = cutout.data * source_header['PHOTFNU'] * 1e6
        cutout_header = cutout.wcs.to_header()
        if index:
            reprojected_data, reprojected_header = reproject_interp((cutout_data, cutout_header), target_header, shape_out=image_arr[0].shape)
            image_arr.append(reprojected_data)
        else:
            target_header = cutout_header
            image_arr.append(cutout_data)

        if not hide_filter_names:
            ax.text(0.05, 0.95 - index * 0.1, f'{filter}', c=col_arr[index], fontsize=fontsize, ha='left', va='top', transform=ax.transAxes, zorder=10)
            if filter == 'F150W': ax.add_patch(plt.Circle((0.35, 0.93 - index * 0.1), PSF_FWHM_arcsec_dict[filter]/2, color=col_arr[index], fill=True, lw=0.5, transform=ax.transAxes)) # for over-plotting PSF size

    # -------normalising each image to exptime---------
    for index in range(len(exptime_arr)):
        image_arr[index] = image_arr[index] * exptime_arr[0] / exptime_arr[index]

    # -------create RGB image---------
    pctl, maximum = 99.9, 0.
    for img in image_arr:
        val = np.percentile(img, pctl)
        if val > maximum: maximum = val

    rgb_image = make_rgb(image_arr[0], image_arr[1], image_arr[2], interval=ManualInterval(vmin=0, vmax=maximum), stretch=SqrtStretch())

    # -------plot RGB image---------
    extent = (-arcsec_limit - offset, arcsec_limit - offset, -arcsec_limit - offset, arcsec_limit - offset)
    p = ax.imshow(rgb_image, origin='lower', extent=extent, alpha=1)
    ax.set_aspect('auto') 
    ax.scatter(0, 0, marker='x', s=10, c='grey')

    # ----------annotate axis---------------
    if not hide_pa:
        ax = annotate_PAs(pa_arr, ax, fontsize=fontsize, color='w')

    ax.set_xlim(-arcsec_limit, arcsec_limit)  # arcsec
    ax.set_ylim(-arcsec_limit, arcsec_limit)  # arcsec

    if not skip_annotate: ax = annotate_axes(ax, 'arcsec', 'arcsec', fontsize=fontsize, hide_xaxis=False, hide_yaxis=False, hide_cbar=True)

    plot_output_dir = Path(images_dir) / 'figs'
    plot_output_dir.mkdir(exist_ok=True, parents=True)
    figname = plot_output_dir / f'{field}_{objid:05d}_RGB.png'
    fig.savefig(figname, transparent=make_transparent)
    print(f'\nSaved figure as {figname}')
    plt.show(block=False)

    return fig

# --------------------------------------------------------------------------------------------------------------------
def get_passage_masses_from_cosmos(df, args, id_col='objid', field_col='field', cosmos_idcol='id'):
    '''
    Derives stellar masses of PASSAGE galaxies present in the given dataframe from COSMOS-Web catalog
    Returns dataframe
    '''
    passage_fields = [item for item in np.unique(df[field_col]) if 'Par' in item]
    
    df_cosmos = pd.DataFrame()
    for index, thisfield in enumerate(passage_fields):
        cosmosfilename = args.input_dir / 'COSMOS' /  f'cosmoswebb_objects_in_{thisfield}.fits'
        df_cosmos_thisfield = Table.read(cosmosfilename).to_pandas()
        sed_cols_to_extract = ['passage_id', cosmos_idcol, 'ra', 'dec', 'mass_med', 'sfr_med', 'ssfr_med']
        df_cosmos = pd.concat([df_cosmos, df_cosmos_thisfield[sed_cols_to_extract]])
    
    df_cosmos = df_cosmos.rename(columns={cosmos_idcol: 'cosmos_id'})
    df_cosmos['passage_id'] = df_cosmos['passage_id'].astype(str)

    df['passage_id'] = df[field_col].astype(str) + '-' + df[id_col].astype(str)  # making a unique combination of field and object id
    df = pd.merge(df, df_cosmos, on=['passage_id'], how='right')
    df = df.drop('passage_id', axis=1)
    df = df.rename(columns={'mass_med': 'log_mass', 'sfr_med': 'log_sfr', 'ssfr_med': 'log_ssfr'})
    df = df[(df['log_mass'] > 0) & (df['log_sfr'] > 0)]
    df = df.dropna(subset=['log_mass', 'log_sfr'], axis=0)

    return df

# --------------------------------------------------------------------------------------------------------------------
def save_fig(fig, fig_dir, figname, args):
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
    print(f'\nSaved figure as {figname}')
    plt.show(block=False)

    return

