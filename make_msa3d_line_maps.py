'''
    Filename: make_msa3d_line_maps.py
    Notes: Performs emission line fitting on MSA-3D datacubes and stores them in fits files
    Author : Ayan
    Created: 06-02-26
    Example: run make_msa3d_line_maps.py --do_all_obj
             run make_msa3d_line_maps.py --id 2111
'''

from header import *
from util import *

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def get_ifu_cube(cube_fits_file, err_fits_file=None):
    '''
    Reads in the IFU datacube from a given filename
    Also reads in the corresponding errorcube, if that filename is provided
    Returns both cubes as 3D numpy arrays as well as the wavelength array
    '''
    print(f'\tReading in IFU cube from {cube_fits_file}..')
    cube = fits.open(cube_fits_file)[0].data
    if err_fits_file is None: cube_err = None
    else: cube_err = fits.open(err_fits_file)[0].data

    

    return cube, cube_err, wave_arr

# --------------------------------------------------------------------------------------------------------------------
def linefit_cube(cube, wave_arr, args, cube_err=None):
    '''
    Runs spaxel-by-spaxel emission line fitting of the emission line list provided (as args.line_list)
    Accounts for the corresponding errorcube, if provided
    Returns N x 2D emission line maps, including corresponding uncertainties
    '''
    print(f'\Fitting cube for ID {args.id}..')
    img_shape = np.shape(cube)[1:]
    linemap_arr = np.zeros(len(args.line_list), img_shape[0], img_shape[1])
    wave_to_fit = [rest_wave_dict[item] for item in args.line_list]

    rest_wave = wave_arr / (1 + args.z) # converting to rest-frame wavelength, in nm

    for i in range(img_shape[0]):
        for j in range(img_shape[1]):
            if not args.silent: print(f'\t\tFitting flux in cell ({i},{j})..')
            flux = cube[:, i, j]
            if cube_err is None: flux_err = np.zeros(np.shape(flux))
            else: flux_err = cube_err[:, i, j]

            lineflux_arr = linefit_spaxel(flux, rest_wave, wave_to_fit, flux_err=flux_err)
            linemap_arr[:, i, j] = lineflux_arr

    return linemap_arr

# --------------------------------------------------------------------------------------------------------------------
def linefit_spaxel(flux, rest_wave, wave_to_fit, flux_err=None):
    '''
    Runs emission line fitting on one spectrum (spaxel) of the emission line wavelengths provided (wave_to_fit)
    Accounts for the corresponding error spectrum, if provided
    Returns N emission line fluxes, including corresponding uncertainties
    '''

    return 0#lineflux_arr

# --------------------------------------------------------------------------------------------------------------------
def save_line_maps(linemap_arr, maps_fits_file, args):
    '''
    Saves the N x 2D emission line maps as N-extension fits files, along with other relevant header info
    '''

    return

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    args.fontfactor = 1.5
    args.id_arr = args.id

    # -------------setup directories and global variables----------------
    cube_fits_dir = args.input_dir / 'fluxcubes'
    err_fits_dir = args.input_dir / 'errcubes'
    maps_fits_dir = args.output_dir / 'maps'

    catalog_file = args.input_dir / 'redshifts.dat'

    # ----------------reading in catalog---------------------
    columns = ['id', 'ra', 'dec', 'type', 'redshift', 'mass', 'sfr', 'ssfr', 'Zgrad_kpc', 'Zgrad_kpc_u', 'Zgrad_re', 'Zgrad_re_u', 'num']
    df = pd.read_csv(catalog_file, delim_whitespace=True, comment='#', names=columns)

    if not args.do_all_obj: df = df[df['id'].isin(args.id_arr)]

    # ----------------looping over the objects in this chunk-------------
    for index, obj in df.iterrows():
        start_time2 = datetime.now()
        print(f'Commencing ({index + 1}/{len(df)}) ID {obj["id"]}..')
        args.id =obj['id']
        args.redshift = obj['redshift']

        # ------determining directories and filenames---------
        args.cube_fits_file = cube_fits_dir / f'cube_square_medians_{args.id:d}_hdr_rect.fits'
        args.err_fits_file = err_fits_dir / f'cube_square_medians_{args.id:d}_hdr_rect_err.fits'
        args.maps_fits_file = maps_fits_dir / f'{args.id:05d}.maps.fits'

        # -----------read in the cube--------------
        cube, cube_err, wave_arr = get_ifu_cube(args.cube_fits_file, err_fits_file=args.err_fits_file)

        # -----------emission line fitting--------------
        linemap_arr = linefit_cube(cube, wave_arr, args, cube_err=cube_err)

        # -----------save the emission maps in fits file-------------
        save_line_maps(linemap_arr, args.maps_fits_file, args)


        print(f'\nCompleted ID {args.id} in {timedelta(seconds=(datetime.now() - start_time2).seconds)}, {len(df) - index - 1} to go!')

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
