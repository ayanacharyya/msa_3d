'''
    Filename: compare_kinematic_maps.py
    Notes: Compares kinematic maps from Ayan's emission line fitting vs others'
    Author : Ayan
    Created: 06-02-26
    Example: run compare_kinematic_maps.py --do_all_obj --snr_cut 3
             run compare_kinematic_maps.py --id 8512 --snr_cut 0 --plot_snr
             run compare_kinematic_maps.py --id 8512 --snr_cut 3
'''

from header import *
from util import *
from make_msa3d_line_maps import read_line_maps_fits, read_msa3d_catalog, plot_2D_map

start_time = datetime.now()

# ---------------------------------------------------------------------------------------------------
def read_mengting_velocity_map(filename):
    '''
    Reads in the fits file containing kinematic maps prepared by Mengting
    Returns dictionary with numpy arrays: velocity and velocity disperson maps
    '''
    print(f'Reading in {filename}..')
    hdul = fits.open(filename)
    vel_map = hdul['VELOCITY HA'].data
    sigma_map = hdul['INTRINSIC SIGMA HA'].data
    kin_dict = {'vel': vel_map, 'sigma': sigma_map}

    return kin_dict

# ---------------------------------------------------------------------------------------------------
def read_takafumi_velocity_map(filename):
    '''
    Reads in the fits file containing kinematic maps prepared by Takafumi
    Returns dictionary with numpy arrays: velocity and velocity disperson maps, along with uncertainties
    '''
    print(f'Reading in {filename}..')
    hdul = fits.open(filename)
    mask = hdul['MASK'].data

    vel_map = hdul['MEAN_KMS'].data
    vel_map_err = hdul['MEAN_KMS_ERR'].data

    sigma_map = hdul['SIGMA_KMS'].data
    sigma_map_err = hdul['SIGMA_KMS_ERR'].data
    
    vel_map = np.where(mask, np.nan, vel_map)
    vel_map_err = np.where(mask, np.nan, vel_map_err)
    sigma_map = np.where(mask, np.nan, sigma_map)
    sigma_map_err = np.where(mask, np.nan, sigma_map_err)

    kin_dict = {'vel': vel_map, 'sigma': sigma_map, 'vel_err': vel_map_err, 'sigma_err': sigma_map_err}

    return kin_dict


# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    args.fontfactor = 1.2
    args.id_arr = args.id

    label_dict = {'vel':r'v$_{los}$', 'sigma': r'$\sigma$_v', 'snr':'SNR'}
    unit_dict = {'vel': ' (km/s)', 'sigma': ' km/s', 'snr':''}
    cmin_dict = {'vel': -100, 'sigma': 0, 'snr': 0}
    cmax_dict = {'vel': 100, 'sigma': 200, 'snr': 20}
    cmap_dict = {'vel':'PuOr', 'sigma': 'plasma', 'snr': 'viridis'}
    
    # -------------setup directories and global variables----------------
    maps_fits_dir = args.output_dir / 'maps'
    quants_fits_dir = args.output_dir / 'quants'
    args.fig_dir = args.output_dir / 'plots'
    mengting_fits_dir = args.input_dir / 'Mengting_velocity_maps'
    takafumi_fits_dir = args.input_dir / 'Takafumi_velocity_maps'

    catalog_file = args.input_dir / 'redshifts.dat'
    tie_vdisp_text = '_tie_vdisp' if args.tie_vdisp else ''
    snr_cut_text = f'_snr{args.snr_cut}' if args.snr_cut is not None else ''
    dered_text = f'_nodered' if args.nodered else ''

    # ----------------reading in catalog---------------------
    df = read_msa3d_catalog(catalog_file)

    # ---------matching with available objects from others--------------
    mengting_ids = [int(Path(item).stem[15:-10]) for item in glob.glob(str(mengting_fits_dir) + '/*kinematic.fits')]
    df = df[df['id'].isin(mengting_ids)]
    
    takafumi_ids = [int(Path(item).stem[:-9]) for item in glob.glob(str(takafumi_fits_dir) + '/*velfield.fits')]
    df = df[df['id'].isin(takafumi_ids)].reset_index(drop=True)
    
    print(f'{len(df)} objects remain overlapping between Mengting and Takafumis objects')
    if not args.do_all_obj: df = df[df['id'].isin(args.id_arr)].reset_index(drop=True)

    # ----------------looping over the objects in this chunk-------------
    for index, obj in df.iterrows():
        print(f'Commencing ({index + 1}/{len(df)}) ID {obj["id"]}..')
        args.id =obj['id']
        args.z = obj['redshift']
        args.distance = cosmo.luminosity_distance(args.z)
        args.upto_arcsec = args.upto_kpc * cosmo.arcsec_per_kpc_proper(args.z).value # arcsec
        args.upto_pix = args.upto_arcsec / msa_pix_size_arcsec

        # ------determining directories and filenames---------
        maps_fits_file = maps_fits_dir / f'{args.id:05d}{tie_vdisp_text}.maps.fits'
        mengting_fits_file = mengting_fits_dir / f'cube_jan-22-24_{args.id}_kinematic.fits'
        takafumi_fits_file = takafumi_fits_dir / f'{args.id}_velfield.fits'

        # -----------read kinematic maps--------------
        fit_results_ayan, _ = read_line_maps_fits(maps_fits_file, args)
        fit_results_mengting = read_mengting_velocity_map(mengting_fits_file)
        fit_results_takafumi = read_takafumi_velocity_map(takafumi_fits_file)

        line = 'H-alpha'
        params = ['vel', 'sigma']

        # -------------setup the figure-----------------
        nrow, ncol = 2, 3
        fig, axes = plt.subplots(nrow, ncol, figsize=(7, 5.5))
        fig.subplots_adjust(left=0.08, right=0.93, top=0.93, bottom=0.1, wspace=0., hspace=0.1)

        if args.plot_snr:
            fig_snr, axes_snr = plt.subplots(nrow, ncol - 1, figsize=(8, 5))
            fig_snr.subplots_adjust(left=0.1, right=0.93, top=0.93, bottom=0.12, wspace=0., hspace=0.1)

        # ----------looping through fitted parameters--------------------
        for index, param in enumerate(['vel', 'sigma']):
            map_ayan = fit_results_ayan[line][param]
            map_mengting = fit_results_mengting[param]
            map_takafumi = fit_results_takafumi[param]

            axes[index][0] = plot_2D_map(map_ayan, axes[index][0], f'{param} (Ayan)', args, cmap=cmap_dict[param], takelog=False, vmin=cmin_dict[param], vmax=cmax_dict[param], hide_xaxis=index < nrow - 1, hide_yaxis=False, hide_cbar=True)            
            axes[index][1] = plot_2D_map(map_takafumi, axes[index][1], f'{param} (Takafumi)', args, cmap=cmap_dict[param], takelog=False, vmin=cmin_dict[param], vmax=cmax_dict[param], hide_xaxis=index < nrow - 1, hide_yaxis=True, hide_cbar=True)
            axes[index][2] = plot_2D_map(map_mengting, axes[index][2], f'{param} (Mengting)', args, cmap=cmap_dict[param], takelog=False, vmin=cmin_dict[param], vmax=cmax_dict[param], hide_xaxis=index < nrow - 1, hide_yaxis=True, hide_cbar=False)

            if args.plot_snr:
                snr_ayan = map_ayan / fit_results_ayan[line][f'{param}_err']
                snr_takafumi = map_takafumi / fit_results_takafumi[f'{param}_err']

                axes_snr[index][0] = plot_2D_map(snr_ayan, axes_snr[index][0], f'{param} SNR (Ayan)', args, cmap=cmap_dict['snr'], takelog=False, vmin=cmin_dict['snr'], vmax=cmax_dict['snr'], hide_xaxis=index < nrow - 1, hide_yaxis=False, hide_cbar=True)            
                axes_snr[index][1] = plot_2D_map(snr_takafumi, axes_snr[index][1], f'{param} SNR (Takafumi)', args, cmap=cmap_dict['snr'], takelog=False, vmin=cmin_dict['snr'], vmax=cmax_dict['snr'], hide_xaxis=index < nrow- 1, hide_yaxis=True, hide_cbar=False)

        fig.text(0.1, 0.98, f'ID {args.id}', fontsize=args.fontsize, c='k', ha='left', va='top')
        if args.plot_snr:
            fig_snr.text(0.1, 0.98, f'ID {args.id}', fontsize=args.fontsize, c='k', ha='left', va='top')

        # ---------saving figures-------------------
        figname = f'{args.id}_compare_velocity_maps{tie_vdisp_text}{snr_cut_text}{dered_text}.png'
        save_fig(fig, args.fig_dir, figname, args)

        if args.plot_snr:
            save_fig(fig_snr, args.fig_dir, figname.replace('velocity', 'velocity_snr'), args)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
