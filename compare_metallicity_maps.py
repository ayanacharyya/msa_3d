'''
    Filename: compare_metallicity_maps.py
    Notes: Compares metallicity maps from Ayan's emission line fitting vs others'
    Author : Ayan
    Created: 06-02-26
    Example: run compare_metallicity_maps.py --do_all_obj --Zdiag NB --snr_cut 3
             run compare_metallicity_maps.py --id 8512 --Zdiag NB --snr_cut 0 --plot_snr
             run compare_metallicity_maps.py --id 8512 --Zdiag R3 --snr_cut 3
             run compare_metallicity_maps.py --do_all_obj --Zdiag N2 --use_P04 --snr_cut 0
'''

from header import *
from util import *
from make_msa3d_line_maps import read_msa3d_catalog, plot_2D_map
from make_metallicity_sfr_maps import read_quant_maps_fits, get_quant_from_quant_maps

start_time = datetime.now()

# ---------------------------------------------------------------------------------------------------
def read_mengting_metallicity_map(filename):
    '''
    Reads in the fits file containing metallicity maps prepared by Mengting
    Returns 2D numpy arrays: metallicity map and corresponding uncertainty
    '''
    print(f'Reading in {filename}..')
    label = 'METAL'
    hdul = fits.open(filename)
    logOH_map = hdul[f'{label}'].data
    logOH_map_err = hdul[f'{label}_ERR'].data

    return logOH_map, logOH_map_err

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    args.fontfactor = 1.2
    args.id_arr = args.id

    logOH_min, logOH_max, logOH_cmap = 7.0, 9.0, 'viridis'
    snr_min, snr_max, snr_cmap = 0, 20, 'rainbow'

    # -------------setup directories and global variables----------------
    quants_fits_dir = args.output_dir / 'quants'
    args.fig_dir = args.output_dir / 'plots'
    mengting_fits_dir = args.input_dir / 'Mengting_metallicity_maps'

    catalog_file = args.input_dir / 'redshifts.dat'
    tie_vdisp_text = '_tie_vdisp' if args.tie_vdisp else ''
    snr_cut_text = f'_snr{args.snr_cut}' if args.snr_cut is not None else ''
    dered_text = f'_nodered' if args.nodered else ''
    P04_text = '_P04' if args.use_P04 else ''
    
    # ----------------reading in catalog---------------------
    df = read_msa3d_catalog(catalog_file)

    # ---------matching with available objects from others--------------
    mengting_ids = [int(Path(item).stem[15:-6]) for item in glob.glob(str(mengting_fits_dir) + '/*metal.fits')]
    df = df[df['id'].isin(mengting_ids)]
    print(f'{len(df)} objects remain overlapping between Mengtings objects')
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
        quants_fits_file = quants_fits_dir / f'{args.id:05d}_Zdiag_{args.Zdiag}{tie_vdisp_text}{snr_cut_text}{dered_text}{P04_text}.quants.fits'
        mengting_fits_file = mengting_fits_dir / f'cube_jan-22-24_{args.id}_metal.fits'

        # -----------read metallicity map--------------
        try:
            quant_maps, spatial_header = read_quant_maps_fits(quants_fits_file)
            logOH_map_ayan, logOH_map_err_ayan, _, logOH_mask_ayan = get_quant_from_quant_maps('logOH', quant_maps)
            logOH_map_ayan = np.where(logOH_mask_ayan, np.nan, logOH_map_ayan)
            logOH_map_err_ayan = np.where(logOH_mask_ayan, np.nan, logOH_map_err_ayan)

            logOH_map_mengting, logOH_map_err_mengting = read_mengting_metallicity_map(mengting_fits_file)
        except Exception as e:
            print(f'Skipping id {args.id} due to {e}')
            continue
        
        # -------------setup the figure-----------------
        fig, axes = plt.subplots(2 if args.plot_snr else 1, 2, figsize=(6.5, 6 if args.plot_snr else 4.5))
        fig.subplots_adjust(left=0.1, right=0.93, top=0.93, bottom=0.12, wspace=0.)
        axes = np.atleast_2d(axes)

        # ----------plotting the metallicity maps--------------------
        axes[0][0] = plot_2D_map(logOH_map_ayan, axes[0][0], f'log O/H + 12 ({args.Zdiag})', args, cmap=logOH_cmap, takelog=False, vmin=logOH_min, vmax=logOH_max, hide_xaxis=args.plot_snr, hide_yaxis=False, hide_cbar=True)
        axes[0][1] = plot_2D_map(logOH_map_mengting, axes[0][1], f'log O/H + 12 (Mengintg)', args, cmap=logOH_cmap, takelog=False, vmin=logOH_min, vmax=logOH_max, hide_xaxis=args.plot_snr, hide_yaxis=True, hide_cbar=False)

        if args.plot_snr:
            snr_map_ayan = logOH_map_ayan / logOH_map_err_ayan
            snr_map_mengting = logOH_map_mengting / logOH_map_err_mengting

            axes[1][0] = plot_2D_map(snr_map_ayan, axes[1][0], f'SNR ({args.Zdiag})', args, cmap=snr_cmap, takelog=False, vmin=snr_min, vmax=snr_max, hide_xaxis=False, hide_yaxis=False, hide_cbar=True)
            axes[1][1] = plot_2D_map(snr_map_mengting, axes[1][1], f'SNR (Mengintg)', args, cmap=snr_cmap, takelog=False, vmin=snr_min, vmax=snr_max, hide_xaxis=False, hide_yaxis=True, hide_cbar=False)

        fig.text(0.1, 0.98, f'ID {args.id}', fontsize=args.fontsize, c='k', ha='left', va='top')

        # ---------saving figures-------------------
        figname = f'{args.id}_compare_metallicity_maps_Zdiag_{args.Zdiag}{tie_vdisp_text}{snr_cut_text}{dered_text}.png'
        save_fig(fig, args.fig_dir, figname, args)    

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
