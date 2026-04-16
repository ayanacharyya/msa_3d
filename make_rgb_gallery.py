'''
    Filename: make_rgb_gallery.py
    Notes: Make a gallery of RGB images based on the MSA-3D emission line maps
    Author : Ayan
    Created: 30-03-26
    Example: run make_rgb_gallery.py --snr_cut 3
'''

from header import *
from util import *
setup_plot_style()
from make_msa3d_line_maps import read_line_maps_fits, read_msa3d_catalog, compute_rgb

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def plot_rgb(image, ax, label, args, hide_xaxis=False, hide_yaxis=False):
    '''
    Plots a given 2D image in a given axis
    Returns the axis handle
    '''
    ax.set_facecolor('k')
    extent = (-args.upto_arcsec, args.upto_arcsec, -args.upto_arcsec, args.upto_arcsec)
    p = ax.imshow(image, origin='lower', extent=extent)
        
    ax.scatter(0, 0, marker='x', s=10, c='grey')
    ax.set_aspect('auto') 
    
    ax = annotate_axes(ax, 'arcsec', 'arcsec', fontsize=args.fontsize / args.fontfactor, label=label, hide_xaxis=hide_xaxis, hide_yaxis=hide_yaxis, label_color='w', bbox=False)

    return ax

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    args.fontfactor = 1.4
    args.id_arr = args.id
    args.upto_arcsec = 0.8
    
    # -------------setup directories and global variables----------------
    maps_fits_dir = args.output_dir / 'maps'
    args.fig_dir = args.output_dir / 'plots'
    
    catalog_file = args.input_dir / 'redshifts.dat'
    tie_vdisp_text = '_tie_vdisp' if args.tie_vdisp else ''
    snr_cut_text = f'_snr{args.snr_cut}'

    # ----------------reading in catalog---------------------
    df = read_msa3d_catalog(catalog_file)

    # -----------setting up figure--------------------
    nrow, ncol = 6, 6
    fig, axes = plt.subplots(nrow, ncol, figsize=(8.5, 7.5))
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.06, top=0.98, wspace=0., hspace=0.)

    # ----------------looping over the objects in this chunk-------------
    for index, obj in df.iterrows():
        print(f'Commencing ({index + 1}/{len(df)}) ID {obj["id"]}..')
        row = index // ncol
        col = index % ncol
        
        args.id =obj['id']
        args.z = obj['redshift']
        args.upto_pix = args.upto_arcsec / msa_pix_size_arcsec

        # ------determining directories and filenames---------
        args.maps_fits_file = maps_fits_dir / f'{args.id:05d}{tie_vdisp_text}.maps.fits'

        # ---------reading in existing line fitted file-----------
        if not os.path.exists(args.maps_fits_file):
            print(f'Line fitted file for {args.id} does not exist')
            axes[row, col].remove()
            continue
        fit_results, _ = read_line_maps_fits(args.maps_fits_file, args)

        # -----------gettign RGB image------------------
        rgb_image = compute_rgb(fit_results, args, rlines='SII-6717,SII-6730', glines='H-alpha', blines='OIII-5007')

        if type(rgb_image) == str:
            axes[row, col] = plot_rgb(np.zeros((msa_npix_x, msa_npix_y)) * np.nan, axes[row, col], f'{args.id}', args, hide_xaxis=row < nrow - 1, hide_yaxis=col > 0)
            axes[row, col].text(0.5, 0.5, f'{rgb_image}\nunavailable', color='w', ha='center', va='center', transform=axes[row, col].transAxes, fontsize=args.fontsize / args.fontfactor)
        else:
            rgb_image = cut_2Dmap(rgb_image, args.upto_pix)
            axes[row, col] = plot_rgb(rgb_image, axes[row, col], f'{args.id}', args, hide_xaxis=row < nrow - 1, hide_yaxis=col > 0)

        for spine in axes[row, col].spines.values(): spine.set_edgecolor('white')
        axes[row, col].tick_params(color='white', which='both')
    
    # ----------to remove empty subplots--------------------
    for ind in range(index + 1, nrow * ncol):
        row = ind // ncol
        col = ind % ncol
        axes[row, col].remove()

    # ----------save figure---------------
    figname = f'all_RGB_{tie_vdisp_text}{snr_cut_text}.png'
    save_fig(fig, args.fig_dir, figname, args)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
