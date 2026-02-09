'''
    Filename: make_msa3d_line_maps.py
    Notes: Performs emission line fitting on MSA-3D datacubes and stores them in fits files
    Author : Ayan
    Created: 06-02-26
    Example: run make_msa3d_line_maps.py --do_all_obj
             run make_msa3d_line_maps.py --id 2145 --plot_line_maps
'''

from header import *
from util import *

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def get_ifu_cube(cube_fits_file, err_fits_file=None):
    '''
    Reads in the IFU datacube from a given filename
    Also reads in the corresponding errorcube, if that filename is provided
    Returns both cubes as 3D numpy arrays as well as the wavelength array and the wcs info
    '''
    print(f'\tReading in IFU cube from {cube_fits_file}..')
    hdu = fits.open(cube_fits_file)[0]
    cube = hdu.data
    if err_fits_file is None: cube_err = None
    else: cube_err = fits.open(err_fits_file)[0].data

    # ----------getting the wavelength array----------
    header = hdu.header
    n_wave = header['NAXIS3']
    wcs = pywcs.WCS(header)
    pixel_coords = np.stack([np.zeros(n_wave), np.zeros(n_wave), np.arange(n_wave)], axis=1)
    wavelengths = wcs.all_pix2world(pixel_coords, 0)[:, 2]
    wavelengths *= 1e10 # from meters to Angstroms

    return cube, cube_err, wavelengths, wcs

# --------------------------------------------------------------------------------------------------------------------
def linefit_cube(cube, wavelengths, args, cube_err=None, flam_col='flam', flam_u_col='flam_u', wave_col='rest_wave', label_col='labels', cont_col='cont'):
    '''
    Runs spaxel-by-spaxel emission line fitting of the emission line list provided (as args.line_list)
    Accounts for the corresponding errorcube, if provided
    Returns N x 2D emission line maps, including corresponding uncertainties
    '''
    img_shape = np.shape(cube)[1:]
    rest_wave_arr = wavelengths / (1 + args.z) # converting to rest-frame wavelength, in Angstrom

    # ---------initialising lines to fit dataframe-----------------
    df_lines = pd.DataFrame({wave_col: [rest_wave_dict[item] for item in args.line_list], label_col: args.line_list})
    df_lines = df_lines[df_lines[wave_col].between(np.min(rest_wave_arr), np.max(rest_wave_arr))]
    print(f'\tFitting cube for ID {args.id} for {len(df_lines)} lines..')

   # ---------initialising fit result dict-----------------
    fit_results = {}
    params = ['flux', 'flux_err', 'vel', 'vel_err', 'sigma', 'sigma_err']
    for label in df_lines[label_col]:
        fit_results[label] = {p: np.full(img_shape, np.nan) for p in params}

    # --------------looping over all spaxels---------------
    start_time3 = datetime.now()
    for i in range(img_shape[0]):
        for j in range(img_shape[1]):
            if not args.silent: print(f'\t\tFitting spectrum in cell ({(i + 1) * (j + 1)}/{img_shape[0] * img_shape[1]}) ({i},{j})..')
            flux = cube[:, i, j]
            if cube_err is None: flux_err = np.zeros(np.shape(flux))
            else: flux_err = cube_err[:, i, j]

            # -------initiliasing and filtering spectrum dataframe------------
            df_spec = pd.DataFrame({wave_col: rest_wave_arr, 
                                    flam_col: flux.astype(np.float64),
                                    flam_u_col: flux_err.astype(np.float64),
                                    })
            check_cols = [flam_col,flam_u_col]
            df_spec = df_spec.dropna(subset=check_cols)
            df_spec = df_spec[np.isfinite(df_spec[check_cols]).all(axis=1)]
            df_spec = df_spec[df_spec[flam_u_col] > 0]

            # ------------calling line fitter for this spaxel---------------
            if len(df_spec) > 0:
                fit_results_spaxel = linefit_spaxel(df_spec, df_lines, args, i=i, j=j, flam_col=flam_col, flam_u_col=flam_u_col, wave_col=wave_col, label_col=label_col, cont_col=cont_col)
                
                # ---------populating the fitted dict------------
                for label in df_lines[label_col]:
                    for param in params:
                        fit_results[label][param][i, j] = fit_results_spaxel[label][param]
            elif not args.silent:
                print(f'\t\t..no valid spectrum in this cell. Moving on..')
            
    print(f'\nCompleted fitting ID {args.id} in {timedelta(seconds=(datetime.now() - start_time3).seconds)}')

    return fit_results

# --------------------------------------------------------------------------------------------------------------------
def gaussian(x, amp, mean, sigma):
    return amp * np.exp(-(x - mean)**2 / (2 * sigma**2))

# --------------------------------------------------------------------------------------------------------------------
def global_gaussian_model(x, *popt, rest_waves=None, tie_vdisp=False):
    '''
    Computes n-component Gaussian model for n emission lines
    Returns np.array
    '''
    model = np.zeros_like(x)

    if tie_vdisp: # Structure: [v_los, shared_sigma, amp1, amp2, ... ampN]
        shared_sigma = popt[1]
        for i, rest_wave in enumerate(rest_waves):
            amp = popt[i + 2]
            mu = rest_wave * (1 + popt[0] / c_km_s)
            sigma_wave = (shared_sigma * mu) / c_km_s
            model +=  gaussian(x, amp, mu, sigma_wave)
    else: # Structure: [v_los, amp1, sig1, amp2, sig2, ... ampN, sigN]
        for i, rest_wave in enumerate(rest_waves):
            amp = popt[2*i + 1]
            sig_kms = popt[2*i + 2]
            mu = rest_wave * (1 + popt[0] / c_km_s)
            sigma_wave = (sig_kms * mu) / c_km_s
            model += gaussian(x, amp, mu, sigma_wave)
            
    return model

# --------------------------------------------------------------------------------------------------------------------
def linefit_spaxel(df_spec, df_lines, args, i=None, j=None, flam_col='flam', flam_u_col='flam_u', wave_col='rest_wave', label_col='labels', cont_col='cont'):
    '''
    Runs emission line fitting on one spectrum (spaxel) of the emission line wavelengths provided (df_lines)
    Accounts for the corresponding error spectrum, if provided
    Returns N emission line fluxes, including corresponding uncertainties
    '''
    if i is None: i = 'X'
    if j is None: j = 'X'
    # ---------initialising fit result dict-----------------
    fit_results_spaxel = {}
    params = ['flux', 'flux_err', 'vel', 'vel_err', 'sigma', 'sigma_err']
    for label in df_lines[label_col]:
        fit_results_spaxel[label] = {p: np.full(len(df_lines), np.nan) for p in params}

    # -----------continuum fitting (with masking)-----------
    cont_mask = np.ones(len(df_spec), dtype=bool)
    for line_wave in df_lines[wave_col]:
        cont_mask &= (np.abs(df_spec[wave_col].values - line_wave) > args.mask_window)
    df_spec_masked = df_spec[cont_mask]
    
    spline = UnivariateSpline(df_spec_masked[wave_col], df_spec_masked[flam_col], w=1/df_spec_masked[flam_u_col], s=len(df_spec_masked))
    df_spec[cont_col] = spline(df_spec[wave_col])
    contsub_col = 'flam_contsub'
    df_spec[contsub_col] = df_spec[flam_col] - df_spec[cont_col] # Subtract continuum

    # -----------initialising the parameters-----------------------
    amp_guess = np.max(df_spec[contsub_col])
    sig_guess = 2.0 # Roughly 100-200 km/s depending on resolution
    vel_guess = 100. # km/s

    if args.tie_vdisp: # Structure: [v_los, shared_sigma, amp1, amp2, ... ampN]
        p0 = [vel_guess, sig_guess] + [amp_guess] * len(df_lines)
        lbounds = [-500, 0.5] + [0] * len(df_lines)
        ubounds = [500, 15] + [np.inf] * len(df_lines)
    
    else: # Structure: [v_los, amp1, sig1, amp2, sig2, ... ampN, sigN]
        p0 = [vel_guess] + [amp_guess, sig_guess] * len(df_lines)
        lbounds = [-500] + [0, 0.5] * len(df_lines)
        ubounds = [500] + [np.inf, 15] * len(df_lines)
        
    # -----------------fitting the spectrum---------------------
    model_func = lambda x, *params: global_gaussian_model(x, *params, rest_waves=df_lines['rest_wave'].values, tie_vdisp=args.tie_vdisp)
    try:
        popt, pcov = curve_fit(model_func, df_spec[wave_col], df_spec[contsub_col], p0=p0, sigma=df_spec[flam_u_col], bounds=(lbounds, ubounds))
        perr = np.sqrt(np.diag(pcov))
        
        v_fit, v_err = popt[0], perr[0]

        for ind in range(len(df_lines)):
            if args.tie_vdisp: # Structure: [v_los, shared_sigma, amp1, amp2, ... ampN]
                sig_fit, sig_err = popt[1], perr[1]
                amp_fit, amp_err = popt[ind + 2], perr[ind + 2]
            else: # Structure: [v_los, amp1, sig1, amp2, sig2, ... ampN, sigN]
                amp_fit, amp_err = popt[2 * ind + 1], perr[2 * ind + 1]
                sig_fit, sig_err = popt[2 * ind + 2], perr[2 * ind + 2]
            
            # Integrated Flux Calculation
            mu_w = df_lines.iloc[ind][wave_col] * (1 + ufloat(v_fit, v_err) / c_km_s)
            sig_w = (ufloat(sig_fit, sig_err) * mu_w) / c_km_s
            flux = ufloat(amp_fit, amp_err) * sig_w * np.sqrt(2 * np.pi)
            
            fit_results_spaxel[df_lines.iloc[ind][label_col]] = {
                'flux': flux.n, 'flux_err': flux.s,
                'vel': v_fit, 'vel_err': v_err,
                'sigma': sig_fit, 'sigma_err': sig_err
            }

        # --------------plotting the spectrum------------------------
        fig, ax = plt.subplots(1, figsize=(12, 3), layout='constrained')

        ax = plot_line_fit_spaxel(ax, df_spec, df_lines, args, popt=popt, flam_col=flam_col, flam_u_col=flam_u_col, wave_col=wave_col, label_col=label_col, cont_col=cont_col)
        ax = plot_linelist(ax, df_lines, fontsize=args.fontsize / args.fontfactor, color='cornflowerblue')

        figname = f'{args.id}_pixel_{i}-{j}_linefit.png'
        save_fig(fig, args.linefit_fig_dir, figname, args, silent=True)
        if args.debug_linefit: sys.exit(f'Stopping at pixel ({i},{j}) since --debug_linefit is turned on..')
        else: plt.close()
    
    except Exception as e:
        print(f'Fit failed for ({i},{j}) due to: {e}')

    return fit_results_spaxel

# --------------------------------------------------------------------------------------------------------------------
def plot_line_fit_spaxel(ax, df_spec, df_lines, args, popt=None, color='orangered', flam_col='flam', flam_u_col='flam_u', wave_col='rest_wave', label_col='labels', cont_col='cont'):
    '''
    Plots the spectrum along one spaxel in the given axis handle, for debugging purposes
    Returns the axis handle
    '''
    norm_factor = 1 #if args.show_log_flux else 1e-19

    # ---------------plot the observed spectrum---------------
    ax.step(df_spec[wave_col], df_spec[flam_col] / norm_factor, lw=1, c=color, alpha=1, where='mid')
    ax.plot(df_spec[wave_col], df_spec[cont_col] / norm_factor, lw=1, c='grey', alpha=1)
    ax.fill_between(df_spec[wave_col], (df_spec[flam_col] - df_spec[flam_u_col]/2) / norm_factor, (df_spec[flam_col] + df_spec[flam_u_col]/2) / norm_factor, lw=0, color=color, alpha=0.5, step='pre')#, drawstyle='steps')

    # -----------plot the fitted spectrum--------------------
    if popt is not None:
        x_arr = np.linspace(df_spec[wave_col].min(), df_spec[wave_col].max(), 1000)
        total_model = global_gaussian_model(x_arr, *popt, rest_waves=df_lines['rest_wave'].values, tie_vdisp=args.tie_vdisp)
        ax.plot(x_arr, total_model / norm_factor, color='k', lw=1.5, linestyle='--', label='Fitted model', zorder=5)

    if args.show_log_flux: ax.set_yscale('log')

    # ----------------annotate axis labels---------------
    ax.set_xlabel(r'Rest-frame wavelength ($\AA$)', fontsize=args.fontsize)
    #ylabel = r'f$_{\lambda}$ ' + '(%.0e ' % norm_factor + r'ergs/s/cm$^2$/A)'
    ax.set_ylabel(r'f$_{\lambda}$ ergs/s/cm$^2$/A)', fontsize=args.fontsize)
    #if not args.show_log_flux: ax.set_ylim(0, args.flam_max) # flam_max should be in units of 1e-19 ergs/s/cm^2/A
    ax.tick_params(axis='both', which='major', labelsize=args.fontsize)

    # ---observed wavelength axis-------
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=5, prune='both'))
    ax2.set_xticklabels(['%.2F' % (item * (1 + args.z) / 1e4) for item in ax2.get_xticks()], fontsize=args.fontsize)
    ax2.set_xlabel(r'Observed wavelength ($\mu$)', fontsize=args.fontsize)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_linelist(ax, df_lines, fontsize=10, color='cornflowerblue'):
    '''
    Plots a list of emission line wavelengths on the given axis
    Returns axis handle
    '''
    last_flipped = True # flip switch for determining if last label was placed to the left or right of the vertical line

    for index in range(len(df_lines)):
        this_wave = df_lines.iloc[index]['rest_wave']
        ax.axvline(this_wave, c=color, lw=1, alpha=0.5)
        if last_flipped:
            xpos = this_wave - np.diff(ax.get_xlim())[0] * 0.04
        else:
            xpos = this_wave + np.diff(ax.get_xlim())[0] * 0.01
        last_flipped = not last_flipped
        ypos = ax.get_ylim()[1] * 0.98
        ax.text(xpos, ypos, df_lines.iloc[index]['labels'].strip(), rotation=90, va='top', ha='left', fontsize=fontsize)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def save_line_maps_fits(fit_results, maps_fits_file, wcs, args):
    '''
    Saves the N x 2D emission line maps as N-extension fits files, along with other relevant header info
    '''
    # ------------setting up primary header---------------
    w2d = wcs.dropaxis(2)
    spatial_header = w2d.to_header()
    primary_hdu = fits.PrimaryHDU(header=spatial_header)

    primary_hdu.header['CONTENT'] = 'Emission Line Maps'
    primary_hdu.header['ID'] = args.id
    primary_hdu.header['REDSHIFT'] = args.z
    
    hdul = fits.HDUList([primary_hdu])

    params = [
        ('flux', 'FLUX', 'erg/s/cm2'),
        ('flux_err', 'FLUX_ERR', 'erg/s/cm2'),
        ('v_los', 'VEL', 'km/s'),
        ('v_err', 'VEL_ERR', 'km/s'),
        ('sigma', 'SIGMA', 'km/s'),
        ('sigma_err', 'SIGMA_ERR', 'km/s')
    ]

    # ------------looping over fitted lines---------------
    for label in list(fit_results.keys()):   
        line_data = fit_results[label]

        for key, ext_suffix, unit in params:
            data = line_data.get(key)
            if data is None: continue
                
            hdu = fits.ImageHDU(data=data.astype(np.float32), header=spatial_header)
            hdu.header['EXTNAME'] = f'{label}_{ext_suffix}'
            hdu.header['LINE'] = label
            hdu.header['BUNIT'] = unit
            hdul.append(hdu)

    hdul.writeto(maps_fits_file, overwrite=True)
    print(f'Successfully saved {len(hdul)-1} extensions to {maps_fits_file}"')
    return

# --------------------------------------------------------------------------------------------------------------------
def plot_2D_map(image, ax, label, args, cmap='cividis', clabel='', takelog=True, vmin=None, vmax=None, hide_xaxis=False, hide_yaxis=False, hide_cbar=True, hide_cbar_ticks=False, cticks_integer=True):
    '''
    Plots a given 2D image in a given axis
    Returns the axis handle
    '''

    if takelog: image =  np.log10(image.data)

    offset = args.pix_size_arcsec / 2 # half a pixel offset to make sure cells in 2D plot are aligned with centers and not edges
    args.extent = (-args.arcsec_limit - offset, args.arcsec_limit - offset, -args.arcsec_limit - offset, args.arcsec_limit - offset)

    p = ax.imshow(image, cmap=cmap, origin='lower', extent=args.extent, vmin=vmin, vmax=vmax)
    
    ax.scatter(0, 0, marker='x', s=10, c='grey')
    ax.set_aspect('auto') 
    
    ax = annotate_axes(ax, '', '', fontsize=args.fontsize / args.fontfactor, label=label, clabel=clabel, hide_xaxis=hide_xaxis, hide_yaxis=hide_yaxis, hide_cbar=hide_cbar, p=p, hide_cbar_ticks=hide_cbar_ticks, cticks_integer=cticks_integer)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_line_maps(fit_results, args):
    '''
    Plots the N emission line maps in one figure and saves it to file
    Returns the figure handle
    '''
    # --------------setting up figures------------------------
    nrow, ncol = 3, 3
    cmin, cmax, ncbins = -19, -16, 5
    cmin_snr, cmax_snr = 0, 10
    cmap = 'cividis'
    fig, axes = plt.subplots(nrow, ncol, figsize=(10, 10), layout='constrained')

    if args.plot_snr:
        fig_snr, axes_snr = plt.subplots(nrow, ncol, figsize=(10, 10), layout='constrained')

    # ----------plot line maps--------------------
    for index, fitted_line in enumerate(list(fit_results.keys())):
        row = index // ncol
        col = index % ncol
        
        # ------------extract linemap from fit_results dict--------------
        fluxmap = fit_results[fitted_line]['flux']
        axes[row, col] = plot_2D_map(fluxmap, axes[row, col], f'{args.id}: {args.line_list[index]}', args, cmap=cmap, takelog=True, vmin=cmin, vmax=cmax, hide_xaxis=index < row - 1, hide_yaxis=col > 0)

        if args.plot_snr:
            errmap = fit_results[fitted_line]['flux_err']
            snrmap = fluxmap / errmap
            axes_snr[row, col] = plot_2D_map(snrmap, axes_snr[row, col], f'{args.id}: {args.line_list[index]}: SNR', args, cmap=cmap, takelog=False, vmin=cmin_snr, vmax=cmax_snr, hide_xaxis=index < row - 1, hide_yaxis=col > 0)

    # --------common colorbar------------------
    fig = make_colorbar_top(fig, axes, 'Flux', cmap, cmin, cmax, ncbins, args.fontsize, aspect=60)
    if args.plot_snr: fig_snr = make_colorbar_top(fig_snr, axes_snr, 'SNR', cmap, cmin_snr, cmax_snr, ncbins, args.fontsize, aspect=60)

    # ---------saving figures-------------------
    figname = f'{args.id}_fluxmaps.png'
    save_fig(fig, args.fig_dir, figname, args)    
    if args.plot_snr: save_fig(fig_snr, args.fig_dir, Path(str(figname).replace('flux', 'snr')), args)

    return fig

# --------------------------------------------------------------------------------------------------------------------
def read_line_maps_fits(filename):
    '''
    Reads a multi-extension FITS file and reconstructs the fit_results dictionary.
    Returns:
    fit_results: dict. Format {label: {param: 2D array}}
    spatial_header: astropy.io.fits.Header. The WCS/spatial info
    '''
    fit_results = {}
    
    with fits.open(filename) as hdul:
        spatial_header = hdul[0].header
        suffix_to_key = {
            'FLUX': 'flux',
            'FLUX_ERR': 'flux_err',
            'VEL': 'v_los',
            'VEL_ERR': 'v_err',
            'SIGMA': 'sigma',
            'SIGMA_ERR': 'sigma_err'
        }
        
        # Loop through all extensions starting from index 1
        for i in range(1, len(hdul)):
            extname = hdul[i].name
            
            if '_' not in extname:
                continue
                
            # Split the EXTNAME into the line label and the parameter suffix
            parts = extname.rsplit('_', 1)
            label = parts[0]
            suffix = parts[1]
            
            if suffix not in suffix_to_key:
                continue
            
            if label not in fit_results:
                fit_results[label] = {}
            
            dict_key = suffix_to_key[suffix]
            fit_results[label][dict_key] = hdul[i].data.copy()

    return fit_results, spatial_header

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
    args.fig_dir = args.output_dir / 'plots'

    catalog_file = args.input_dir / 'redshifts.dat'

    # ----------------reading in catalog---------------------
    columns = ['id', 'ra', 'dec', 'type', 'redshift', 'mass', 'sfr', 'ssfr', 'Zgrad_kpc', 'Zgrad_kpc_u', 'Zgrad_re', 'Zgrad_re_u', 'num']
    df = pd.read_csv(catalog_file, delim_whitespace=True, comment='#', names=columns)
    ndf = len(df)

    # ----------curtail dataframe to relevant redshifts------------
    msa_3d_wave_lim = [970, 1820] # in nm
    hbeta_zlim = get_zrange_for_line('H-beta', obs_wave_range=msa_3d_wave_lim)
    halpha_zlim = get_zrange_for_line('H-alpha', obs_wave_range=msa_3d_wave_lim)
    zlim = (hbeta_zlim[0], halpha_zlim[1])
    df = df[df['redshift'].between(zlim[0], zlim[1])]
    print(f'{len(df)} of {ndf} galaxies remain after redshift cut of {zlim[0]:.2f} < z < {zlim[1]:.2f}')

    if not args.do_all_obj: df = df[df['id'].isin(args.id_arr)]

    # ----------------looping over the objects in this chunk-------------
    for index, obj in df.iterrows():
        start_time2 = datetime.now()
        print(f'Commencing ({index + 1}/{len(df)}) ID {obj["id"]}..')
        args.id =obj['id']
        args.z = obj['redshift']

        # ------determining directories and filenames---------
        args.cube_fits_file = cube_fits_dir / f'cube_square_medians_{args.id:d}_hdr_rect.fits'
        args.err_fits_file = err_fits_dir / f'cube_square_medians_{args.id:d}_hdr_rect_err.fits'
        args.maps_fits_file = maps_fits_dir / f'{args.id:05d}.maps.fits'

        args.linefit_fig_dir = args.fig_dir / f'{args.id}_linefit_plots'
        args.linefit_fig_dir.mkdir(exist_ok=True, parents=True)
        
        if not os.path.exists(args.maps_fits_file) or args.clobber:
            # -----------read in the cube--------------
            cube, cube_err, wavelengths, wcs = get_ifu_cube(args.cube_fits_file, err_fits_file=args.err_fits_file)

            # -----------emission line fitting--------------
            fit_results = linefit_cube(cube, wavelengths, args, cube_err=cube_err)

            # -----------save the emission maps in fits file-------------
            save_line_maps_fits(fit_results, args.maps_fits_file, wcs, args)
        else:
            fit_results = read_line_maps_fits(args)
        
        # --------plot the emission line maps-------------
        if args.plot_line_maps: fig = plot_line_maps(fit_results, args)

        print(f'\nCompleted ID {args.id} in {timedelta(seconds=(datetime.now() - start_time2).seconds)}, {len(df) - index - 1} to go!')

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
