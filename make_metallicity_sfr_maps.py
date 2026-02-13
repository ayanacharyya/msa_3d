'''
    Filename: make_metallicity_sfr_maps.py
    Notes: Performs metallicity and SFR maps from MSA-3D emission line maps and stores them in fits files
    Author : Ayan
    Created: 06-02-26
    Example: run make_metallicity_sfr_maps.py --do_all_obj --Zdiag NB --plot_met_sfr
             run make_metallicity_sfr_maps.py --id 2145 --Zdiag NB --plot_metallicity
             run make_metallicity_sfr_maps.py --id 2145 --Zdiag N2 --plot_met_sfr
             run make_metallicity_sfr_maps.py --id 2145 --Zdiag R3 --Zbranch low --use_C25 --plot_met_sfr
'''

from header import *
from util import *
from make_msa3d_line_maps import read_line_maps_fits, read_msa3d_catalog, plot_2D_map, get_emission_line_map

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def compute_EB_V(Ha_flux, Hb_flux,verbose=False):
    '''
    Calculates and returns the color excess given observed H alpha and H beta fluxes
    Based on Eqn 4 of Dominguez+2013 (https://iopscience.iop.org/article/10.1088/0004-637X/763/2/145/pdf)
    '''
    # -----handling masks separately because uncertainty package cannot handle masks---------
    if np.ma.isMaskedArray(Ha_flux):
        net_mask = Ha_flux.mask | Hb_flux.mask
        Ha_flux = Ha_flux.data
        Hb_flux = Hb_flux.data
    else:
        net_mask = False

    theoretical_ratio = 2.86

    if hasattr(Hb_flux, "__len__"): # if it is an array
        # --------computing the ratio and appropriate errors------------
        new_mask = unp.nominal_values(Hb_flux) == 0
        Hb_flux[new_mask] = -1 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        obs_ratio = Ha_flux / Hb_flux
        obs_ratio = np.ma.masked_where(new_mask, obs_ratio)

        # --------computing the log of the ratio and appropriate errors------------
        new_mask = unp.nominal_values(obs_ratio.data) <= 0
        obs_ratio[new_mask] = 1e-9 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        EB_V = 1.97 * unp.log10(obs_ratio.data / theoretical_ratio)
        EB_V[obs_ratio.data < theoretical_ratio] = 0
        EB_V = np.ma.masked_where(new_mask | obs_ratio.mask | net_mask, EB_V)

    else: # if it is scalar
        try:
            obs_ratio = Ha_flux / Hb_flux
            if obs_ratio < theoretical_ratio:
                EB_V = ufloat(0, 0)
                if verbose: print(f'Based on integrated fluxes Ha = {Ha_flux:.2e}, Hb = {Hb_flux:.2e}, the observed ratio is LOWER than theoretical ratio, so E(B-V) would be unphysical, so just assuming {EB_V}')
            else:
                EB_V = 1.97 * unp.log10(obs_ratio / theoretical_ratio)
                if verbose: print(f'Based on integrated fluxes Ha = {Ha_flux:.2e}, Hb = {Hb_flux:.2e}, determined E(B-V) =', EB_V)
        except:
            EB_V = ufloat(np.nan, np.nan)

    return EB_V

def get_EB_V(fit_results, args):
    '''
    Computes and returns the spatially resolved as well as integrated dust extinction map from a given dictionary of line flux maps
    Based on Eqn 4 of Dominguez+2013 (https://iopscience.iop.org/article/10.1088/0004-637X/763/2/145/pdf)
    '''

    Ha_map, Ha_int = get_emission_line_map('H-alpha', fit_results, args, dered=False) # do not need to deredden the lines when we are fetching the flux in order to compute reddening
    Hb_map, Hb_int = get_emission_line_map('H-beta', fit_results, args, dered=False)

    EB_V_map = compute_EB_V(Ha_map, Hb_map)
    EB_V_int = compute_EB_V(Ha_int, Hb_int)

    return EB_V_map, EB_V_int
# --------------------------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------------------------
def compute_Z_P25_SFR(O3S2, N2S2):
    '''
    Calculates and returns the SFR metallicity given observed O3S2 and N2S2 line ratios
    Conversion factor is from Eq 1 of Peluso+2025
    '''
    return 8.78 + 0.97 * N2S2 - 0.11 * O3S2 - 0.39 * N2S2 ** 2 + 0.09 * N2S2 * O3S2 + 0.02 * O3S2 ** 2

def compute_Z_P25_AGN(O3S2, N2S2):
    '''
    Calculates and returns the AGN metallicity given observed O3S2 and N2S2 line ratios
    Conversion factor is from Eq 2 of Peluso+2025
    '''
    return 8.85 + 1.10 * N2S2 - 0.04 * O3S2

def compute_Z_P25_Comp(O3S2, N2S2):
    '''
    Calculates and returns the Composite region metallicity given observed O3S2 and N2S2 line ratios
    Conversion factor is from Eq 3 of Peluso+2025
    '''
    return 8.83 + 1.07 * N2S2 + 0.10 * O3S2

def compute_Z_P25(OIII5007_flux, NII6584_flux, SII6717_flux, AGN_map, args):
    '''
    Calculates and returns the SF/AGN metallicity given observed line fluxes
    Conversion factor is from Peluso+2025
    '''
    # -----handling masks separately because uncertainty package cannot handle masks---------
    if np.ma.isMaskedArray(OIII5007_flux):
        net_mask = SII6717_flux.mask | OIII5007_flux.mask | NII6584_flux.mask
        OIII5007_flux = OIII5007_flux.data
        SII6717_flux = SII6717_flux.data
        NII6584_flux = NII6584_flux.data
    else:
        net_mask = False

    if hasattr(SII6717_flux, "__len__"): # if it is an array
        # --------computing the ratio and appropriate errors------------
        O3S2 = take_safe_log_ratio(OIII5007_flux, SII6717_flux)
        N2S2 = take_safe_log_ratio(NII6584_flux, SII6717_flux)

        # --------computing the polynomial and appropriate errors------------
        log_OH = []
        if args.debug_Zdiag:
            fig, ax = plt.subplots(2, 1, figsize=(6, 8), sharex=True)
            ax[0].set_ylabel('O3S2')
            ax[1].set_ylabel('N2S2')
            ax[1].set_xlabel('log(O/H)+12')
            ax[0].set_xlim(7.5, 9.5)

        O3S2_flat = O3S2.data.flatten()
        N2S2_flat = N2S2.data.flatten()
        AGN_map_flat = AGN_map.flatten()
        for index in range(len(O3S2_flat)):
            try:
                if AGN_map_flat[index] > 0: this_log_OH = compute_Z_P25_AGN(O3S2_flat[index], N2S2_flat[index])
                else: this_log_OH = compute_Z_P25_SFR(O3S2_flat[index], N2S2_flat[index])
                log_OH.append(this_log_OH)
                if args.debug_Zdiag:
                    ax[0].scatter(unp.nominal_values(this_log_OH), unp.nominal_values(O3S2_flat[index]), lw=0, s=50, c='r' if AGN_map_flat[index] > 0 else 'b')
                    ax[0].errorbar(unp.nominal_values(this_log_OH), unp.nominal_values(O3S2_flat[index]), xerr=unp.std_devs(this_log_OH), yerr=unp.std_devs(O3S2_flat[index]), fmt='none', lw=0.5, c='r' if AGN_map_flat[index] > 0 else 'b')
                    ax[1].scatter(unp.nominal_values(this_log_OH), unp.nominal_values(N2S2_flat[index]), lw=0, s=50, c='r' if AGN_map_flat[index] > 0 else 'b')
                    ax[1].errorbar(unp.nominal_values(this_log_OH), unp.nominal_values(N2S2_flat[index]), xerr=unp.std_devs(this_log_OH), yerr=unp.std_devs(N2S2_flat[index]), fmt='none', lw=0.5, c='r' if AGN_map_flat[index] > 0 else 'b')
            except:
                log_OH.append(ufloat(np.nan, np.nan))
        log_OH = np.ma.masked_where(O3S2.mask | N2S2.mask, np.reshape(log_OH, np.shape(O3S2)))

    else: # if it is scalar
        try:
            O3S2 = unp.log10(OIII5007_flux / SII6717_flux)
            N2S2 = unp.log10(NII6584_flux / SII6717_flux)
            if AGN_map > 0: log_OH = compute_Z_P25_AGN(O3S2, N2S2)
            else: log_OH = compute_Z_P25_SFR(O3S2, N2S2)
        except:
            log_OH = ufloat(np.nan, np.nan)

    return log_OH

def get_Z_P25(fit_results, args):
    '''
    Computes and returns the spatially resolved as well as intregrated Peluso+2025 metallicity from a given line flux dictionary
    '''
    SII6717_map, SII6717_int = get_emission_line_map('SII-6717', fit_results, args)
    SII6730_map, SII6730_int = get_emission_line_map('SII-6730', fit_results, args)
    OIII5007_map, OIII5007_int = get_emission_line_map('OIII-5007', fit_results, args)
    NII6584_map, NII6584_int = get_emission_line_map('NII-6584', fit_results, args)

    SII_map = take_safe_log_sum(SII6717_map, SII6730_map)
    SII_int = unp.log10(SII6717_int + SII6730_int)
    
    logOH_map = compute_Z_P25(OIII5007_map, NII6584_map, SII_map, args.distance_from_AGN_line_map, args)
    logOH_int = compute_Z_P25(OIII5007_int, NII6584_int, SII_int, args.distance_from_AGN_line_int, args)

    return logOH_map, logOH_int
# --------------------------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------------------------
def compute_Z_C19(ratio, coeff, ax=None, branch='high', use_C25=False, silent=False, fontsize=15):
    '''
    Calculates and returns the metallicity given observed line fluxes ratio and coefficient, according to Curti+2019
    '''
    # -----handling turnover situations, where measured ratio is beyond model peak ratio---------
    metallicity_offset = 0 if use_C25 else 8.69
    reasonable_Z_limit = [7.0, 8.4] if use_C25 else [7.6, 8.9] # Z limits within which each calibration is valid
    #reasonable_Z_limit = [6.5, 9.1]
    model = np.poly1d(coeff)
    model_diff = np.polyder(model, m=1)
    model_turnovers = np.roots(model_diff) + metallicity_offset
    possible_model_turnovers = model_turnovers[(model_turnovers > reasonable_Z_limit[0]) & (model_turnovers < reasonable_Z_limit[1])]

    if len(possible_model_turnovers) > 0:
        logOH_turnover = np.max(possible_model_turnovers) # based on solving the differential of polynomial with above coefficients, this is the value of Z where the relation peaks
        ratio_turnover = model(logOH_turnover - metallicity_offset)
    else:
        logOH_turnover, ratio_turnover = np.nan, np.nan

    if ax is not None and np.isfinite(logOH_turnover):
        for index, ax1 in enumerate(ax):
            ax1.axvline(logOH_turnover, ls='--', c='k', lw=1, label='Turnover location' if index == 0 else None)
            ax1.axhline(ratio_turnover, ls='--', c='k', lw=1)
            #ax1.fill_betweenx([-5, 5], reasonable_Z_limit[0], reasonable_Z_limit[1], color='cyan', alpha=0.1, lw=0)

    # ---------getting the line ratio limits of this diagnostics-----------
    if not silent:
        print(f'\nRatio at turnover = {ratio_turnover}')
        print(f'Ratio at min valid metallicity = {model(reasonable_Z_limit[0] - metallicity_offset)}')
        print(f'Ratio at max valid metallicity = {model(reasonable_Z_limit[1] - metallicity_offset)}')

     # --------determining data and masks------------
    if np.ma.isMaskedArray(ratio):
        ratio_arr = np.atleast_1d(ratio.data).flatten()
        mask = ratio.mask
    else:
        ratio_arr = np.atleast_1d(ratio).flatten()
        mask = None

    # --------computing the metallicitities------------
    log_OH = []
    no_solution_not_labelled = True

    for index2, this_ratio in enumerate(ratio_arr):
        this_ratio = unp.nominal_values(this_ratio)
        if np.isfinite(ratio_turnover) and this_ratio > ratio_turnover:
            log_OH.append(ufloat(logOH_turnover, 0.))
            if ax is not None:
                ax[0].axhline(this_ratio, lw=0.5, ls='solid', c='sienna', alpha=0.5, label='No solution' if no_solution_not_labelled else None)
                no_solution_not_labelled = False
        else:
            try:
                poly_to_solve = np.hstack([coeff[:-1], [coeff[-1] - this_ratio]])
                roots = np.roots(poly_to_solve)
                real_roots = np.sort(roots[np.isreal(roots)]) + metallicity_offset # see Table 1 caption in Curti+19
                possible_roots = real_roots[(real_roots > reasonable_Z_limit[0]) & (real_roots < reasonable_Z_limit[1])]
                
                impossible_roots = list(set(real_roots) - set(possible_roots))
                if ax is not None:
                    for index, real_root in enumerate(impossible_roots): ax[0].scatter(real_root, this_ratio, lw=0, s=10, c='grey')

                if branch == 'high': # THIS IS WHERE MAKING THE CHOICE TO GO WITH THE HIGHER METALLICITY BRANCH, WHEREVER THE CHOICE NEEDS TO BE MADE
                    this_log_OH = np.max(possible_roots)
                elif branch == 'low':
                    this_log_OH = np.min(possible_roots)
                log_OH.append(ufloat(this_log_OH, 0.))
            except:
                this_log_OH = np.nan
                log_OH.append(ufloat(this_log_OH, 0))
                if ax is not None:
                    ax[0].axhline(this_ratio, ls='solid', lw=0.5, c='sienna', alpha=0.5, label='No solution' if no_solution_not_labelled else None)
                    no_solution_not_labelled = False

            if np.isfinite(this_log_OH) and ax is not None:
                col_arr = ['cornflowerblue', 'limegreen', 'k']
                for index, real_root in enumerate(possible_roots): ax[0].scatter(real_root, this_ratio, ec='k', lw=0.5, s=50, c=col_arr[index] if len(real_roots) > 1 else 'salmon', zorder=100)
                try:
                    this_model_turnover = model_turnovers[-2]
                    ratio_turnover2 = max(model(this_model_turnover - metallicity_offset), model(reasonable_Z_limit[0] - metallicity_offset))
                    ax[1].scatter(this_log_OH, this_ratio, ec='k', lw=0.5, s=50, c='limegreen' if branch == 'low' and this_ratio > ratio_turnover2 else 'cornfloweblue' if branch == 'high' and this_ratio > ratio_turnover2 else 'salmon', zorder=100)
                except:
                    pass
    log_OH = np.reshape(log_OH, np.shape(ratio))
    
    if mask is None: log_OH = log_OH.tolist()
    else: log_OH = np.ma.masked_where(mask, log_OH)
    
    if ax is not None:
        handles, labels = ax[0].get_legend_handles_labels()
        fig = ax[0].figure
        fig.legend(handles, labels, bbox_to_anchor=(0.9, 0.1), loc='lower right', fontsize=fontsize)
        #ax[0].legend(loc='lower left')

    return log_OH

def get_emission_line_maps(fit_results, line_labels, args, dered=True):
    '''
    Computes and returns the spatially resolved as well as intregrated metallicity based on a given Curti+2019 calibration, from a given dictionary of emission lines
    '''
    # ----------getting line list from fits file-------------
    if not set(line_labels) <= set(args.available_lines):
        print(f'All lines in {line_labels} not available in this stack')
        return None, None

    line_map_arr, line_int_arr = [], []
    for line in line_labels:
        line_map, line_int = get_emission_line_map(line, fit_results, args, dered=dered)
        line_map_arr.append(line_map)
        line_int_arr.append(line_int)
    
    return line_map_arr, line_int_arr

def get_Z_C19(fit_results, args):
    '''
    Computes and returns the spatially resolved as well as intregrated metallicity based on a given Curti+2019 and Cataldi+2025 calibrations, from a given dictionary of emission lines
    '''
    # ------getting appropriate emission lines and calibration coefficients--------------
    if args.Zdiag == 'R3':
        line_map_arr, line_int_arr = get_emission_line_maps(fit_results, ['OIII-5007', 'H-beta'], args)
        if line_map_arr is None: return None, ufloat(np.nan, np.nan)
        
        ratio_map = take_safe_log_ratio(line_map_arr[0], line_map_arr[1])
        try: ratio_int = unp.log10(line_int_arr[0] / line_int_arr[1])
        except ValueError: ratio_int = ufloat(np.nan, np.nan)
        
        if args.use_C25: coeff = [-3.587, -4.475, 1.345, -0.08951] # from Table 4 of Cataldi+25
        else: coeff = [-0.277, -3.549, -3.593, -0.981]  # c0-3 parameters from Table 2 of Curti+19 2nd row (R3)

    elif args.Zdiag == 'N2':
        line_map_arr, line_int_arr = get_emission_line_maps(fit_results, ['NII-6584', 'H-alpha'], args)
        if line_map_arr is None: return None, ufloat(np.nan, np.nan)
        
        ratio_map = take_safe_log_ratio(line_map_arr[0], line_map_arr[1])
        try: ratio_int = unp.log10(line_int_arr[0] / line_int_arr[1])
        except ValueError: ratio_int = ufloat(np.nan, np.nan)
        
        if args.use_C25: coeff = [-6.481, 0.8163] # from Table 4 of Cataldi+25
        else: coeff = [-0.489, 1.513, -2.554, -5.293, -2.867]  # c0-4 parameters from Table 2 of Curti+19 5th row (N2)

    elif args.Zdiag == 'O3S2':
        line_map_arr, line_int_arr = get_emission_line_maps(fit_results, ['OIII-5007', 'H-beta', 'SII-6717', 'SII-6730', 'H-alpha'], args)
        if line_map_arr is None: return None, ufloat(np.nan, np.nan)
        
        SII_map = take_safe_log_sum(line_map_arr[2], line_map_arr[3], skip_log=True)
        try: SII_int = line_int_arr[2] + line_int_arr[3]
        except ValueError: SII_int = ufloat(np.nan, np.nan)

        R1_map = take_safe_log_ratio(line_map_arr[0], line_map_arr[1], skip_log=True)
        try: R1_int = line_int_arr[0] / line_int_arr[1]
        except ValueError: R1_int = ufloat(np.nan, np.nan)
                
        R2_map = take_safe_log_ratio(SII_map, line_map_arr[4], skip_log=True)
        try: R2_int = SII_int / line_int_arr[4]
        except ValueError: R2_int = ufloat(np.nan, np.nan)
        
        ratio_map = take_safe_log_ratio(R1_map, R2_map)
        try: ratio_int = unp.log10(R1_int / R2_int)
        except ValueError: ratio_int = ufloat(np.nan, np.nan)
        
        coeff = [0.191, -4.292, -2.538, 0.053, 0.332]  # c0-4 parameters from Table 2 of Curti+19 last row (O3S2)

    elif args.Zdiag == 'S2':
        line_map_arr, line_int_arr = get_emission_line_maps(fit_results, ['SII-6717', 'SII-6730', 'H-alpha'], args)
        if line_map_arr is None: return None, ufloat(np.nan, np.nan)
        
        SII_map = take_safe_log_sum(line_map_arr[0], line_map_arr[1], skip_log=True)
        try: SII_int = line_int_arr[0] + line_int_arr[1]
        except ValueError: SII_int = ufloat(np.nan, np.nan)
        
        ratio_map = take_safe_log_ratio(SII_map, line_map_arr[2])
        try: ratio_int = unp.log10(SII_int / line_int_arr[2])
        except ValueError: ratio_int = ufloat(np.nan, np.nan)
        
        coeff = [-0.442, -0.360, -6.271, -8.339, -3.559]  # c0-3 parameters from Table 2 of Curti+19 3rd-to-last row (S2)

    elif args.Zdiag == 'RS32':
        line_map_arr, line_int_arr = get_emission_line_maps(fit_results, ['OIII-5007', 'H-beta', 'SII-6717', 'SII-6730', 'H-alpha'], args)
        if line_map_arr is None: return None, ufloat(np.nan, np.nan)

        SII_map = take_safe_log_sum(line_map_arr[2], line_map_arr[3], skip_log=True)
        try: SII_int = line_int_arr[2] + line_int_arr[3]
        except ValueError: SII_int = ufloat(np.nan, np.nan)
        
        R1_map = take_safe_log_ratio(line_map_arr[0], line_map_arr[1], skip_log=True)
        try: R1_int = line_int_arr[0] / line_int_arr[1]
        except ValueError: R1_int = ufloat(np.nan, np.nan)
        
        R2_map = take_safe_log_ratio(SII_map, line_map_arr[4], skip_log=True)
        try: R2_int = SII_int / line_int_arr[4]
        except ValueError: R2_int = ufloat(np.nan, np.nan)
        
        ratio_map = take_safe_log_sum(R1_map, R2_map)
        try: ratio_int = unp.log10(R1_int + R2_int)
        except ValueError: ratio_int = ufloat(np.nan, np.nan)
        
        coeff = [-0.054, -2.546, -1.970, 0.082, 0.222]  # c0-3 parameters from Table 2 of Curti+19 2nd-to-last row (RS32)

    else:
        print(f'Could not apply any of the metallicity diagnostics, so returning NaN metallicities')
        return None, ufloat(np.nan, np.nan)

    coeff = coeff[::-1] # because in Curti+2019 the coefficients are listed in the reverse order compared to what np.poly1d prefers

    # ----setting up Z diag debugging plots---------
    if args.debug_Zdiag:
        fig, ax = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
        ax[0].set_xlabel('Possible values of log(O/H)+12', fontsize=args.fontsize)
        if args.Zbranch == 'high': ax[1].set_xlabel('log(O/H)+12 = max(solution)', fontsize=args.fontsize)
        elif args.Zbranch == 'low': ax[1].set_xlabel('log(O/H)+12 = min(solution)', fontsize=args.fontsize)
        ax[0].set_ylabel(f'Observed log {args.Zdiag}', fontsize=args.fontsize)
        allowed_Z_limit = [7.0, 8.4] if args.use_C25 else [7.6, 8.9] # Z limits within which each calibration is valid
        Z_limits = [6.5, 9.5]
        ratio_limits = [-0., 1.5]

        metallicity_offset = 0 if args.use_C25 else 8.69
        xarr_valid = np.linspace(allowed_Z_limit[0], allowed_Z_limit[1], 100)
        ax[0].plot(xarr_valid, np.poly1d(coeff)(xarr_valid - metallicity_offset), lw=3, c='k', ls='solid', label='Valid C25 calibration' if args.use_C25 else 'Valid C20 calibration')
        ax[1].plot(xarr_valid, np.poly1d(coeff)(xarr_valid - metallicity_offset), lw=3, c='k', ls='solid')

        xarr_full = np.linspace(Z_limits[0], Z_limits[1], 100)
        ax[0].plot(xarr_full, np.poly1d(coeff)(xarr_full - metallicity_offset), lw=0.5, c='k', ls='dotted', label='Extrapolated C25 calibration' if args.use_C25 else 'Extrapolated C20 calibration')
        ax[1].plot(xarr_full, np.poly1d(coeff)(xarr_full - metallicity_offset), lw=0.5, c='k', ls='dotted')

        # for ax1 in ax:
        #     ax1.fill_betweenx([-5, 5], Z_limits[0], allowed_Z_limit[0], color='grey', lw=0, alpha=0.2, label='Calibration invalid')
        #     ax1.fill_betweenx([-5, 5], allowed_Z_limit[1], Z_limits[1], color='grey', lw=0, alpha=0.2)

        ax[0].set_ylim(ratio_limits[0], ratio_limits[1])
        ax[0].set_xlim(Z_limits[0], Z_limits[1])
        ax[1].set_xlim(Z_limits[0], Z_limits[1])
    else:
        ax = None

    # -------estimating the metallicities---------------
    logOH_map = compute_Z_C19(ratio_map, coeff, ax=ax, branch=args.Zbranch, use_C25=args.use_C25, silent=True, fontsize=args.fontsize)
    logOH_int = compute_Z_C19(ratio_int, coeff, ax=ax, branch=args.Zbranch, use_C25=args.use_C25, silent=True, fontsize=args.fontsize)

    # -------saving the debugging plots---------------
    if ax is not None:
        Zbranch_text = '' if args.Zdiag in ['NB', 'P25', 'Te', 'N2'] else f'-{args.Zbranch}'
        figname = args.output_dir / 'plots' / f'stacked_metallicity_debug_Zdiag_{args.Zdiag}{Zbranch_text}.png'
        fig.savefig(figname, transparent=args.fortalk, dpi=200)
        print(f'\nSaved figure at {figname}')

    return logOH_map, logOH_int
# --------------------------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------------------------
def compute_Z_NB(line_label_array, line_waves_array, line_flux_array, args):
    '''
    Calculates and returns the NebulaBayes metallicity given a list of observed line fluxes
    '''
    line_flux_array = [np.atleast_1d(item) for item in line_flux_array]
    
    map_shape = np.shape(line_flux_array)[1:]
    if len(map_shape) == 1: npixels = map_shape[0]
    else: npixels = map_shape[0] * map_shape[1]

    IDs_array = np.arange(npixels).flatten()
    unique_IDs_array = np.unique(IDs_array)
    
    print(f'\nAbout to start running NB, with {len(line_label_array)} lines: {line_label_array}..\n')
    if len(unique_IDs_array) > 60: print(f'This might take ~{int(len(unique_IDs_array) / 60)} min')

    # -----making a "net" mask array and separating out the line fluxes form the input unumpy arrays---------
    net_mask = np.zeros(np.shape(line_flux_array[0]), dtype=bool)
    obs_flux_array, obs_err_array = [], []

    for index in range(len(line_flux_array)):
        if np.ma.isMaskedArray(line_flux_array[index]):
            net_mask = net_mask | line_flux_array[index].mask
            obs_flux_array.append(unp.nominal_values(line_flux_array[index].data).flatten())
            obs_err_array.append(unp.std_devs(line_flux_array[index].data).flatten())
        else:
            obs_flux_array.append(unp.nominal_values(line_flux_array[index]).flatten())
            obs_err_array.append(unp.std_devs(line_flux_array[index]).flatten())

    obs_flux_array = np.array(obs_flux_array)
    obs_err_array = np.array(obs_err_array)
    net_mask_array = net_mask.flatten()

    # -----loading the NB HII region model grid---------
    NB_Model_HII = NB_Model('HII', line_list=line_label_array)

    out_dir = args.output_dir / 'NB_results' / f'{args.id}_NB_results'
    out_subdirs = [out_dir / 'prior_plots', out_dir / 'likelihood_plots', out_dir / 'posterior_plots', out_dir / 'best_model_catalogs', out_dir / 'param_estimates_catalogs']
    for this_out_subdir in out_subdirs: this_out_subdir.mkdir(exist_ok=True, parents=True)

    # -----looping over each pixel to calculate NB metallicity--------
    logOH_array = []
    logOH_dict_unique_IDs = {}
    counter = 0
    start_time3 = datetime.now()

    for index in range(npixels):
        this_ID = IDs_array[index]
        if net_mask_array[index]: # no need to calculate for those pixels that are already masked
            print(f'Skipping NB for masked pixel {this_ID} ({index + 1}/{npixels})..')
            logOH = ufloat(np.nan, np.nan)
        elif this_ID in logOH_dict_unique_IDs.keys():
            print(f'Skipping NB due to existing measurement from unique ID {this_ID} ({index + 1}/{npixels})..')
            logOH = logOH_dict_unique_IDs[this_ID]
        else:
            start_time4 = datetime.now()

            # ------getting all line fluxes-------------
            obs_fluxes = obs_flux_array[:,index]
            obs_errs = obs_err_array[:, index]

            # ------discarding lines with negative fluxes-------------
            good_obs = obs_fluxes >= 0
            obs_fluxes = obs_fluxes[good_obs]
            obs_errs = obs_errs[good_obs]
            line_labels = list(np.array(line_label_array)[good_obs])
            line_waves = np.array(line_waves_array)[good_obs]

            if len(line_labels) > 1:
                # -------setting up NB parameters----------
                dered = 'Hbeta' in line_labels and 'Halpha' in line_labels and args.dered_in_NB
                norm_line = 'Hbeta' if 'Hbeta' in line_labels else 'OIII5007' if 'OIII5007' in line_labels else 'NII6583_Halpha' if 'NII6583_Halpha' in line_labels else 'Halpha' if 'Halpha' in line_labels else line_labels[0]
                kwargs = {'prior_plot': os.path.join(out_dir, 'prior_plots', f'{this_ID}_HII_prior_plot.pdf'),
                        'likelihood_plot': os.path.join(out_dir, 'likelihood_plots', f'{this_ID}_HII_likelihood_plot.pdf'),
                        'posterior_plot': os.path.join( out_dir, 'posterior_plots', f'{this_ID}_HII_posterior_plot.pdf'),
                        'estimate_table': os.path.join(out_dir, 'best_model_catalogs', f'{this_ID}_HII_param_estimates.csv'),
                        'best_model_table': os.path.join(out_dir, 'param_estimates_catalogs', f'{this_ID}_HII_best_model.csv'),
                        'verbosity': 'ERROR',
                        'norm_line':norm_line,
                        'deredden': dered,
                        'propagate_dered_errors': dered,
                        'obs_wavelengths': line_waves if dered else None
                        }

                # -------running NB--------------
                print(f'Deb1576: binID {this_ID} ({index + 1}/{npixels}): nlines={len(obs_fluxes)}, {dict(zip(line_labels, obs_fluxes))}, norm_line = {norm_line}, dereddening on the fly? {dered}') ##
                Result = NB_Model_HII(obs_fluxes, obs_errs, line_labels, **kwargs)

                # -------estimating the resulting logOH, and associated uncertainty-----------
                df_estimates = Result.Posterior.DF_estimates # pandas DataFrame
                logOH_est = df_estimates.loc['12 + log O/H', 'Estimate']
                logOH_low = df_estimates.loc['12 + log O/H', 'CI68_low']
                logOH_high = df_estimates.loc['12 + log O/H', 'CI68_high']
                logOH_err = np.mean([logOH_est - logOH_low, logOH_high - logOH_est])
                logOH = ufloat(logOH_est, logOH_err)
                
                print(f'Ran NB for unique ID {this_ID} ({counter + 1} out of {len(unique_IDs_array)}) with {len(obs_fluxes)} good fluxes in {timedelta(seconds=(datetime.now() - start_time4).seconds)}')            
            else:
                logOH = ufloat(np.nan, np.nan)
                print(f'Could not run NB for unique ID {this_ID} ({counter + 1} out of {len(unique_IDs_array)}) with only {len(obs_fluxes)} good fluxes')

            counter += 1
            logOH_dict_unique_IDs.update({this_ID: logOH}) # updating to unique ID dictionary once logOH has been calculated for this unique ID

        logOH_array.append(logOH)
    print(f'\nRan NB for total {counter} unique pixels out of {len(obs_flux_array[0])}, in {timedelta(seconds=(datetime.now() - start_time3).seconds)}\n')

    # ---------collating all the metallicities computed---------
    log_OH = np.ma.masked_where(net_mask, np.reshape(logOH_array, np.shape(line_flux_array[0])))
    if len(log_OH) == 1: log_OH = log_OH.data[0]

    return log_OH

def get_Z_NB(fit_results, args):
    '''
    Computes and returns the spatially resolved as well as intregrated metallicity from a given line flux dictionary, based on Bayesian
    statistics, using NebulaBayes
    '''
    # -----dict for converting line label names to those acceptable to NB---------
    if args.use_original_NB_grid: line_label_dict = {'H-beta':'Hbeta', 'OIII-5007':'OIII5007', 'OIII-4363':'OIII4363', 'H-alpha':'Halpha', 'NII-6584':'NII6583', 'SII-6717':'SII6716', 'SII-6730':'SII6731'}
    else: line_label_dict = {'H-beta':'Hbeta', 'OIII-5007':'OIII5007', 'OIII-4363':'OIII4363', 'H-alpha':'Halpha', 'NII-6584':'NII6583', 'SII-6717':'SII6716', 'SII-6730':'SII6731'}

    line_map_array, line_int_array, line_label_array, line_waves_array = [], [], [], []
    for line in args.available_lines:
        line_map, line_int = get_emission_line_map(line, fit_results, args, dered=not args.dered_in_NB)

        if line in line_label_dict and line not in args.exclude_lines:
            line_map_array.append(line_map)
            line_int_array.append(line_int)
            line_label_array.append(line_label_dict[line])
            line_waves_array.append(rest_wave_dict[line]) # in Angstroms

    # ----------calling NB----------------------
    logOH_map = compute_Z_NB(line_label_array, line_waves_array, line_map_array, args)
    logOH_int = compute_Z_NB(line_label_array, line_waves_array, line_int_array, args)

    return logOH_map, logOH_int
# --------------------------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------------------------
def get_metallicity_map(fit_results, args):
    '''
    Computes 2D metallicity and integrated metallicity from a given dictionary of several stacked lines 
    Returns metallicity map (2D array) along with uncertainty
    '''
    if args.Zdiag == 'NB':
        logOH_map, logOH_int = get_Z_NB(fit_results, args)
    elif args.Zdiag == 'P25' and set(['OIII-5007', 'H-alpha', 'SII-6717', 'SII-6730']) <= set(args.available_lines):
        logOH_map, logOH_int = get_Z_P25(fit_results, args)
    else:
        logOH_map, logOH_int = get_Z_C19(fit_results, args)

    return logOH_map, logOH_int

def compute_SFR(Ha_flux, distance):
    '''
    Calculates and returns the SFR given observed H alpha fluxes
    Conversion factor is from Kennicutt 1998 (Eq 2 of https://ned.ipac.caltech.edu/level5/Sept01/Rosa/Rosa3.html)
    '''
    # -----handling masks separately because uncertainty package cannot handle masks---------
    if np.ma.isMaskedArray(Ha_flux):
        net_mask = Ha_flux.mask
        Ha_flux = Ha_flux.data
    else:
        net_mask = False

    Ha_lum = Ha_flux * 4 * np.pi * (distance.to('cm').value) ** 2 # converting to ergs/s (luminosity)
    #sfr = Ha_lum * 7.9e-42 # luminosity in ergs/s; SFR in Msun/yr; from Kennicutt+1998
    #sfr = Ha_lum * 10**(-41.67) # luminosity in ergs/s; SFR in Msun/yr; from Reddy+2022
    sfr = Ha_lum * ufloat(7.5, 1.3) * 10 ** (-42) # luminosity in ergs/s; SFR in Msun/yr; from Shivaei+2015

    if hasattr(Ha_flux, "__len__"): # if it is an array
        sfr = np.ma.masked_where(net_mask, sfr)

    return sfr
# --------------------------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------------------------
def update_quant_dict(quant, quant_map, quant_int):
    '''
    Add the spatially resolved map and integrated values (and corresponding uncertainties) of a given quant into a given quant dictionary
    '''
    if np.ma.isMaskedArray(quant_map):
        mask = quant_map.mask
        quant_map = quant_map.data
    else:
        mask = False

    new_dict = {}
    for this_quant in [f'{quant}_mask', f'{quant}', f'{quant}_err']: new_dict.update({this_quant:{}})
    
    new_dict[f'{quant}_mask']['map'] = mask
    new_dict[f'{quant}']['map'] = unp.nominal_values(quant_map)
    new_dict[f'{quant}_err']['map'] = unp.std_devs(quant_map)
    new_dict[f'{quant}']['int'] = quant_int.n
    new_dict[f'{quant}_err']['int'] = quant_int.s

    return new_dict

def compute_quant_maps(fit_results, args):
    '''
    Computes metallicity, SFR and other parameters from a given set of emission line maps
    Returns dictionary containing the quantity maps (and their corresponding uncertainty maps)
    '''
    quant_maps = {}

    # -----------get metallicity map-------------
    logOH_map, logOH_int = get_metallicity_map(fit_results, args)
    quant_maps.update(update_quant_dict('logOH', logOH_map, logOH_int))

    # -----------get SFR map----------
    Halpha_map, Halpha_int = get_emission_line_map('H-alpha', fit_results, args)
    sfr_map = compute_SFR(Halpha_map, args.distance)
    sfr_int = compute_SFR(Halpha_int, args.distance)
    quant_maps.update(update_quant_dict('sfr', sfr_map, sfr_int))

    return quant_maps
# --------------------------------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------------------------------
def plot_quant_map(ax, quant, quant_maps, args, cmap='cividis', clabel='', takelog=True, vmin=None, vmax=None, hide_xaxis=False, hide_yaxis=False, hide_cbar=True, hide_cbar_ticks=False, cticks_integer=True):
    '''
    Extracts the quantity map from the given dict of quant_maps, and plots it on the given axis handle
    Returns axis handle
    '''
    quant_map, _, quant_int, quant_mask = get_quant_from_quant_maps(quant, quant_maps)
    quant_masked = np.ma.masked_where(quant_mask, quant_map)
    ax = plot_2D_map(quant_masked, ax, f'{quant}: {quant_int:.2f}', args, cmap=cmap, clabel=clabel, takelog=takelog, vmin=vmin, vmax=vmax, hide_xaxis=hide_xaxis, hide_yaxis=hide_yaxis, hide_cbar=hide_cbar, hide_cbar_ticks=hide_cbar_ticks, cticks_integer=cticks_integer)

    return ax

def plot_met_sfr_corr(ax, quant_maps, args, color='salmon', colorby=None, cmap='cividis', hide_xaxis=False, hide_yaxis=False, clabel='', hide_cbar=True, hide_cbar_ticks=False, cticks_integer=True, log_sfr_min=None, log_sfr_max=None, logOH_min=None, logOH_max=None):
    '''
    Extracts the metallicity and SFR maps from the given dict of quant_maps, and plots their correlation on the given axis handle
    Returns axis handle
    '''
    logOH_map, logOH_map_err, _, logOH_mask = get_quant_from_quant_maps('logOH', quant_maps)
    sfr_map, sfr_map_err, _, sfr_mask = get_quant_from_quant_maps('sfr', quant_maps)

    # ----------culling only area under galaxy-------------------
    ny, nx = np.shape(logOH_map)
    crop = np.s_[int(ny/2 - args.upto_pix): int(ny/2 + args.upto_pix), int(nx/2 - args.upto_pix): int(nx/2 + args.upto_pix)]

    logOH_map = logOH_map[crop]
    logOH_map_err = logOH_map_err[crop]
    logOH_mask = logOH_mask[crop]

    sfr_map = sfr_map[crop]
    sfr_map_err = sfr_map_err[crop]
    sfr_mask = sfr_mask[crop]

    # -----------making and curtailing dataframe-------------
    df = pd.DataFrame({'logOH': logOH_map.flatten().astype(np.float64), 'logOH_u': logOH_map_err.flatten().astype(np.float64), 'logOH_mask': logOH_mask.flatten().astype(bool),
                       'sfr': sfr_map.flatten().astype(np.float64), 'sfr_u': sfr_map_err.flatten().astype(np.float64), 'sfr_mask': sfr_mask.flatten().astype(bool),
                       })
    ndf = len(df)

    df = df[~(df['logOH_mask'] | df['sfr_mask'])]
    log_quant = unp.log10(unp.uarray(df['sfr'], df['sfr_u']))
    df['log_sfr'] = unp.nominal_values(log_quant)
    df['log_sfr_u'] = unp.std_devs(log_quant)
    df = df.drop(columns=['logOH_mask', 'sfr_mask', 'sfr', 'sfr_u'])
    df = df.dropna(subset=['logOH', 'log_sfr']).reset_index(drop=True)
    
    for col in ['log_sfr_u']: df = df[df[col] < np.percentile(df[col], 99)] # curtailing extreme outliers in uncertainty
    

    print(f'{len(df)} of {ndf} spaxels remain with finite logOH and SFR')

    # -----------plotting scatter----------------------
    if colorby is None: color = color
    else: color = df[colorby]

    p = ax.scatter(df['log_sfr'], df['logOH'], s=10, lw=1, ec='k', color=color)
    ax.errorbar(df['log_sfr'], df['logOH'], xerr=df['log_sfr_u'], yerr=df['logOH_u'], c=color, fmt='none', lw=0.5, alpha=0.7, zorder=-5)

    # -----------------annotating axes---------------------
    ax.set_xlim(log_sfr_min, log_sfr_max)
    ax.set_ylim(logOH_min, logOH_max)

    ax = annotate_axes(ax, r'Log SFR (M$_\odot$/yr)', r'$\log$ O/H + 12', fontsize=args.fontsize / args.fontfactor, clabel=clabel, hide_xaxis=hide_xaxis, hide_yaxis=hide_yaxis, hide_cbar=hide_cbar, p=p, hide_cbar_ticks=hide_cbar_ticks, cticks_integer=cticks_integer)

    return ax
# -------------------------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------------------------
def save_quant_maps_fits(quant_maps, quants_fits_file, wcs, args):
    '''
    Saves the N x 2D spatially resolved quantity maps as N-extension fits files, along with other relevant header info
    '''
    # ------------setting up primary header---------------
    spatial_header = wcs.to_header()
    primary_hdu = fits.PrimaryHDU(header=spatial_header)

    primary_hdu.header['CONTENT'] = 'Various physical quantity maps'
    primary_hdu.header['ID'] = args.id
    primary_hdu.header['REDSHIFT'] = args.z
    primary_hdu.header['Z_DIAGNOSTIC'] = args.Zdiag
    
    hdul = fits.HDUList([primary_hdu])

    params_dict = {
        'logOH':{'label':'LOG_OH', 'unit':''},
        'logOH_err':{'label':'LOG_OH_ERR', 'unit':''},
        'logOH_mask':{'label':'LOG_OH_MASK', 'unit':''},
        'sfr':{'label':'SFR', 'unit':'Msun/s'},
        'sfr_err':{'label':'SFR_ERR', 'unit':'Msun/s'},
        'sfr_mask':{'label':'SFR_MASK', 'unit':''},
    }

    # ------------looping over fitted lines---------------
    for label in list(quant_maps.keys()):   
        quant_data = quant_maps[label]

        data = quant_data['map']
        if data is None: continue
            
        hdu = fits.ImageHDU(data=data.astype(np.float32), header=spatial_header)
        hdu.header['EXTNAME'] = params_dict[label]['label']
        hdu.header['BUNIT'] = params_dict[label]['unit']
        if 'int' in quant_data:
            quant_int = quant_data['int']
            if np.isnan(quant_int): hdu.header[f'INT_{params_dict[label]["label"]}'] = 'NAN'
            else: hdu.header[f'INT_{params_dict[label]["label"]}'] = quant_int
        hdul.append(hdu)

    hdul.writeto(quants_fits_file, overwrite=True)
    print(f'Successfully saved {len(hdul)-1} extensions to {quants_fits_file}"')
    return

def read_quant_maps_fits(filename):
    '''
    Reads a multi-extension FITS file and reconstructs the quant_maps dictionary.
    Returns:
    quant_maps: dict. Format {label: {param: 2D array}}
    spatial_header: astropy.io.fits.Header. The WCS/spatial info
    '''
    print(f'Reading existing maps fits file from {filename}')
    quant_maps = {}
    
    hdul = fits.open(filename)
    spatial_header = hdul[0].header
    suffix_to_key = {
        'LOG_OH': 'logOH',
        'LOG_OH_ERR': 'logOH_err',
        'LOG_OH_MASK': 'logOH_mask',
        'SFR': 'sfr',
        'SFR_ERR': 'sfr_err',
        'SFR_MASK': 'sfr_mask',
    }
    
    # Loop through all extensions starting from index 1
    for index in range(1, len(hdul)):
        hdu = hdul[index]
        extname = hdu.name
        if extname not in suffix_to_key: continue
        header = hdu.header

        label = suffix_to_key[extname]
        quant_maps.update({label:{}})
        quant_maps[label]['map'] = hdu.data.copy()
        if f'INT_{extname}' in header:
            quant_int = header[f'INT_{extname}']
            if quant_int == 'NAN': quant_maps[label]['int'] = np.nan
            else: quant_maps[label]['int'] = quant_int

    return quant_maps, spatial_header

def get_quant_from_quant_maps(quant, quant_maps):
    '''
    Extracts maps and integrated values for a given quant from a given quant dictionary
    '''
    mask = quant_maps[f'{quant}_mask']['map']
    map = quant_maps[f'{quant}']['map']
    map_err = quant_maps[f'{quant}_err']['map']
    integrated = ufloat(quant_maps[f'{quant}']['int'], quant_maps[f'{quant}_err']['int'])

    return map, map_err, integrated, mask
   
# --------------------------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    args.fontfactor = 1.2
    args.id_arr = args.id

    #logOH_min, logOH_max = 7.0, 8.0
    logOH_min, logOH_max = None, None
    log_sfr_min, log_sfr_max = 16, 21
    #log_sfr_min, log_sfr_max = None, None
    
    # -------------setup directories and global variables----------------
    maps_fits_dir = args.output_dir / 'maps'
    quants_fits_dir = args.output_dir / 'quants'
    args.fig_dir = args.output_dir / 'plots'
    quants_fits_dir.mkdir(exist_ok=True, parents=True)

    catalog_file = args.input_dir / 'redshifts.dat'

    # ----------------reading in catalog---------------------
    df = read_msa3d_catalog(catalog_file)
    if not args.do_all_obj: df = df[df['id'].isin(args.id_arr)].reset_index(drop=True)

    # ----------------looping over the objects in this chunk-------------
    for index, obj in df.iterrows():
        start_time2 = datetime.now()
        print(f'Commencing ({index + 1}/{len(df)}) ID {obj["id"]}..')
        args.id =obj['id']
        args.z = obj['redshift']
        args.distance = cosmo.luminosity_distance(args.z)
        args.upto_arcsec = args.upto_kpc * cosmo.arcsec_per_kpc_proper(args.z).value # arcsec
        args.upto_pix = args.upto_arcsec / msa_pix_size_arcsec

        # ------determining directories and filenames---------
        args.maps_fits_file = maps_fits_dir / f'{args.id:05d}.maps.fits'
        args.quants_fits_file = quants_fits_dir / f'{args.id:05d}_Zdiag_{args.Zdiag}.quants.fits'

        # ---------measuring the various quantitites--------
        if not os.path.exists(args.quants_fits_file) or args.clobber:
            # -----------read in the emission line maps--------------
            fit_results, spatial_header = read_line_maps_fits(args.maps_fits_file)
            wcs = pywcs.WCS(spatial_header)

            args.available_lines = list(fit_results.keys())
            if 'H-alpha' not in args.available_lines:
                print(f'ID {args.id} (z={args.z:.2f}) does not have H-alpha, so skipping it..')
                continue
            if 'H-beta' not in args.available_lines:
                print(f'ID {args.id} (z={args.z:.2f}) does not have H-beta, so skipping it..')
                continue

            EB_V_map, args.EB_V_int = get_EB_V(fit_results, args)
            args.EB_V_map = EB_V_map.data

            # -----------computing various quantities--------------
            quant_maps = compute_quant_maps(fit_results, args)

            # -----------save the quanttity maps in fits file-------------
            save_quant_maps_fits(quant_maps, args.quants_fits_file, wcs, args)
        else:
            quant_maps, spatial_header = read_quant_maps_fits(args.quants_fits_file)
        
        # --------plot the qunatity maps-------------
        if args.plot_metallicity:
            fig, ax = plt.subplots(1, figsize=(8, 8), layout='constrained')
            
            ax = plot_quant_map(ax, 'logOH', quant_maps, args, cmap='cividis', takelog=False, vmin=logOH_min, vmax=logOH_max, hide_cbar=False)
            
            fig.text(0.1, 0.98, f'ID {args.id}', fontsize=args.fontsize, c='k', ha='left', va='top')
            save_fig(fig, args.fig_dir, f'{args.id}_metallicity_{args.Zdiag}_maps.png', args)    

        if args.plot_met_sfr:
            fig, axes = plt.subplots(1, 3, figsize=(10, 4.5))
            fig.subplots_adjust(left=0.06, right=0.98, bottom=0.12, top=0.9, wspace=0.5)
            
            axes[0] = plot_quant_map(axes[0], 'logOH', quant_maps, args, cmap='cividis', takelog=False, vmin=logOH_min, vmax=logOH_max, hide_cbar=False)
            axes[1] = plot_quant_map(axes[1], 'sfr', quant_maps, args, cmap='Blues', takelog=True, vmin=log_sfr_min, vmax=log_sfr_max, hide_yaxis=True, hide_cbar=False)
            axes[2] = plot_met_sfr_corr(axes[2], quant_maps, args, color='salmon', log_sfr_min=log_sfr_min, log_sfr_max=log_sfr_max, logOH_min=logOH_min, logOH_max=logOH_max)
            
            fig.text(0.05, 0.95, f'ID {args.id}', fontsize=args.fontsize, c='k', ha='left', va='top')
            save_fig(fig, args.fig_dir, f'{args.id}_metallicity_{args.Zdiag}_SFR_corr.png', args)    

        print(f'\nCompleted ID {args.id} in {timedelta(seconds=(datetime.now() - start_time2).seconds)}, {len(df) - index - 1} to go!')

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
