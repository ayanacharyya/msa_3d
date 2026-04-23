'''
    Filename: plot_relations.py
    Notes: Performs metallicity and SFR maps from MSA-3D emission line maps and stores them in fits files
    Author : Ayan
    Created: 30-03-26
    Example: run plot_relations.py --Zdiag NB --snr_cut 3 --fit_correlation
'''

from header import *
from util import *
from make_metallicity_sfr_maps import odr_fit, plot_fitted_line

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def make_plot(xcol, ycol, df, ax, args, colorcol=None, cmap='plasma', color='cornflowerblue', size=50, do_fit=False):
    '''
    Accepts dataframe, and column names to plot, and plots those on a given axis handle, then annotates the axes based on label dictionaries
    Returns the axis
    '''
    label_dict = {'redshift':'Redshift', 'log_mass':'Log stellar mass', 'log_sfr':'Log SFR', 'Zgrad_kpc':r'$\nabla_r$Z', 'Zgrad_re':r'$\nabla_r$Z', 'logZ_logSFR_slope':r'$\nabla_r$Z - $\Sigma_{\rm SFR}$ slope', 'vdisp_mean':r'Mean $\sigma_{\rm disp}$', 'vdisp_50':r'Median $\sigma_{\rm disp}$', 't_mix':r'$t_{\rm mix}$'}
    unit_dict = {'redshift':'', 'log_mass':r'M$_{\odot}$', 'log_sfr':r'M$_{\odot}$/yr', 'Zgrad_kpc':'dex/kpc', 'Zgrad_re':r'dex/r$_e$', 'logZ_logSFR_slope':'', 'vdisp_mean':'km/s', 'vdisp_50':'km/s', 't_mix':'yr'}

    if args.Zdiag == 'NB':
        vmin_dict = {'redshift':0.9, 'log_mass':8.0, 'log_sfr':-4.0, 'Zgrad_kpc':-2, 'Zgrad_re':-2, 'logZ_logSFR_slope':-60, 'vdisp_mean':50, 'vdisp_50':50, 't_mix':0}
        vmax_dict = {'redshift':1.9, 'log_mass':11.0, 'log_sfr':1.0, 'Zgrad_kpc':2, 'Zgrad_re':2, 'logZ_logSFR_slope':20, 'vdisp_mean':150, 'vdisp_50':150, 't_mix':350}
        
    elif args.Zdiag == 'R23':
        vmin_dict = {'redshift':0.9, 'log_mass':8.0, 'log_sfr':-4.0, 'Zgrad_kpc':-2, 'Zgrad_re':-2, 'logZ_logSFR_slope':-0.5, 'vdisp_mean':50, 'vdisp_50':50, 't_mix':-250}
        vmax_dict = {'redshift':1.9, 'log_mass':11.0, 'log_sfr':1.0, 'Zgrad_kpc':2, 'Zgrad_re':2, 'logZ_logSFR_slope':0.7, 'vdisp_mean':150, 'vdisp_50':150, 't_mix':250}

    # vmin_dict = {'redshift':0.9, 'log_mass':8.0, 'log_sfr':-4.0, 'Zgrad_kpc':-2, 'Zgrad_re':-2, 'logZ_logSFR_slope':-0.3, 'vdisp_mean':None, 'vdisp_50':None, 't_mix':-70}
    # vmax_dict = {'redshift':1.9, 'log_mass':11.0, 'log_sfr':1.0, 'Zgrad_kpc':2, 'Zgrad_re':2, 'logZ_logSFR_slope':0.4, 'vdisp_mean':None, 'vdisp_50':None, 't_mix':100}

    p = ax.scatter(df[xcol], df[ycol], s=size, lw=1, ec='k', c=df[colorcol] if colorcol is not None else color)

    xerr = [df['vdisp_50'] - df['vdisp_16'], df['vdisp_84'] - df['vdisp_50']] if xcol == 'vdisp_50' else df[xcol + '_u'] if xcol + '_u' in df else None
    yerr = [df['vdisp_50'] - df['vdisp_16'], df['vdisp_84'] - df['vdisp_50']] if ycol == 'vdisp_50' else df[ycol + '_u'] if ycol + '_u' in df else None
    ax.errorbar(df[xcol], df[ycol], yerr=yerr, xerr=xerr, fmt='none', lw=1, color=color, alpha=0.7, zorder=-5)

    # -----------fitting scatter plot----------------------
    if do_fit:
        linefit_odr = odr_fit(df, quant_x=xcol, quant_y=ycol, recenter_y=False)
        xarr = np.linspace(vmin_dict[xcol] * 1.10 if vmin_dict[xcol] is not None else df[xcol].min() * 0.8, 
                           vmax_dict[xcol] * 0.95 if vmax_dict[xcol] is not None else df[xcol].max() * 1.2, 
                           20)
        ax = plot_fitted_line(ax, linefit_odr, xarr, 'salmon', args, quant=ycol, label=f'Slope = {linefit_odr[0]: .2f}')

    # ----------annotating plot----------------------------------
    ax = annotate_axes(ax, f'{label_dict[xcol]} ({unit_dict[xcol]})', f'{label_dict[ycol]} ({unit_dict[ycol]})', xlim=[vmin_dict[xcol], vmax_dict[xcol]], ylim=[vmin_dict[ycol], vmax_dict[ycol]], args=args, hide_cbar=colorcol is None, p=p, clabel=label_dict[colorcol] if colorcol is not None else '')

    return ax

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    args.fontfactor = 1.2
    args.id_arr = args.id
    
    # -------------setup directories and global variables----------------
    args.fig_dir = args.output_dir / 'plots'
    tie_vdisp_text = '_tie_vdisp' if args.tie_vdisp else ''
    snr_cut_text = f'_snr{args.snr_cut}'
    Z_SFR_slope_file = args.output_dir / 'catalogs' / f'Z_{args.Zdiag}_SFR_slopes{tie_vdisp_text}{snr_cut_text}.csv'

    # -----------read dataframe-------------------
    df = pd.read_csv(Z_SFR_slope_file)
    def_omit_ids = [9527, 7561]
    maybe_omit_ids = [9337, 8512, 7314, 2145]
    omit_ids = def_omit_ids + maybe_omit_ids
    #df = df[~df['id'].isin(omit_ids)]

    # -------make plots-----------
    fig, axes = plt.subplots(2, 2, figsize=(10, 6))
    fig.subplots_adjust(left=0.1, right=0.90, bottom=0.12, top=0.98, wspace=0.6, hspace=0.3)
    axes = axes.flatten()

    axes[0] = make_plot('log_mass', 'logZ_logSFR_slope', df, axes[0], args, colorcol='redshift', do_fit=False)
    axes[1] = make_plot('log_mass', 't_mix', df, axes[1], args, colorcol='redshift', do_fit=False)
    axes[2] = make_plot('vdisp_mean', 't_mix', df, axes[2], args, colorcol='log_mass', do_fit=args.fit_correlation)
    axes[3] = make_plot('vdisp_50', 't_mix', df, axes[3], args, colorcol='log_mass', do_fit=args.fit_correlation)

    # ----------save figure---------------
    figname = f'all_relations_Z_{args.Zdiag}{tie_vdisp_text}{snr_cut_text}.png'
    save_fig(fig, args.fig_dir, figname, args)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
