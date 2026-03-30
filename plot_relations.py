'''
    Filename: plot_relations.py
    Notes: Performs metallicity and SFR maps from MSA-3D emission line maps and stores them in fits files
    Author : Ayan
    Created: 30-03-26
    Example: run plot_relations.py --Zdiag NB
'''

from header import *
from util import *

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def make_plot(xcol, ycol, df, ax, args, colorcol=None, cmap='plasma', color='cornflowerblue', size=50):
    '''
    Accepts dataframe, and column names to plot, and plots those on a given axis handle, then annotates the axes based on label dictionaries
    Returns the axis
    '''
    label_dict = {'redshift':'Redshift', 'log_mass':'Log stellar mass', 'log_sfr':'Log SFR', 'Zgrad_kpc':r'$\nabla_r$Z', 'Zgrad_re':r'$\nabla_r$Z', 'logZ_logSFR_slope':r'$\nabla_r$Z - $\Sigma_{\rm SFR}$ slope', 'vdisp_mean':r'Mean $\sigma_{\rm disp}$', 'vdisp_50':r'Median $\sigma_{\rm disp}$', 't_mix':r'$t_{\rm mix}$'}
    unit_dict = {'redshift':'', 'log_mass':r'M$_{\odot}$', 'log_sfr':r'M$_{\odot}$/yr', 'Zgrad_kpc':'dex/kpc', 'Zgrad_re':r'dex/r$_e$', 'logZ_logSFR_slope':'', 'vdisp_mean':'km/s', 'vdisp_50':'km/s', 't_mix':'yr'}
    vmin_dict = {'redshift':0.9, 'log_mass':8.0, 'log_sfr':-4.0, 'Zgrad_kpc':-2, 'Zgrad_re':-2, 'logZ_logSFR_slope':-10, 'vdisp_mean':None, 'vdisp_50':None, 't_mix':0}
    vmax_dict = {'redshift':1.9, 'log_mass':11.0, 'log_sfr':1.0, 'Zgrad_kpc':2, 'Zgrad_re':2, 'logZ_logSFR_slope':10, 'vdisp_mean':None, 'vdisp_50':None, 't_mix':300}

    p = ax.scatter(df[xcol], df[ycol], s=size, lw=1, ec='k', c=df[colorcol] if colorcol is not None else color)

    xerr = [df['vdisp_50'] - df['vdisp_16'], df['vdisp_84'] - df['vdisp_50']] if xcol == 'vdisp_50' else df[xcol + '_u'] if xcol + '_u' in df else None
    yerr = [df['vdisp_50'] - df['vdisp_16'], df['vdisp_84'] - df['vdisp_50']] if ycol == 'vdisp_50' else df[ycol + '_u'] if ycol + '_u' in df else None
    ax.errorbar(df[xcol], df[ycol], yerr=yerr, xerr=xerr, fmt='none', lw=1, color=color, alpha=0.7, zorder=-5)

    ax.set_xlim(vmin_dict[xcol], vmax_dict[xcol])
    ax.set_ylim(vmin_dict[ycol], vmax_dict[ycol])
    ax = annotate_axes(ax, f'{label_dict[xcol]} ({unit_dict[xcol]})', f'{label_dict[ycol]} ({unit_dict[ycol]})', args=args, hide_cbar=colorcol is None, p=p, clabel=label_dict[colorcol] if colorcol is not None else '')

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
    Z_SFR_slope_file = args.output_dir / 'catalogs' / f'Z_{args.Zdiag}_SFR_slopes{tie_vdisp_text}.csv'

    # -----------read dataframe-------------------
    df = pd.read_csv(Z_SFR_slope_file)

    # -------make plots-----------
    fig, axes = plt.subplots(2, 2, figsize=(10, 6))
    fig.subplots_adjust(left=0.08, right=0.92, bottom=0.12, top=0.98, wspace=0.5, hspace=0.3)
    axes = axes.flatten()

    axes[0] = make_plot('log_mass', 'logZ_logSFR_slope', df, axes[0], args, colorcol='redshift')
    axes[1] = make_plot('log_mass', 't_mix', df, axes[1], args, colorcol='redshift')
    axes[2] = make_plot('vdisp_mean', 't_mix', df, axes[2], args, colorcol='log_mass')
    axes[3] = make_plot('vdisp_50', 't_mix', df, axes[3], args, colorcol='log_mass')

    # ----------save figure---------------
    figname = f'all_relations_Z_{args.Zdiag}_{tie_vdisp_text}.png'
    save_fig(fig, args.fig_dir, figname, args)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
