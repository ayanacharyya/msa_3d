'''
    Filename: plot_relations.py
    Notes: Performs metallicity and SFR maps from MSA-3D emission line maps and stores them in fits files
    Author : Ayan
    Created: 30-03-26
    Example: run plot_relations.py --Zdiag NB
'''

from header import *
from util import *

# --------------------------------------------------------------------------------------------------------------------
def make_plot(xcolname, ycolname, df, ax, args, color='cornflowerblue', size=50):
    '''
    Accepts dataframe, and column names to plot, and plots those on a given axis handle, then annotates the axes based on label dictionaries
    Returns the axis
    '''
    label_dict = {'log_mass':'Log stellar mass', 'log_sfr':'Log SFR', 'Zgrad_kpc':r'$\nabla_r$Z', 'Zgrad_re':r'$\nabla_r$Z', 'logZ_logSFR_slope':r'$\nabla_r$Z - $\Sigma_{\rm SFR}$ slope', 'vdisp_mean':r'Mean $\sigma_{\rm disp}$', 'vdisp_50':r'Median $\sigma_{\rm disp}$', 't_mix':r'$t_{\rm mix}$'}
    unit_dict = {'log_mass':r'M$_{\odot}$', 'log_sfr':r'M$_{\odot}$/yr', 'Zgrad_kpc':'dex/kpc', 'Zgrad_re':r'dex/r$_e$', 'logZ_logSFR_slope':'', 'vdisp_mean':'km/s', 'vdisp_50':'km/s', 't_mix':'yr'}
    vmin_dict = {'log_mass':8.0, 'log_sfr':-4.0, 'Zgrad_kpc':-2, 'Zgrad_re':-2, 'logZ_logSFR_slope':-10, 'vdisp_mean':0, 'vdisp_50':0, 't_mix':0}
    vmax_dict = {'log_mass':11.0, 'log_sfr':1.0, 'Zgrad_kpc':2, 'Zgrad_re':2, 'logZ_logSFR_slope':10, 'vdisp_mean':100, 'vdisp_50':100, 't_mix':300}

    ax.scatter(df[xcolname], df[ycolname], s=size, lw=1, ec='k', c=color)

    xerr = [df['vdisp_50'] - df['vdisp_16'], df['vdisp_84'] - df['vdisp_50']] if xcolname == 'vdisp_50' else df[xcolname + '_u'] if xcolname + '_u' in df else None
    yerr = [df['vdisp_50'] - df['vdisp_16'], df['vdisp_84'] - df['vdisp_50']] if ycolname == 'vdisp_50' else df[ycolname + '_u'] if ycolname + '_u' in df else None
    ax.errorbar(df[xcolname], df[ycolname], yerr=yerr, xerr=xerr, fmt='none', lw=1, color=color, alpha=0.7)

    ax.set_xlim(vmin_dict[xcolname], vmax_dict[xcolname])
    ax.set_ylim(vmin_dict[ycolname], vmax_dict[ycolname])
    ax = annotate_axes(ax, f'{label_dict[xcolname]} ({unit_dict[xcolname]})', f'{label_dict[ycolname]} ({unit_dict[ycolname]})', args=args)

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
    fig, axes = plt.subplots(2, 2, figsize=(8, 6), layout='constrained')
    axes = axes.flatten()

    axes[0] = make_plot('logZ_logSFR_slope', 'log_mass', df, axes[0], args)
    axes[1] = make_plot('t_mix', 'log_mass', df, axes[1], args)
    axes[2] = make_plot('vdisp_mean', 't_mix', df, axes[2], args)
    axes[3] = make_plot('vdisp_50', 't_mix', df, axes[3], args)

    # ----------save figure---------------
    figname = f'all_relations_Z_{args.Zdiag}_{tie_vdisp_text}.png'
    save_fig(fig, args.fig_dir, figname, args)
