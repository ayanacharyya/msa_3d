'''
    Filename: plot_ohno_oceans.py
    Notes: Compares OHNO and BPT diagrams using 2D NIRSpec data from OCEANS data
    Author : Ayan
    Created: 04-06-26
    Example: run plot_ohno_oceans.py --snr_cut 3
'''

from header import *
from util import *
start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    args.fontfactor = 1.2
    
    # -------------setup directories and global variables----------------
    data_dir = args.input_dir / 'oceans/kinematics'
    speccat_file = args.output_dir / 'catalogs/oceans_fluxes.csv'
    fig_dir = args.output_dir / 'plots'
    snr_cut_text = f'_snr{args.snr_cut}' if args.snr_cut is not None else ''

    # ----------------reading in spectrosocpic catalog---------------------
    if not os.path.exists(speccat_file) or args.clobber:
        # ----------------reading in input catalog---------------------
        files=glob.glob(str(data_dir) + '/*kinfit*.pkl')

        lines = ['Hbeta', '[OIII]4960', '[OIII]5008', 'Halpha']
        df_spec = pd.DataFrame(columns=np.hstack([['id', 'redshift'], lines]))

        # -------------setup the figure-----------------
        nrow, ncol = 1, 2
        fig, axes = plt.subplots(nrow, ncol, figsize=(7, 5.5))
        fig.subplots_adjust(left=0.08, right=0.93, top=0.93, bottom=0.1, wspace=0., hspace=0.1)
        
        # ----------------looping over the objects-------------
        for index, file in enumerate(files):
            print(f'\\nCommencing ({index + 1}/{len(files)})..')
            this_dict = {}

            # -------reading in file-------------
            f = open(file, 'rb')
            fit_results = pickle.load(f)['fit_result']

            # -------reading in basic parameters-------------
            id = fit_results['oceans_id']
            redshift = fit_results['z']
            available_lines = list(fit_results['line_fits'].keys())
            print(f'\tObject {id} has lines: {available_lines}')

            this_dict.update({'id': id,
                              'redshift': redshift,
                              })

            # ------reading in 2D spectra---------
            for line in lines:
                print(f'\n\tDoing {line} for ID {id}..')

                if line in available_lines:
                    line_dict = fit_results['line_fits'][line]
                    rest_wave = line_dict['obs_wave'] / (1 + redshift) # micron

                    # -----------summing 2D spectra to 1D--------------
                    line_flux = 0
                    for index, this_row in enumerate(line_dict['row_fits']):
                        if 'map_1g' in this_row:
                            print(f'\t\tExtracting fluxes from row {index+1}/{len(line_dict["row_fits"])}..')
                            
                            # -----------deriving line fluxes--------------
                            amplitude_row = this_row['map_1g'][0] # MJy
                            sigma_row = this_row['map_1g'][2] # km/s
                            flux_row = np.sqrt(2 * np.pi) * 1e-8 * (amplitude_row * sigma_row) / rest_wave # ergs/s/cm^2
                            line_flux += flux_row
                        else:
                            print(f'\t\tRow {index+1}/{len(line_dict["row_fits"])} does not have map_1g component, so skipping this')
                    if line_flux == 0:
                        print(f'\t1-component fits {line} unavailable for ID {id}')
                        this_dict.update({line: np.nan})
                    else:
                        print(f'\tIntegrated {line} flux for ID {id} is {line_flux} ergs/s//cm^2')
                        this_dict.update({line: line_flux})

                else:
                    print(f'\t{line} unavailable for ID {id}')
                    this_dict.update({line: np.nan})

            # ----------appneding to speccat--------------------
            df_spec.loc[len(df_spec)] = this_dict

            # ----------memory cleanup--------------------
            del f
            del fit_results
            gc.collect()
        
        # -----------writing out speccat---------------
        df_spec.to_csv(speccat_file, index=None)
        print(f'Saved speccat for {len(df_spec)} objects as {speccat_file}')
    
    else:
        print(f'Reading from existing {speccat_file}')

    # -------------making the plot-------------------
    #axes[0][0] = plot_bpt(df_spec, axes[0][0], args, ynum='OIII', yden='Hb', xnum='NeIII', xden='OII', hide_xaxis=False, hide_yaxis=False, hide_cbar=True)            
    #axes[0][1] = plot_bpt(df_spec, axes[0][1], args, ynum='OIII', yden='Hb', xnum='NII', xden='Ha', hide_xaxis=False, hide_yaxis=True, hide_cbar=True)            

    fig.text(0.1, 0.98, f'#galaxies = {len(df_spec)}', fontsize=args.fontsize, c='k', ha='left', va='top')

    # ---------saving figures-------------------
    figname = f'oceans_compare_ohno{snr_cut_text}.png'
    save_fig(fig, fig_dir, figname, args)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')