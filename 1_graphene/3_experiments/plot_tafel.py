

import numpy as np 
import xlrd
import click 
import os
from pathlib import Path
import matplotlib
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"
plt.rc('axes', labelsize=32)    # fontsize of the x and y labels
plt.rcParams['xtick.labelsize'] = 28
plt.rcParams['ytick.labelsize'] = 28
from useful_functions import get_fit_from_points
from matplotlib import cm


@click.command()
@click.option('--f')
def main(f):

    output = 'output'
    Path(output).mkdir(parents=True, exist_ok=True)

    loc = f

    fig, ax = plt.subplots(1, 1, figsize=(6, 8))

    wb = xlrd.open_workbook(loc)
    fit_lim = -1.

    sheet_names = wb.sheet_names()
    cmap = matplotlib.cm.get_cmap('viridis', len(sheet_names))

    for i, label in enumerate(sheet_names):
        sheet = wb.sheet_by_index(i)
        pot_sheet = sheet.col_slice(0, start_rowx=2, end_rowx=None)
        current_sheet = sheet.col_slice(-1, start_rowx=2, end_rowx=None)

        potential = []
        current = []

        for j in range(len(pot_sheet)):
            try:
                potential.append(float(pot_sheet[j].value))
            except ValueError:
                continue
            current.append(float(current_sheet[j].value))
        
        ax.plot(potential, np.log10(current), 'o', label=label, color=cmap(i))
        potential = np.array(potential)
        current = np.array(current)

        # Fit the Tafel slope
        fit_index = [i for i in range(len(potential)) if potential[i] > fit_lim ]
        fit_potential = potential[fit_index]
        fit_log_current = np.log10(current[fit_index])

        fit_tafel = get_fit_from_points(fit_log_current, fit_potential, 1) 
        fit_plot = get_fit_from_points(fit_potential, fit_log_current, 1)
        ax.plot(potential, fit_plot['p'](potential), color=cmap(i), alpha=0.5)
        print(label, -1 * fit_tafel['fit'][0]* 1000)


    # ax.set_yscale('log')
    ax.set_ylabel(r'$\text{log}(j_{CO})$ / $mAcm^{-2}$')
    ax.set_xlabel(r'$U_{NHE}$ / V')
    ax.axvline(fit_lim, color='k', alpha=0.5, ls='--')

    ax.tick_params(direction='out', length=6, width=2, colors='k',
                grid_color='r', grid_alpha=0.5)

    ax.legend(loc='best')
    fig.tight_layout()
    fig.savefig(os.path.join(output, 'tafel_plot_'+str(fit_lim).replace('-', 'm')+'.pdf')) 

        




if __name__ == "__main__":
    main()
