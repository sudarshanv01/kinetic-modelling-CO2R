import click
from pprint import pprint
import json
from useful_functions import create_output_directory
import numpy as np
from ase.data import atomic_numbers
from pprint import pprint
from useful_functions import get_fit_from_points
from scipy.interpolate import griddata, interp2d, Rbf
import matplotlib.pyplot as plt
from ase.data.colors import jmol_colors
from matplotlib import ticker
from plot_params import get_plot_params
import xlrd
import string


def main(ax, low):


    sheet_names_low = low.sheet_names()
    potential_sheet = low.sheet_by_index(0) #sheet_names_low[0]
    current_sheet = low.sheet_by_index(1)#sheet_names_low[1]
    axt = ax.twinx()

    time = potential_sheet.col_slice(0, start_rowx=2, end_rowx=None)
    rhe_potential = potential_sheet.col_slice(1, start_rowx=2, end_rowx=None)
    time = [float(a.value) for a in time]
    rhe_potential = [float(a.value) for a in rhe_potential]
    axt.plot(time, rhe_potential, 'o-', color='k', alpha=0.25, label='Potential')

    time_c = current_sheet.col_slice(0, start_rowx=2, end_rowx=None)
    Fe_co = current_sheet.col_slice(1, start_rowx=2, end_rowx=None)
    Ni_co = current_sheet.col_slice(2, start_rowx=2, end_rowx=None) 
    time_c = [float(a.value) for a in time_c]
    Fe_co = [float(a.value) for a in Fe_co]
    Ni_co = [float(a.value) for a in Ni_co]
    
    ax.plot(time_c, Ni_co, 'o-', color='tab:blue')
    ax.plot(time_c, Fe_co, 'o-', color='tab:red')
    ax.set_xlabel(r'Time / min')
    axt.set_ylabel(r'Potential / V vs. RHE')


if __name__ == '__main__':
    get_plot_params()

    low = xlrd.open_workbook('experiments/degradation_small.xlsx')
    high = xlrd.open_workbook('experiments/degradation_large.xlsx')
    fig, ax = plt.subplots(2, 1, figsize=(10,10), constrained_layout=True)
    main(ax[0],low)
    main(ax[1],high)
    ax[0].plot([],[], color='tab:blue', label='NiNC')
    ax[0].plot([],[], color='tab:red', label='FeNC')
    ax[0].legend(loc='best', frameon=False, fontsize=14)
    ax[0].set_ylabel(r'$j_{\mathregular{CO}}$ / mAcm$^{-2}$')
    ax[1].plot([],[], color='tab:red', label=r'H$_{2}$')
    ax[1].plot([],[], color='tab:blue', label=r'CO')
    ax[1].set_ylabel(r'$j$ / mAcm$^{-2}$')
    ax[1].legend(loc='best', frameon=False, fontsize=14)
    alphabet = list(string.ascii_lowercase)
    for i, a in enumerate(ax):
        a.annotate(alphabet[i]+')', xy=(0.0, 1.1), xycoords='axes fraction', fontsize=20)
    fig.savefig('output_figure/SI_degradation_data.pdf')