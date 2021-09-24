
"""Create Figure 4 of the manuscript."""

import click
import json
import string
from pathlib import Path
import xlrd
import numpy as np
from pprint import pprint
from scipy.interpolate import Rbf
from matplotlib import ticker
import matplotlib.pyplot as plt
from ase.data.colors import jmol_colors
from ase.data import atomic_numbers
from ase import build
from ase.thermochemistry import IdealGasThermo, HarmonicThermo
from useful_functions import get_fit_from_points
from plot_params import get_plot_params

def get_electronic_energy_points(energy_file, surfaces, facet):
    """Get the electronic energy for all the points from the CatMAP input file."""
    data = [a.split('\t') for a in energy_file.split('\n') ]
    results = {}
    for i, dat in enumerate(data):
        if i == 0:
            continue
        if dat[0] in surfaces and dat[1] == facet:
            dE = float(dat[3])
            results.setdefault(dat[0],{})[dat[2]+'_s'] = dE
        if dat[0] in ['Ni', 'Fe']:
            dE = float(dat[3])
            results.setdefault(dat[0]+'_'+dat[1],{})[dat[2]+'_s'] = dE
    return results

class FreeEnergiesForFigure4:
    """Get the free energies to plot in Figure 4."""
    def __init__(self):
        self.frequencies = {}
        self.frequencies['CO'] = [1832.373539, 481.555294, 467.482512, 425.061714, 73.09318, 68.233697]
        self.frequencies['COOH'] = [3579.438347, 1573.100374, 1241.519251, 932.338544, 624.650858, 619.94479, 468.883838, 246.531351, 245.478053, 75.740879, 73.602391, 15.378853]
        self.frequencies['CO2'] = [1846.027301, 1179.574916, 546.070088, 534.269426, 168.712357, 136.301896, 66.113421, 51.400916, 36.766228]

        self.frequencies['COg'] = [ 30.574661, 30.574661, 2102.17 ]
        self.frequencies['CO2g'] =  [ 8.210836, 8.986016, 626.201759, 626.207647, 1303.431414, 2337.103721 ]
        self.frequencies['H2g'] = [ 101.851377, 101.851377, 4357.748816 ]
        self.frequencies['H2Og'] = [ 10, 76.593339, 84.905187, 1599.342801, 3715.405649, 3823.990900 ] 
    
    def get_free_energies(self):
        """Get the free energies for the different species at standard conditions."""
        frequencies = self.frequencies

        # Pressure to plot the free energy diagram
        pressure = {'CO2g': 0.2, 'COg': 0.1, 'H2g': 1., 'H2Og': 1.}
        print('Pressure of molecules (in atm):')
        pprint(pressure)

        cmtoeV = 0.00012 
        COg_G = IdealGasThermo(
            cmtoeV * np.array(frequencies['COg']),
            geometry='linear', 
            atoms=build.molecule('CO'),
            symmetrynumber=1,
            spin=0,
        ).get_gibbs_energy(298.15, pressure['COg']*101325, verbose=False)

        CO2g_G = IdealGasThermo(
            cmtoeV * np.array(frequencies['CO2g']),
            geometry='linear', 
            atoms=build.molecule('CO2'),
            symmetrynumber=2,
            spin=0,
        ).get_gibbs_energy(298.15, pressure['CO2g']*101325, verbose=False) 

        H2g_G = IdealGasThermo(
            cmtoeV * np.array(frequencies['H2g']),
            geometry='linear', 
            atoms=build.molecule('H2'),
            symmetrynumber=2,
            spin=0,
        ).get_gibbs_energy(298.15, pressure['H2g']*101325, verbose=False)

        H2Og_G = IdealGasThermo(
            cmtoeV * np.array(frequencies['H2Og']),
            geometry='nonlinear', 
            atoms=build.molecule('H2O'),
            symmetrynumber=2,
            spin=0,
        ).get_gibbs_energy(298.15, pressure['H2Og']*101325, verbose=False)

        COOH_ads_G = HarmonicThermo(cmtoeV * np.array(self.frequencies['COOH'])).get_helmholtz_energy(298.15, verbose=False)
        CO2_ads_G = HarmonicThermo(cmtoeV * np.array(self.frequencies['CO2'])).get_helmholtz_energy(298.15, verbose=False)

        dG_CO2 = CO2_ads_G - CO2g_G
        dG_COOH = COOH_ads_G - (CO2g_G + 0.5 * H2g_G)

        print('Standard condition free energy corrections of CO2 and COOH are')
        print(f'dG(CO2)     :{dG_CO2}')
        print(f'dG(COOH)    :{dG_COOH}')
        print('\n')
        return dG_CO2, dG_COOH

def plot_map(fig, ax, maps, descriptors, points, potential, pH, \
            plot_single_atom, plot_metal, coverage_index=-1, min_val=1e-20, log_scale=True,\
                cmapname='Blues_r', annotate_rate_limiting=False,\
                    coverage_plot=False, plot_cmap=True, inten=1, plot_legend=False):
    """Generate the main kinetic plot; this function keeps getting called for each axis

    :param fig: matplotlib figure
    :type fig: fig
    :param ax: axis to plot on
    :type ax: ax
    :param maps: map to plot on axis
    :type maps: list
    :param descriptors: descriptor names
    :type descriptors: list
    :param points: data points to be plotted on graph
    :type points: dict
    :param potential: SHE potential
    :type potential: float
    :param pH: pH that was used in the run
    :type pH: float
    :param plot_metal: Should the metal points be plotted?
    :type plot_metal: bool
    :param coverage_index: index of the species used to plot the coverage - same indexing as CatMAP, defaults to -1
    :type coverage_index: int, optional
    :param min_val: minimum value to be used, defaults to 1e-20
    :type min_val: float, optional
    :param log_scale: Should a log scale be used?, defaults to True
    :type log_scale: bool, optional
    :param annotate_rate_limiting: should annotations be used?, defaults to False
    :type annotate_rate_limiting: bool, optional
    :param plot_cmap: Plot the colormap?, defaults to True
    :type plot_cmap: bool, optional
    :param inten: Intensity of plot - alpha setting in matplotlib, defaults to 1
    :type inten: int, optional
    """

    ## CHE correction for COOH*
    species = [r'CO_{2}', 'COOH', 'CO']
    CHE_COOH = potential + 0.059 * pH

    x = [] ; y = [] ; z = [] 
    z_cov = [] ## for the coverage plot
    for row in maps:
        descriptor, production = row
        dGa, dGb = descriptor
        x.append(dGa) ; y.append(dGb)
        if not coverage_plot:
            rate_CO = np.max(production)
            z.append(rate_CO)
        if coverage_plot:
            z.append(production[coverage_index])

    x = [a + CHE_COOH for a in x]
    z = np.array(z).transpose()
    x = np.array(x) ; y = np.array(y)

    free_energies = FreeEnergiesForFigure4() 
    dG_CO2, dG_COOH = free_energies.get_free_energies()

    ## Add free energy contributions 
    x = x + dG_COOH 
    y = y + dG_CO2

    # Remove any rates that are very low
    if plot_cmap:
        z1 = [min_val if a_ < min_val else a_ for a_ in z]
        if log_scale:
            x_dense_val = np.linspace(min(x), max(x))
            y_dense_val = np.linspace(min(y), max(y))
            x_dense, y_dense = np.meshgrid(x_dense_val, y_dense_val)
            z_smooth = Rbf(x, y, np.log10(z1), function='linear', smooth=1)
            z_smooth_dense = z_smooth(x_dense, y_dense)
            tcf = ax.contourf(x_dense, y_dense, 10**z_smooth_dense, cmap=cmapname, locator=ticker.LogLocator())
        else:
            x_dense_val = np.linspace(min(x), max(x))
            y_dense_val = np.linspace(min(y), max(y))
            x_dense, y_dense = np.meshgrid(x_dense_val, y_dense_val)
            z_smooth = Rbf(x, y, z, function='linear', smooth=1)
            z_smooth_dense = z_smooth(x_dense, y_dense)
            z_smooth_dense = z_smooth_dense.clip(min=0, max=1.)
            tcf = ax.contourf(x_dense, y_dense, z_smooth_dense, cmap=cmapname, levels=np.arange(0,1.1,0.1) )

        cb = fig.colorbar(tcf, ax=ax)
        if coverage_plot:
            cb.ax.set_title(r'$\theta_{\mathregular{%s}}$ / ML '%species[coverage_index])
        else:
            cb.ax.set_title(r'TOF / $\mathregular{s}^{-1}$')
        ax.plot([],[], 'o',  color='r', label='Fe(SV)')
        ax.plot([],[], 'o', color='g',  label='Fe(DV)')
        ax.plot([],[], 'o',  color='b', label='Ni(SV)')
        ax.plot([],[], 'o', color='tab:brown',  label='Ni(DV)')

        if plot_legend:
            ax.legend(loc='lower right', frameon=False, fontsize=16)
        ax.set_xlabel(r'$\Delta G_{\mathregular{'+descriptors[0].replace('_s','').replace('2', '_{2}')+'}}(%1.1f\mathregular{V}_{\mathregular{SHE}})$ / eV'%potential)
        ax.set_ylabel(r'$\Delta G_{\mathregular{'+descriptors[1].replace('_s','').replace('2', '_{2}')+'}}(%1.1f\mathregular{V}_{\mathregular{SHE}})$ / eV'%potential)
    
    xbounds = ax.get_xbound()
    ybounds = ax.get_ybound()

    all_Ga = []
    all_Gb = []
    for metal in points:
        try:
            color = jmol_colors[atomic_numbers[metal]]
            color='k'
        except KeyError:
            color = jmol_colors[atomic_numbers[metal.split('_')[0]]]
        try:
            if points[metal][descriptors[0]] > 1.5:
                continue 
            # Elements to scale with
            x = points[metal][descriptors[0]] + CHE_COOH
            y = points[metal][descriptors[1]]

            x = x + dG_COOH
            y = y + dG_CO2

            if 'Ni' in metal and plot_single_atom == False:
                continue
            if 'Fe' in metal and plot_single_atom == False:
                continue
            if 'Ni' in metal or 'Fe' in metal and plot_single_atom == True:
                label_info = metal.split('_')
                if int(label_info[1])==1 and 'Ni' in metal: color_vac='b' 
                elif int(label_info[1])==2 and 'Ni' in metal: color_vac='tab:brown'
                elif int(label_info[1])==1 and 'Fe' in metal: color_vac='r'
                elif int(label_info[1])==2 and 'Fe' in metal: color_vac='g'
                ax.plot(x , y, marker='o', markersize=12, color=color_vac,)
            elif plot_metal:
                ax.plot(x , y, 'o', markersize=12, color=color)

        except KeyError:
            continue
        
        if plot_metal:
            if 'Ni' in metal or 'Fe' in metal:
                # continue
                met, vac, dop = metal.split('_')
                ax.annotate('$\mathregular{%sN_{%s}}$'%(met,dop), fontsize=12,color='k', xy=(x+0.05,y+0.05), alpha=0.5).draggable()
            else:
                met = metal
                ax.annotate(metal,
                            xy=(x+0.05,y+0.05), color=color, fontsize=14).draggable()
        
        if 'Ni' not in metal and 'Fe' not in metal:
            all_Ga.append(x)
            all_Gb.append(y)
    
    fit = get_fit_from_points(all_Ga, all_Gb, 1)
    bounds = ax.get_xbound()

    if plot_metal:
        ax.plot(bounds, fit['p'](bounds), color='k', ls='--', lw=3)
        ax.fill_between(bounds, fit['p'](bounds)+0.1, fit['p'](bounds)-0.1, color='k', alpha=0.1)

    ax.plot(bounds, bounds, color='k', ls='-')
    ax.fill_between(bounds, np.array(bounds)+0.1, np.array(bounds)-0.1, color='k', alpha=0.1)

    
    if annotate_rate_limiting:
        ax.arrow(-1.4, 0.2, 0, 1, width=0.05,color='white'  )
        ax.arrow(0, -1., 1, 0, width=0.05,color='white'  )
        ax.annotate(r'$\mathregular{CO_2}^* \to \mathregular{COOH}^* $ limited', xy=(0.2, 0.05), xycoords='axes fraction', color='white').draggable()
        ax.annotate(r'$ \mathregular{CO_2}_{(\mathregular{g})} \to \mathregular{CO_2}^*$ limited', xy=(0.1, 0.8), xycoords='axes fraction', color='white').draggable()
        ax.annotate('Parity Line', xy=(0.4, 0.2), xycoords='axes fraction', rotation=37, color='k')
        ax.annotate('(211) Scaling', xy=(0.1, 0.4), xycoords='axes fraction', rotation=37, color='k')
    ax.set_xlim(xbounds)
    ax.set_ylim(ybounds)

@click.command()
@click.option('--kfiles', type=str, default='aiida_output/kinetic_model_data.json')
def main(kfiles):

    kineticspk = json.load(open('chosen_pk.json'))['kinetics_pk']

    with open(kfiles, 'r') as handle:
        data_tot = json.load(handle)
    data = data_tot[kineticspk]

    fig = plt.figure(constrained_layout=True, figsize=(8,9))
    figs = plt.figure(constrained_layout=True, figsize=(10,9))
    figt = plt.figure(constrained_layout=True, figsize=(8,6))

    gs = fig.add_gridspec(6,1)
    gsS = figs.add_gridspec(2,2)
    gsT = figt.add_gridspec(1,1)
    axc = [ fig.add_subplot(gs[0:3,0]), fig.add_subplot(gs[3:,0]) ]
    axt = figt.add_subplot(gsT[0,0]) 
    axco2 = figs.add_subplot(gsS[1,:]) 
    ax = axc

    # Experiments
    tpd = xlrd.open_workbook('experiments/TPD.xls')

    # Plot the TPD graph
    colors = ['tab:red', 'tab:blue', 'tab:green']
    sheet_names_tpd = tpd.sheet_names()
    for i, label in enumerate(sheet_names_tpd):
        sheet = tpd.sheet_by_index(i)
        temp_sheet = sheet.col_slice(0, start_rowx=2, end_rowx=None)
        signal_sheet = sheet.col_slice(-1, start_rowx=2, end_rowx=None)
        temperature = []
        signal = []
        for j in range(len(temp_sheet)):
            temperature.append(float(temp_sheet[j].value)+273.15)
            signal.append(float(signal_sheet[j].value))
        index = [a for a in range(len(temperature)) if temperature[a] > 300]
        temperature = np.array(temperature)
        signal = np.array(signal)
        axt.plot(temperature[index], signal[index], '-', lw=4, color=colors[i], label=label)
    axt.set_ylabel(r'Signal / arb. units')
    axt.set_yticks([])
    axt.set_xlabel(r'Temperature / K')
    axt.legend(loc='best', frameon=False)

    ## computational plot
    data_points = get_electronic_energy_points(data['energy_file'], data['surfaces'], data['facet'][0])

    for a in axc:
        a.set_xlim([-1.75, 1.6])
        a.set_ylim([-1.75, 1.6])

    plot_map(
        fig=fig,
        ax=axc[0],
        maps=data['production_rate'],
        descriptors=data['descriptors'],
        points=data_points,
        potential=data['potential'],
        pH=data['pH'],
        plot_single_atom=True,
        cmapname='coolwarm',
        plot_metal=True,
        annotate_rate_limiting=True,
    )
    plot_map(
        fig=fig,
        ax=axc[1],
        maps=data['coverage_map'],
        descriptors=data['descriptors'],
        points=data_points,
        potential=data['potential'],
        pH=data['pH'],
        plot_single_atom=True,
        cmapname='Oranges',
        plot_metal=True,
        log_scale=False,
        coverage_plot=True,
        coverage_index=-1,
        plot_legend=True,
    )
    axc[1].annotate(r'CO$^*$ poisoned', xy=(0.03, 0.9), color='white', xycoords='axes fraction')

    ## Label the diagram
    alphabet = list(string.ascii_lowercase)
    for i, a in enumerate(ax):
        a.annotate(alphabet[i]+')', xy=(0.0, 1.1), xycoords='axes fraction', fontsize=20).draggable()
    plt.show()
    fig.savefig('output_figure/figure_kinetics.pdf')
    figt.savefig('output_figure/TPD.pdf')

if __name__ == '__main__':
    Path('./output_figure').mkdir(parents=True, exist_ok=True)
    get_plot_params()
    main()
    
