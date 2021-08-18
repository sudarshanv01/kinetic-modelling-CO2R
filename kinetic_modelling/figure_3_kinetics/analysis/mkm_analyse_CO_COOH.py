
"""

For a given CalcJobNode plot out the coverage map
the rate map and the production rate map 

Needs the following inputs 

* Node of the CalcJobNode where CatMAP was the parser

"""

import sys
import os
import numpy as np 
from pathlib import Path
from useful_classes import bcolors
import matplotlib
from ase.data.colors import jmol_colors
from ase.data import atomic_numbers
from pprint import pprint
from useful_functions import get_fit_from_points
import click
from ase.thermochemistry import HarmonicThermo
from matplotlib import ticker
from useful_functions import create_output_directory
from plot_params import get_plot_params
import json
import string
from scipy.interpolate import griddata, interp2d, Rbf

import matplotlib.pyplot as plt
# plt.style.use('science')


def plot_CO2_vs_COOH(fig,ax, maps, descriptors, points, potential, pH, \
            plot_single_atom, plot_metal, min_val=1e-20, log_scale=True,\
                cmapname='Blues_r', annotate_rate_limiting=False,\
                    coverage_plot=False, plot_cmap=True, inten=1):


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
            z.append(production[-1])

    z_cov = [0.0 if 0.75 > a > 0.1 else 1. for a in z]

    z = np.array(z).transpose()
    x = np.array(x) ; y = np.array(y)

    ## Add free energy contributions 
    x = x 
    y = y 


    # Remove any rates that are very low
    if plot_cmap:
        z1 = [min_val if a_ < min_val else a_ for a_ in z]
        # tcf = ax.tricontourf(x, y, np.log10(z1), cmap='cividis',  locator=ticker.LogLocator())
        if log_scale:
            tcf = ax.tripcolor(x, y, np.log10(z1), shading='gouraud', cmap=cmapname, alpha=inten, edgecolors='k')
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
            cb.ax.set_title(r'$\theta_{\mathregular{CO}}$ / ML ')
            ax.plot([],[],color='tab:green', lw=6, label='SV')
            ax.plot([],[],color='tab:purple', lw=6, label='DV')
            ax.legend(loc='lower right', frameon=False, fontsize=14)
            # ax.tripcolor(x, y, z_cov, cmap='Greys_r', shading='gouraud', alpha=0.25)
            # cb.set_clim(0.,1.)
        else:
            cb.ax.set_title(r'TOF / $\mathregular{s}^{-1}$')

        ax.set_xlabel(r'$\Delta E_{\mathregular{'+descriptors[0].replace('_s','').replace('2', '_{2}')+'}}(%1.1f\mathregular{V}_{\mathregular{SHE}})$ / eV'%potential)
        ax.set_ylabel(r'$\Delta E_{\mathregular{'+descriptors[1].replace('_s','').replace('2', '_{2}')+'}}(%1.1f\mathregular{V}_{\mathregular{SHE}})$ / eV'%potential)


    all_Ga = []
    all_Gb = []
    for metal in points:
        ## color for points
        try:
            color = jmol_colors[atomic_numbers[metal]]
            color='k'
        except KeyError:
            color = jmol_colors[atomic_numbers[metal.split('_')[0]]]

        try:
            # TODO: Delete this later
            if points[metal][descriptors[0]] > 1.5:
                continue 
            
            # Elements to scale with
            x = points[metal][descriptors[0]] 
            y = points[metal][descriptors[1]]

            x = x 
            y = y 
            if 'Ni' in metal and plot_single_atom == False:
                continue
            if 'Fe' in metal and plot_single_atom == False:
                continue
            if 'Ni' in metal or 'Fe' in metal and plot_single_atom == True:
                label_info = metal.split('_')
                color_vac = 'tab:green' if int(label_info[1])==1 else 'tab:purple'

                ax.plot(x , y, marker='$\mathregular{%s N}_{%d}$'%(label_info[0],\
                        int(label_info[-1])), color=color_vac, markersize=22)
            elif plot_metal:
                ax.plot(x , y, 'o', color=color)

        except KeyError:
            continue
        
        if plot_single_atom:
            if 'Ni' in metal or 'Fe' in metal:

                n_dop = metal.split('_')[-1]
                met = metal.split('_')[0]

                # ax.annotate(met + r'N$_{%s}'%n_dop+'$', 
                #             xy=(x,y), color=color)
        # else:
        if plot_metal:
            if 'Ni' in metal or 'Fe' in metal:
                continue
            met = metal
            ax.annotate(metal,
                        xy=(x,y), color=color, fontsize=12)
        
        if 'Ni' not in metal and 'Fe' not in metal:
            all_Ga.append(x)
            all_Gb.append(y)
    
    fit = get_fit_from_points(all_Ga, all_Gb, 1)
    bounds = ax.get_xbound()
    if plot_metal:
        ax.plot(bounds, fit['p'](bounds), color='k', ls='--', lw=3)
        ax.fill_between(bounds, fit['p'](bounds)+0.1, fit['p'](bounds)-0.1, color='k', alpha=0.1)

    ax.set_xlim([-2,2])
    ax.set_ylim([-2,2])

    
def plot_points(energy_file, surfaces, facet):
    data = [a.split('\t') for a in energy_file.get_content().split('\n') ]
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

@click.command()
@click.option('--pk1')
def main(pk1):

    fig, ax = plt.subplots(1,1,figsize=(8,6), squeeze=False)

    node = load_node(pk1)
    descriptors = node.inputs.descriptor_names.get_list()
    print(descriptors)
    potential = node.inputs.voltage
    pH = node.inputs.pH
    surfaces = node.inputs.surface_names.get_list()
    # get the surface facet calculated 
    species_definitions = node.inputs.species_definitions
    facet = species_definitions['s']['site_names'][0]
    energy_file = node.inputs.energies
    points = plot_points(energy_file, surfaces, facet)
    production_rate = node.outputs.production_rate_map.get_list()
    coverage_map = node.outputs.coverage_map.get_list()

    plot_CO2_vs_COOH(fig, ax[0,0], production_rate, descriptors, points, potential.value, pH.value, \
        plot_single_atom=False, plot_metal=True, cmapname='coolwarm', annotate_rate_limiting=True)

    fig.tight_layout()
    fig.savefig(os.path.join('output','rate_map_'+str(pk1)+'.pdf'))



if __name__ == '__main__':
    create_output_directory()
    get_plot_params()

    main()

