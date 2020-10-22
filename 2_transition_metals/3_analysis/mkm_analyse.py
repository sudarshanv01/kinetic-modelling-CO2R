
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
import json

import matplotlib.pyplot as plt
plt.style.use('science')


def plot_rate_maps(maps, descriptors, points, potential, pH, min_val=1e-20):

    fig, ax = plt.subplots(1, 1, figsize=(8,6))

    x = [] ; y = [] ; z = [] 
    for row in maps:
        descriptor, production = row
        dGa, dGb = descriptor
        rate_CO = production[-1]
        x.append(dGa) ; y.append(dGb)
        z.append(rate_CO)

    z = np.array(z).transpose()
    x = np.array(x) ; y = np.array(y)

    z1 = [min_val if a_ < min_val else a_ for a_ in z]
    tcf = ax.tricontourf(x, y, z1, cmap='viridis',  locator=ticker.LogLocator())

    ax.set_xlabel(r'$\Delta E_{'+descriptors[0]+'}$ / eV')
    ax.set_ylabel(r'$\Delta E_{'+descriptors[1]+'}$ / eV')
    cb = fig.colorbar(tcf, ax=ax)

    all_Ga = []
    all_Gb = []
    for metal in points:
        color = jmol_colors[atomic_numbers[metal]]
        ax.plot(points[metal][descriptors[0]], points[metal][descriptors[1]], 'o', color=color)
        ax.annotate(metal, xy=(points[metal][descriptors[0]], points[metal][descriptors[1]]))
        all_Ga.append(points[metal][descriptors[0]])
        all_Gb.append(points[metal][descriptors[1]])
    
    fit = get_fit_from_points(all_Ga, all_Gb, 1)
    ax.annotate('slope = %1.2f'%fit['fit'][0], color='white', xy=(0.1, 0.8), xycoords='axes fraction')
    bounds = ax.get_xbound()
    ax.plot(bounds, fit['p'](bounds), color='k', ls='-', lw=3)

    return fig, ax

def plot_points(energy_file, surfaces):
    data = [a.split('\t') for a in energy_file.get_content().split('\n') ]
    results = {}
    for dat in data:
        if dat[0] in surfaces:
            # This is a point on the graph 
            dE = float(dat[3])
            results.setdefault(dat[0],{})[dat[2]+'_s'] = dE
    return results

@click.command()
@click.option('--pk')
def main(pk):
    node = load_node(pk)
    descriptors = node.inputs.descriptor_names.get_list()
    potential = node.inputs.voltage
    pH = node.inputs.pH
    surfaces = node.inputs.surface_names.get_list()
    energy_file = node.inputs.energies
    points = plot_points(energy_file, surfaces)

    # rate_map = node.outputs.rate_map
    # reactions = node.inputs.rxn_expressions.get_list()
    # figR, axR = plot_maps(rate_map, descriptors, points, potential=float(potential), type='rates')
    # figR.savefig(os.path.join(output, 'rate_map_'+str(pk)+'.pdf'))

    production_rate = node.outputs.production_rate_map.get_list()
    print(points)
    fig, ax = plot_rate_maps(production_rate, descriptors, points, potential, pH )
    fig.savefig(os.path.join(output, 'rate_map_'+str(pk)+'.pdf'))


if __name__ == '__main__':
    output = 'output'
    Path(output).mkdir( exist_ok=True, parents=True)

    main()

