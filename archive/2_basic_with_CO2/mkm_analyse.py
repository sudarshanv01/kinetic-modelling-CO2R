
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
import matplotlib
from ase.data.colors import jmol_colors
from ase.data import atomic_numbers
from pprint import pprint
from useful_functions import get_fit_from_points
from ase.thermochemistry import HarmonicThermo
from matplotlib import ticker
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"
plt.rc('axes', labelsize=32)    # fontsize of the x and y labels
plt.rcParams['xtick.labelsize'] = 28
plt.rcParams['ytick.labelsize'] = 28

def plot_maps(coverage_map, descriptors, titles, points, type='coverage'):
    a = len(coverage_map[0][1])
    fig, ax = plt.subplots(a, 1, figsize=(3*a, 5*a))
    x = [] ; y = [] ; z = [] 
    for row in coverage_map:
        x.append(row[0][0])
        y.append(row[0][1])
        z.append(row[1])
    z = np.array(z).transpose()
    x = np.array(x) ; y = np.array(y)
    for i, z1 in enumerate(z):
        if type == 'coverage':
            tcf = ax[i].tricontourf(x, y, z1, cmap='coolwarm')
        elif type == 'rates':
            try:
                tcf = ax[i].tricontourf(x, y, z1, cmap='coolwarm',  locator=ticker.LogLocator())
            except ValueError:
                continue
        ax[i].set_xlabel(r'$\Delta E_{'+descriptors[0]+'}$')
        ax[i].set_ylabel(r'$\Delta E_{'+descriptors[1]+'}$')
        # ax[i].set_title('$'+titles[i]+'$')
        cb = fig.colorbar(tcf, ax=ax[i])
    
    # Plot the points 
        dG_COOH = []
        dG_CO = []
        for surface, dG in points.items():
            try:
                metal = surface[0:2]
                colors = jmol_colors[atomic_numbers[metal]]
            except KeyError:
                metal = surface[0]
                colors = jmol_colors[atomic_numbers[metal]]

            ax[i].plot(dG['CO'], dG['COOH'], 'o', color=colors)
            dG_COOH.append(dG['COOH'])
            dG_CO.append(dG['CO'])

        fig.tight_layout()
        fit = get_fit_from_points(dG_CO, dG_COOH, 1)
        print(fit)
        ax[i].plot(ax[i].get_xbound(), fit['p'](ax[i].get_xbound()), '-', color='k', lw=3)


    return [fig, ax]    

def plot_points(energy_file, surfaces):
    data = [a.split('\t') for a in energy_file.get_content().split('\n') ]
    results = {}
    for dat in data:
        if dat[0] in surfaces:
            # This is a point on the graph 
            dE = float(dat[3])
            vibs = np.array(list(map(float,dat[4].strip('][').split(', '))))
            dG = HarmonicThermo(0.00012*vibs, dE).get_helmholtz_energy(temperature=300., verbose=False)
            results.setdefault(dat[0],{})[dat[2]] = dG
    return results


if __name__ == '__main__':
    output = 'output'
    Path(output).mkdir( exist_ok=True, parents=True)
    nodename = int(sys.argv[1])
    node = load_node(nodename)
    descriptors = node.inputs.descriptor_names.get_list()

    surfaces = node.inputs.surface_names.get_list()
    energy_file = node.inputs.energies
    points = plot_points(energy_file, surfaces)

    coverage_map = node.outputs.coverage_map
    # get the values needed for the title 
    figC, axC = plot_maps(coverage_map, descriptors, 
                        titles=descriptors, points=points)
    figC.savefig(os.path.join(output, 'coverage_map.pdf'))

    rate_map = node.outputs.rate_map
    reactions = node.inputs.rxn_expressions.get_list()
    figR, axR = plot_maps(rate_map, descriptors, reactions, points, type='rates')
    figR.savefig(os.path.join(output, 'rate_map.pdf'))

