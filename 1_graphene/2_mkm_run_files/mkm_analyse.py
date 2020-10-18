
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

matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"
plt.rc('axes', labelsize=32)    # fontsize of the x and y labels
plt.rcParams['xtick.labelsize'] = 28
plt.rcParams['ytick.labelsize'] = 28

def plot_maps(coverage_map, descriptors,  points, potential, type='coverage', min_val=1e-20):
    """Generate the 2D heat map from a given CatMAP run

    Args:
        coverage_map (list): Any 2D map 
        descriptors (list): list of len == 2
        points (list): metal points 
        type (str, optional): Type of plot that the data corresponds to. Defaults to 'coverage'.

    Returns:
        fig : matplotlib object of the figure 
        ax : matplotlib axis object
    """

    a = len(coverage_map[0][1])
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    x = [] ; y = [] ; z = [] 
    for row in coverage_map:
        x.append(row[0][0])
        y.append(row[0][1])
        z.append(row[1])
    z = np.array(z).transpose()
    x = np.array(x) ; y = np.array(y)

    z1 = z[1]
    if type == 'coverage':
        tcf = ax.tricontourf(x, y, z1, cmap='viridis')
    elif type == 'rates':
        z1 = [min_val if a_ < min_val else a_ for a_ in z1]
        try:
            tcf = ax.tricontourf(x, y, z1, cmap='viridis',  locator=ticker.LogLocator())
        except ValueError:
            print('')

    ax.set_xlabel(r'$\Delta E_{'+descriptors[0]+'}$')
    ax.set_ylabel(r'$\Delta E_{'+descriptors[1]+'}$')
    cb = fig.colorbar(tcf, ax=ax)
    
    # Plot the points 
    dG_1 = []
    dG_2 = []
    metal_list = []
    # Different possible lists 
    lists = [1]
    description_list = {
        1: 'Scaling of low-spin elements', 
        2: 'Scaling of high-spin elements',
    }
    lines_to_plot = {}
    criteria_surface = {}
    # initialise list
    for l in lists:
        lines_to_plot[l] = []
        criteria_surface[l] = []

    # get spin data
    with open('inputs/magmom.json') as handle:
        spin_data = json.load(handle)
    for surface, dG in points.items():
        try:
            metal = surface[0:2]
            colors = jmol_colors[atomic_numbers[metal]]
            non_metal = surface[2:]
        except KeyError:
            metal = surface[0]
            colors = jmol_colors[atomic_numbers[metal]]
            non_metal = surface[1:]
        n_atoms = non_metal.count('N')
        vac = non_metal.count('C')
        if metal not in metal_list:
            ax.plot([], [], 'o', color=colors, label=metal)
            metal_list.append(metal)

        if np.abs(spin_data[metal][str(vac)][str(n_atoms)]) < 2.1 and metal not in ['Ti', 'V', 'Cr']: 
        # if metal not in ['V', 'Ti', 'Cr']:
            criteria_surface[1].append(surface)
            lines_to_plot[1].append([dG[descriptors[0]]+potential, dG[descriptors[1]]+potential ])
            marker = str(np.abs(round(spin_data[metal][str(vac)][str(n_atoms)],1)))
        ax.plot(dG[descriptors[0]]+potential, dG[descriptors[1]]+potential, color=colors, marker='$'+marker+'$')
        # elif np.abs(spin_data[metal][str(vac)][str(n_atoms)]) > 2.1 and metal not in ['Ti', 'V', 'Cr']:
        # elif np.abs(spin_data[metal][str(vac)][str(n_atoms)]) > 2.1 and metal not in ['Ti', 'V', 'Cr'] :
        #     lines_to_plot[2].append([dG[descriptors[0]]+potential, dG[descriptors[1]]+potential ])
        #     criteria_surface[2].append(surface)
        #     marker = str(np.abs(round(spin_data[metal][str(vac)][str(n_atoms)],1)))
        #     ax.plot(dG[descriptors[0]]+potential, dG[descriptors[1]]+potential, color=colors, marker='$'+marker+'^{2}$')

        # if vac == 2 and n_atoms > 0:
        dG_1.append(dG[descriptors[0]]+potential)
        dG_2.append(dG[descriptors[1]]+potential)
    
    print(f'{bcolors.OKBLUE}The following surfaces fullfill the chosen criteria: {bcolors.ENDC}')
    print(criteria_surface)
    # fit = get_fit_from_points(dG_1, dG_2, 1)
    for s in lines_to_plot:
        data = np.array(lines_to_plot[s]).transpose()
        fit = get_fit_from_points(*data, 1)
        print(f'{bcolors.OKGREEN} Fit of the scaling line {bcolors.ENDC}')
        print(fit['fit'])
        ax.plot(ax.get_xbound(), fit['p'](ax.get_xbound()), '-', color='k', lw=3)
        real_fit = fit['fit']
        # error_p_fit = np.poly1d([fit['fit'][0], fit['fit'][1] +0.4])
        # error_m_fit = np.poly1d([fit['fit'][0], fit['fit'][1] -0.4])
        # ax.plot(ax.get_xbound(), error_p_fit(ax.get_xbound()), '--', color='k', lw=3, alpha=0.5 )
        # ax.plot(ax.get_xbound(), error_m_fit(ax.get_xbound()), '--', color='k', lw=3, alpha=0.5 )

    ax.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=3)

    # Plot the parity line
    # x = ax.get_xbound()
    # ax.plot(x,x, '--', color='k', label='parity line')
    fig.tight_layout()


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
            results.setdefault(dat[0],{})[dat[2]+'_s'] = dG
    return results

@click.command()
@click.option('--pk')
def main(pk):
    node = load_node(pk)
    descriptors = node.inputs.descriptor_names.get_list()
    potential = node.inputs.voltage
#    surfaces = node.inputs.surface_names.get_list()
    surfaces = ['CoNC', 'CoNNC', 'CoNNCC', 'FeNNC', 'FeNNCC', 'FeNNNCC', 'FeCC', 'FeNNNNCC', 'NiNC', 'NiNCC', 'NiNNC', 'NiNNCC', 'NiNNNCC', 'VNCC', 'VNNC', 'VNNCC', 'VNNNCC', 'VNNNNCC', 'MnNC', 'MnNCC', 'MnNNC', 'MnNNCC', 'MnNNNCC', 'MnNNNNCC', 'RhNC', 'RhNCC', 'RhNNC', 'RhNNCC', 'RhNNNC', 'RhNNNCC', 'RhNNNNCC', 'RuNNC', 'RuNNCC', 'RuNNNC', 'RuNNNCC', 'RuNNNNCC', 'TiNC', 'TiNNC', 'TiNNCC', 'TiNNNC', 'TiNNNCC', 'TiNNNNCC', 'CrNC', 'CrNCC', 'CrNNC', 'CrNNCC', 'CrNNNC', 'CrNNNCC', 'CrNNNNCC']
    energy_file = node.inputs.energies
    points = plot_points(energy_file, surfaces)

    # coverage_map = node.outputs.coverage_map
    # # get the values needed for the title 
    # figC, axC = plot_maps(coverage_map, descriptors, 
    #                     titles=descriptors, points=points)
    # figC.savefig(os.path.join(output, 'coverage_map.pdf'))

    rate_map = node.outputs.rate_map
    reactions = node.inputs.rxn_expressions.get_list()
    figR, axR = plot_maps(rate_map, descriptors, points, potential=float(potential), type='rates')
    figR.savefig(os.path.join(output, 'rate_map_'+str(pk)+'.pdf'))

    # production_rate = node.outputs.production_rate_map.get_list()
    # pprint(production_rate)


if __name__ == '__main__':
    output = 'output'
    Path(output).mkdir( exist_ok=True, parents=True)

    main()

