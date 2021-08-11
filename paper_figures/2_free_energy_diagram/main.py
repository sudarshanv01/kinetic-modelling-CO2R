

import numpy as np
from ase.db import connect 
import click
import argparse
from glob import glob
from computation import get_free_energy_diagram_data
from useful_functions import create_output_directory
from plot_params import get_plot_params
# from computation import plot_computational_diagram, plot_variation_with_potential
from computational_panel import FreeEnergyDiagram, plot_variation_with_potential, plot_computational_diagram
import matplotlib.pyplot as plt
from experimental import plot_experimental_data
import csv
from ase.data import atomic_numbers
from ase.data.colors import jmol_colors
import json
import string
from pprint import pprint
from molecule import plot_molecule
import matplotlib as mpl

"""
Create Figure 1 of the paper
"""

def cli_parse():
    ## All the inputs
    parser = argparse.ArgumentParser()
    parser.add_argument('--database_folder', default='../databases/', help='Database folder')
    parser.add_argument('--potential', default=[-0.6, -0.8, -1.], nargs='*', type=float, \
                            help='SHE potential')
    parser.add_argument('--ph', default=2., type=float, help='pH')
    parser.add_argument('--referencedb_name', default='/Users/vijays/Documents/common_databases/gas_phase.db')
    parser.add_argument('--ninc_experiment', default='inputs/pH_effect_NiNC.xlsx')
    parser.add_argument('--fenc_experiment', default='inputs/pH_effect_FeNC.xlsx')
    parser.add_argument('--fenc_xilehu_experiment', default='inputs/pH_effect_FeNC_XileHu.xlsx')
    parser.add_argument('--gold_experiment', default='inputs/pH_effect_Gold.xlsx')
    parser.add_argument('--copc_experiment', default='inputs/pH_effect_CoPc.xlsx')
    parser.add_argument('--molecular_database', default='inputs/molecule_CO2R.db' )
    return parser.parse_args()


def main():
    
    ## Figure specifications
    fig = plt.figure(figsize=(20,12.5))
    gs = fig.add_gridspec(8, 4)
    ax_Au = fig.add_subplot(gs[4:,0])
    ax_Fe = fig.add_subplot(gs[4:,1])
    ax_Ni = fig.add_subplot(gs[4:,2])
    ax_Co = fig.add_subplot(gs[4:,3:])
    cax_Au = fig.add_subplot(gs[2:4,0])
    cax_Fe = fig.add_subplot(gs[2:4,1])
    cax_Ni = fig.add_subplot(gs[2:4,2])
    cax_Co = fig.add_subplot(gs[2:4,3:])
    axf_Au = fig.add_subplot(gs[0:2,0])
    axf_MNC = fig.add_subplot(gs[0:2,1:3])
    axf_CoPc = fig.add_subplot(gs[0:2,3])

    ax = [cax_Au, cax_Fe, cax_Ni, cax_Co, ax_Au, ax_Fe, ax_Ni, ax_Co]

    ## parser information
    parser = cli_parse()
    databases = glob(parser.database_folder + '/*.db')
    data = {}
    plot_experimental_data(parser.gold_experiment,ax_Au, 
                r'pH independent', r'Au(pc)', jmol_colors[atomic_numbers['Au']],
                fit_min=-0.75, fit_lim=-1.05)
    plot_experimental_data(parser.fenc_experiment,ax_Fe, 
                r'pH independent', r'FeNC' ,jmol_colors[atomic_numbers['Fe']], 
                fit_min=-0.5, fit_lim=-0.85)
    plot_experimental_data(parser.fenc_xilehu_experiment,ax_Fe, 
                r'pH independent', r'FeNC' ,jmol_colors[atomic_numbers['Fe']], 
                fit_min=-0.5, fit_lim=-0.85, marker='v', plot_line=False)
    ax_Fe.annotate(r'80 $\frac{\mathregular{mV}}{\mathregular{dec}}$', xy=(0.75,0.8), xycoords='axes fraction')
    plot_experimental_data(parser.ninc_experiment,ax_Ni, 
                r'pH dependent', r'NiNC', jmol_colors[atomic_numbers['Ni']], 
                fit_lim=-1.05, fit_min=-0.75, fit_all=False)
    plot_experimental_data(parser.copc_experiment, ax_Co, 
                r'pH dependent', r'CoPc', jmol_colors[atomic_numbers['Co']],
                fit_lim=-1., fit_min=-0.5, fit_all=False, pH_material='Co')

    ## Get the Free energy diagram
    for potential in parser.potential:
        method = FreeEnergyDiagram(dbnames=databases,\
                                    refdbname=parser.referencedb_name,\
                                    potential=potential, 
                                    pH=parser.ph)
        method.main()
        # data[potential], writeout, writeout_zero, explicit_charge, E0  = diagram_data

        data[potential] = method.diagram
        writeout = method.writeout
        writeout_zero = method.writeout_zero
        explicit_charge = method.explicit_charge
        E0 = method.E0

        filename = 'output/catmap_potential_%1.2f.txt'%potential
        with open(filename, 'w') as handle:
            csvwriter = csv.writer(handle, delimiter='\t')
            for row in writeout:
                csvwriter.writerow(row)
        
    ## do the plot ata for the molecular part
    ## this conforms with the new way of doing finite difference
    ## using the newer implementation, which is why it is a new
    ## class
    plot_molecule([parser.potential[1]], parser.ph, parser.molecular_database, cax_Co, method.references, method.references_E)

    ## save the catmap input file
    filename = 'output/catmap_potential_pzc.txt'
    with open(filename, 'w') as handle:
        csvwriter = csv.writer(handle, delimiter='\t')
        for row in writeout_zero:
            csvwriter.writerow(row)
    pprint(explicit_charge)
    
    with open('../databases/explicit_charge.json', 'w') as handle:
        json.dump(explicit_charge, handle)

    with open('../databases/zero_charge_energies.json', 'w') as handle:
        json.dump(E0, handle)

    plot_computational_diagram(data, [cax_Au, cax_Fe, cax_Ni], SAC_potential=-0.8)

    ## Label the diagram
    alphabet = list(string.ascii_lowercase)
    for i, a in enumerate(ax):
        a.annotate(alphabet[i]+')', xy=(0.05, 0.87), xycoords='axes fraction', fontsize=20)

    ## add in images 

    arr_image = plt.imread('input_images/Au27.020.png', format='png')
    axf_Au.imshow(arr_image)
    axf_Au.axis('off')

    arr_image = plt.imread('input_images/defects.png', format='png')
    axf_MNC.imshow(arr_image)
    axf_MNC.axis('off')

    arr_image = plt.imread('input_images/CoPc.png', format='png')
    axf_CoPc.imshow(arr_image)
    axf_CoPc.axis('off')


    fig.tight_layout()
    fig.savefig('output/figure1.pdf')

    # fig3, ax3 = plot_variation_with_potential(explicit_charge)
    # fig3.savefig('output/figure3.pdf')




if __name__ == '__main__':
    """
    Main set of classes and functions that 
    plot Figure 1 of the manuscript. 
    The plotting script includes: 
    1. Computational free energy diagram for Au, FeNC and NiNC
    2. Replotting experimental data from Wen and Au
    """
    create_output_directory()
    create_output_directory('output_si')
    get_plot_params()
    main()
