"""Create Figure 1 of the paper."""
import json
import string
import argparse
from pprint import pprint
from pathlib import Path
import numpy as np
import argparse
from glob import glob
from plot_params import get_plot_params
import matplotlib.pyplot as plt
import csv
from ase.data import atomic_numbers
from ase.data.colors import jmol_colors
from molecule import plot_molecule
from experimental import plot_experimental_data
from computational_panel import FreeEnergyDiagram, plot_computational_diagram

# Output folder
Path('output').mkdir(parents=True, exist_ok=True)
Path('output_si').mkdir(parents=True, exist_ok=True)


def cli_parse():
    """Parse all the inputs, defaults are the ones used in the paper."""
    parser = argparse.ArgumentParser()
    parser.add_argument('--database_folder', default='../databases/', help='Database folder')
    parser.add_argument('--potential', default=[  -0.4, -0.6, -0.8, ], nargs='*', type=float, \
                            help='SHE potential')
    parser.add_argument('--chosen_potential', default=-0.8, help='Potential to plot all single\
                                            atom and CoPc structures', type=float)
    parser.add_argument('--ph', default=2., type=float, help='pH')
    parser.add_argument('--referencedb_name', default='input_databases/gas_phase.db')
    parser.add_argument('--ninc_experiment', default='inputs/pH_effect_NiNC.xls')
    parser.add_argument('--fenc_experiment', default='inputs/pH_effect_FeNC.xls')
    parser.add_argument('--fenc_xilehu_experiment', default='inputs/pH_effect_FeNC_XileHu.xls')
    parser.add_argument('--gold_experiment', default='inputs/pH_effect_Gold.xls')
    parser.add_argument('--copc_experiment', default='inputs/pH_effect_CoPc.xls')
    parser.add_argument('--molecular_database', default='input_databases/molecule_CO2R.db' )
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

    # parser information
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

    # Get the Free energy diagram
    for potential in parser.potential:
        method = FreeEnergyDiagram(dbnames=databases,\
                                    refdbname=parser.referencedb_name,\
                                    potential=potential, 
                                    pH=parser.ph)
        method.main()

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
        
    # do the plot data for the molecular part
    # this conforms with the new way of doing finite difference
    # using the newer implementation, which is why it is a new class
    plot_molecule([parser.chosen_potential], parser.ph, parser.molecular_database, cax_Co, method.references, method.references_E)

    # save the catmap input file
    filename = 'output/catmap_potential_pzc.txt'
    with open(filename, 'w') as handle:
        csvwriter = csv.writer(handle, delimiter='\t')
        for row in writeout_zero:
            csvwriter.writerow(row)
    print('-------')
    print('Explicit charge is:')
    pprint(explicit_charge)
    
    with open('../databases/explicit_charge.json', 'w') as handle:
        json.dump(explicit_charge, handle)

    with open('../databases/zero_charge_energies.json', 'w') as handle:
        json.dump(E0, handle)

    plot_computational_diagram(data, [cax_Au, cax_Fe, cax_Ni], SAC_potential=parser.chosen_potential)

    # Label the diagram
    alphabet = list(string.ascii_lowercase)
    for i, a in enumerate(ax):
        a.annotate(alphabet[i]+')', xy=(0.05, 0.87), xycoords='axes fraction', fontsize=20)

    # add in images 
    arr_image = plt.imread('input_images/Au27.png', format='png')
    axf_Au.imshow(arr_image)
    axf_Au.axis('off')

    arr_image = plt.imread('input_images/defects.png', format='png')
    axf_MNC.imshow(arr_image)
    axf_MNC.axis('off')

    arr_image = plt.imread('input_images/CoPc.png', format='png')
    axf_CoPc.imshow(arr_image)
    axf_CoPc.axis('off')

    fig.tight_layout()
    fig.savefig('output/figure1.png')



if __name__ == '__main__':
    """
    Main set of classes and functions that 
    plot Figure 1 of the manuscript. 
    The plotting script includes: 
    1. Computational free energy diagram for Au, FeNC and NiNC
    2. Replotting experimental data from Wen and Au
    """
    get_plot_params()
    main()
