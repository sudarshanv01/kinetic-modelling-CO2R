


import numpy as np
from ase.db import connect
import matplotlib.pyplot as plt
import click
from pprint import pprint
from ase.data.colors import jmol_colors
from ase.data import atomic_numbers
from useful_functions import create_output_directory
from plot_params import get_plot_params
from matplotlib.patches import Circle
from ase.data import covalent_radii as radii


def parsedb(database, results):
    for row in database.select():
        results.setdefault(row.proton_conc,{})['pdos'] = row.data.pdos.pdos
        results.setdefault(row.proton_conc,{})['energy'] = row.data.pdos.energies
        results.setdefault(row.proton_conc,{})['integrated_dos'] = row.data.pdos.integrated_dos
        results.setdefault(row.proton_conc,{})['total_dos'] = row.data.pdos.total_dos
        results.setdefault(row.proton_conc,{})['wf'] = row.wf
        results.setdefault(row.proton_conc,{})['atoms'] = row.toatoms()
        results.setdefault(row.proton_conc,{})['magmom'] = row.magmom

def parsedbimp(database, results):
    for row in database.select():
        results.setdefault(row.tot_charge,{})['pdos'] = row.data.pdos.pdos
        results.setdefault(row.tot_charge,{})['energy'] = row.data.pdos.energies
        results.setdefault(row.tot_charge,{})['integrated_dos'] = row.data.pdos.integrated_dos
        results.setdefault(row.tot_charge,{})['total_dos'] = row.data.pdos.total_dos
        results.setdefault(row.tot_charge,{})['wf'] = row.wf
        results.setdefault(row.tot_charge,{})['atoms'] = row.toatoms()
        results.setdefault(row.tot_charge,{})['magmom'] = row.magmom

@click.command()
@click.option('--dbname', default='databases/dos_evolve.db')
@click.option('--dbnameimp', default='databases/dos_evolve_implicit.db')
def main(dbname, dbnameimp):
    results = {}
    parsedb(connect(dbname), results)    
    imp_results = {}
    parsedbimp(connect(dbnameimp), imp_results)

    fig, ax = plt.subplots(2, len(results), figsize=(12,7), squeeze=False)

    for i, conc in enumerate(sorted(results)):
        wf = results[conc]['wf']
        pdos_Fe = np.array(results[conc]['pdos']['Fe']['+']) \
                 +np.array(results[conc]['pdos']['Fe']['-'])
        pdos_C = np.array(results[conc]['pdos']['C']['+']) \
                 +np.array(results[conc]['pdos']['C']['-'])
        pdos_N = np.array(results[conc]['pdos']['N']['+']) \
                 +np.array(results[conc]['pdos']['N']['-']) 
        indices_occupied = [i for i in range(len(results[conc]['energy'])) if results[conc]['energy'][i] < 0.0 ]
        total_dos = np.array(results[conc]['total_dos'][0]) +  np.array(results[conc]['total_dos'][1]) 
        integrated_dos = np.array(results[conc]['integrated_dos'][0]) +  np.array(results[conc]['integrated_dos'][1]) 
        ax[0,i].fill_between(np.array(pdos_Fe[0])[indices_occupied], np.array(results[conc]['energy'])[indices_occupied], color=jmol_colors[atomic_numbers['Fe']], alpha=0.25)
        ax[0,i].fill_between(np.array(pdos_C[0])[indices_occupied], np.array(results[conc]['energy'])[indices_occupied], color=jmol_colors[atomic_numbers['C']], alpha=0.25)
        ax[0,i].fill_between(np.array(pdos_N[0])[indices_occupied], np.array(results[conc]['energy'])[indices_occupied], color=jmol_colors[atomic_numbers['N']], alpha=0.25)
        ax[0,i].plot(pdos_N[0], results[conc]['energy'], color=jmol_colors[atomic_numbers['N']], label='N')
        ax[0,i].plot(pdos_C[0], results[conc]['energy'], color=jmol_colors[atomic_numbers['C']], label='C')
        ax[0,i].plot(pdos_Fe[0], results[conc]['energy'], color=jmol_colors[atomic_numbers['Fe']], label='Fe')
        
        ## get needed integrals for the 
        # indices_occupied = [i for i in range(len(results[conc]['energy'])) if results[conc]['energy'][i] < 0.0 ]
        # spin_up_occupied = np.trapz(np.array(results[conc]['pdos']['Fe']['+'][0])[indices_occupied], np.array(results[conc]['energy'])[indices_occupied])
        # spin_down_occupied = np.trapz(np.array(results[conc]['pdos']['Fe']['-'][0])[indices_occupied], np.array(results[conc]['energy'])[indices_occupied])
        # spin_up = np.trapz(np.array(results[conc]['pdos']['Fe']['+'][0]), np.array(results[conc]['energy']))
        # spin_down = np.trapz(np.array(results[conc]['pdos']['Fe']['-'][0]), np.array(results[conc]['energy']))
        # print(spin_up - spin_down)
        # print(conc, spin_up_occupied / spin_up - spin_down_occupied / spin_down)

        if i == 0:
            ax[0,i].legend(loc='lower right', fontsize=12, frameon=False)
        
        atoms = results[conc]['atoms']
        for atom in atoms:
            color = jmol_colors[atom.number]
            # metal = atom.symbol
            radius = radii[atom.number]
            circle = Circle((atom.x, atom.z), radius, facecolor=color,
                                    edgecolor='k', linewidth=1)
            ax[1,i].add_patch(circle) 
        ax[1,i].axis('equal')
        ax[1,i].set_xticks([])
        ax[1,i].set_yticks([])
        ax[1,i].axis('off')
        ax[1,i].set_title('Mag. mom. = %1.1f'%results[conc]['magmom'], y=-0.2 )
        ax[1,i].annotate(r'$n_{\mathregular{H}^{+}} = $ %s'%conc.replace('proton_conc_',''), \
                        xy=(0.3,0.8), xycoords='axes fraction')

        ax[0,i].set_xticks([])
        ax[0,i].set_ylabel(r'$\epsilon - \epsilon_{f}$ / eV')
        ax[0,i].set_xlabel(r'PDOS / arb units')
        ax[0,i].set_ylim([-10, 5])
        ax[0,i].set_xlim([-0.01, 0.4])
        ax[0,i].set_title(r'$\phi = %1.1f $ eV'%wf)
        ax[0,i].axhline(y=0, color='k', lw=2)

    fig.tight_layout()
    fig.savefig('output/pdos.pdf')

    # plotting implicit results
    fig, ax = plt.subplots(1, len(imp_results)-3, figsize=(20,4), squeeze=False)

    for i, conc in enumerate(sorted(imp_results)):
        wf = imp_results[conc]['wf']
        if wf < 1.9:
            continue
        pdos_Fe = np.array(imp_results[conc]['pdos']['Fe']['+']) \
                 +np.array(imp_results[conc]['pdos']['Fe']['-'])
        pdos_C = np.array(imp_results[conc]['pdos']['C']['+']) \
                 +np.array(imp_results[conc]['pdos']['C']['-'])
        pdos_N = np.array(imp_results[conc]['pdos']['N']['+']) \
                 +np.array(imp_results[conc]['pdos']['N']['-']) 
        indices_occupied = [i for i in range(len(imp_results[conc]['energy'])) if imp_results[conc]['energy'][i] < 0.0 ]
        total_dos = np.array(imp_results[conc]['total_dos'][0]) +  np.array(imp_results[conc]['total_dos'][1]) 
        integrated_dos = np.array(imp_results[conc]['integrated_dos'][0]) +  np.array(imp_results[conc]['integrated_dos'][1]) 
        ax[0,i].fill_between(np.array(pdos_Fe[0])[indices_occupied], np.array(imp_results[conc]['energy'])[indices_occupied], color=jmol_colors[atomic_numbers['Fe']], alpha=0.25)
        ax[0,i].fill_between(np.array(pdos_C[0])[indices_occupied], np.array(imp_results[conc]['energy'])[indices_occupied], color=jmol_colors[atomic_numbers['C']], alpha=0.25)
        ax[0,i].fill_between(np.array(pdos_N[0])[indices_occupied], np.array(imp_results[conc]['energy'])[indices_occupied], color=jmol_colors[atomic_numbers['N']], alpha=0.25)
        ax[0,i].plot(pdos_N[0], imp_results[conc]['energy'], color=jmol_colors[atomic_numbers['N']], label='N')
        ax[0,i].plot(pdos_C[0], imp_results[conc]['energy'], color=jmol_colors[atomic_numbers['C']], label='C')
        ax[0,i].plot(pdos_Fe[0], imp_results[conc]['energy'], color=jmol_colors[atomic_numbers['Fe']], label='Fe')
        ax[0,i].set_xticks([])
        ax[0,i].set_ylabel(r'$\epsilon - \epsilon_{f}$ / eV')
        ax[0,i].set_xlabel(r'PDOS / arb units')
        ax[0,i].set_ylim([-10, 5])
        ax[0,i].set_xlim([-0.01, 0.4])
        ax[0,i].set_title(r'$\phi = %1.1f $ eV'%wf)
        ax[0,i].axhline(y=0, color='k', lw=2)
        ax[0,i].annotate('Mag. mom. = %1.1f'%imp_results[conc]['magmom'], xy=(0.4,0.05),xycoords='axes fraction', fontsize=12 )

    fig.tight_layout()
    fig.savefig('output/imp_pdos.pdf')

if __name__ == "__main__":
    create_output_directory()
    get_plot_params()
    main()