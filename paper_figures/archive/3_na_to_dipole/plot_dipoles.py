

import numpy as np
# from scipy.signal import hilbert
from scipy import signal
import matplotlib.pyplot as plt 
from useful_functions import create_output_directory
from pprint import pprint
import numpy as np
from plot_params import get_plot_params
import click
import json
from scipy.signal import hilbert
from ase.db import connect
from useful_functions import get_fit_from_points
from ase.data import atomic_numbers
from ase.data.colors import jmol_colors
from scipy.integrate import quad, simps
from scipy.optimize import curve_fit

def parsedb(results, database, sac=False, dbconfig={}):

    for row in database.select(*dbconfig):
        if sac: 
            facet = row.vacancy_number.replace('vacancy_','') + '_' + row.dopant_number.replace('dopant_','')
            metal = row.metal_dopant.replace('metal_dopant_','').replace('_nonorth','')
        else:
            facet = row.facets.replace('facet_','')
            metal = row.sampling.replace('sampling_','')
        state = row.states.replace('state_','')
        results.setdefault(facet,{}).setdefault(metal,{}).setdefault(state,{})['dipole'] = row.dipole_field


def main():

    tmdbname = '../databases/transition_metal_vacuum.db'
    sacdbname = '../databases/single_atom_vacuum.db' 
    
    results = {}
    parsedb(results, connect(tmdbname),)
    parsedb(results, connect(sacdbname), sac=True)

    fig, ax = plt.subplots(1, 1, figsize=(9,4.5))

    for facet in results:
        if facet == '111':
            continue
        if '_' in facet:
            facetname = '(' + facet.replace('_',',') + ')'
            color = 'tab:red'
        else:
            facetname = '(' + facet + ')'
            color='tab:blue'
        for metal in results[facet]:
            if metal == 'gas_molecules':
                continue

            xlabel = metal + facetname
            print(metal, facet, results[facet][metal])
            ax.plot(xlabel, results[facet][metal]['CO']['dipole'], 'o', color=color)
            ax.plot(xlabel, results[facet][metal]['COOH']['dipole'], 'v', color=color)
            ax.plot(xlabel, results[facet][metal]['CO2_dos']['dipole'], '*', color=color)
    plt.xticks(rotation=90)
    ax.set_ylabel(r'$\mu$ / e$\AA$')
    ax.plot([],[],'*', color='k', label=r'CO$_2$*')
    ax.plot([],[],'v', color='k', label='COOH*')
    ax.plot([],[],'o', color='k', label='CO*')
    # ax.legend(loc='best', frameon=False)
    # ax.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                # mode="expand", borderaxespad=0, ncol=3, frameon=False, fontsize=14)
    # ax.legend( bbox_to_anchor=(1.04,1), borderaxespad=0, fontsize=12)
    ax.annotate('Transition Metals', xy=(0.1, 0.9), color='tab:blue', xycoords='axes fraction', fontsize=12)
    ax.annotate('MNC', xy=(0.7,0.9), color='tab:red', xycoords='axes fraction', fontsize=12)
    ax.legend(bbox_to_anchor=(1.04,0), loc="lower left", borderaxespad=0, fontsize=14, frameon=False)
    fig.tight_layout()
    fig.savefig('output/dipole_compare.pdf')


    tm_files = '../databases/TM_dos.json'
    sac_files = '../databases/SAC_dos.json'
    gas_files = '../databases/Gas_dos.json'

    tm_data = json.load(open(tm_files, 'r'))
    sac_data = json.load(open(sac_files, 'r'))
    gas_data = json.load(open(gas_files, 'r'))    

    fig, ax = plt.subplots(1, 7, figsize=(14,4))
    i = 0
    for metal in tm_data:
        if metal == 'Al':
            continue
        for j, facet in enumerate(tm_data[metal]):
            if facet == '100':
                continue

            pdos_slab = tm_data[metal][facet]['pdos']['%s(d)'%metal]['+'][0]
            pdos_ads = tm_data[metal][facet]['pdos']['CO2(sp)']['+'][0]
            energies_ads = tm_data[metal][facet]['energies_ads']
            energies_slab = tm_data[metal][facet]['energies_slab']

            ax[i].fill_between(pdos_slab,\
                    energies_slab,\
                    color=jmol_colors[atomic_numbers[metal]])

            ax[i].plot(pdos_ads,\
                    energies_ads,\
                    color='tab:green')
            ax[i].set_ylim([-5,5])
            i += 1

    for metal in sac_data:
        for value in sac_data[metal]:
            vacancy, dopant = value.split('_')
            if vacancy  == '2' and dopant == '4':
                pass
            else:
                continue
            print(metal, vacancy, dopant)
            pdos_slab = []
            pdos_ads = []
            for spin in sac_data[metal][value]['pdos']['%s(d)'%metal]:
                pdos_slab.append(sac_data[metal][value]['pdos']['%s(d)'%metal][spin][0])
                pdos_ads.append(sac_data[metal][value]['pdos']['CO2(sp)'][spin][0])

            energies_ads = sac_data[metal][value]['energies_ads']
            energies_slab = sac_data[metal][value]['energies_slab']

            summed_dos_slab = np.array(pdos_slab).sum(axis=0)
            summed_dos_ads = np.array(pdos_ads).sum(axis=0)
            ax[i].plot( summed_dos_ads, energies_ads, color='tab:green')
            ax[i].fill_between(summed_dos_slab, energies_slab, color=jmol_colors[atomic_numbers[metal]])
            ax[i].set_ylim([-5,5])
            i+= 1

    fig.tight_layout()
    fig.savefig('output/pdos_compare.pdf')


if __name__ == "__main__":
    get_plot_params()
    main()