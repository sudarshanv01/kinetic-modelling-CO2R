
import click
from ase.db import connect 
import matplotlib.pyplot as plt
import json
from pprint import pprint
import numpy as np
from useful_functions import create_output_directory
from plot_params import get_plot_params
from useful_functions import get_fit_from_points

def parsedb(results, database):

    for row in database.select(states='state_slab'):
        d_centre = row.d_centre
        facet = row.facets.replace('facet_','')
        metal = row.sampling.replace('sampling_','')
        results.setdefault(facet,{}).setdefault(metal,{})['d_centre'] = d_centre
        results.setdefault(facet,{}).setdefault(metal,{})['d_width'] = row.d_width
        results.setdefault(facet,{}).setdefault(metal,{})['max_hilbert'] = row.max_hilbert

def hybridization_single(eps_d, metal, adsorbate):
    ## Idealised filling factor
    f = {'Ni':0.9, 'Cu':1., 'Zn':1.1, 'Ag':1, 'Au':1, 'Pd':0.9, 'Pt':0.9}
    # Coupling matrix elements
    V_sd_sq = {'Sc':7.90, 'Ti':4.65, 'V':3.15, 'Cr':2.35, 'Mn':1.94, \
               'Fe': 1.59, 'Co':1.34, 'Ni':1.16, 'Cu':1, 'Zn':0.46,
               'Au':3.35, 'Ag':2.26, 'Pd':2.78, 'Pt':3.90, 'Ru':3.32}
    # epsilon of the state to compare against
    epsilon = {     
            'CO2':-2.00
                    }

    # # The pi states will be some some function the LMTO V_sd
    # alpha = 0.063 #CO: 0.063 #eV-1
    # beta = 1.5 #CO: 1.5 #eV^2
    # V_pi_sq = beta * V_sd_sq[metal]
    # # The overlap matrix is some factor alpha times the coupling matrix elements
    # S_pi = -alpha * np.sqrt(V_pi_sq)
    # # Assume that the sigma coupling is some function of the pi coupling term
    # # 1.3 is from Jens' paper
    # V_sigma_sq = (1.3)**2 * beta * V_sd_sq[metal]
    # # Again the overlap matrix is some function of the coupling matrix elements
    # S_sigma = -alpha * np.sqrt(V_sigma_sq)
    # # Since this NOT is a spin polarized
    # # Degen of pi = 2 and of sigma = 1
    # # Spin up and down multiply by 2
    # eps_homo, eps_lumo = epsilon[adsorbate]


    # V_species = {'C':[V_pi_sq, V_pi_sq],
    #              'CO2':[V_pi_sq, V_sigma_sq],
    #             'CH2':[V_pi_sq, V_sigma_sq],
    #             'CH3':[V_pi_sq, V_sigma_sq],
    #             'CO':[V_pi_sq, V_sigma_sq],
    #             'N':[V_pi_sq, V_pi_sq],
    #             'NH':[V_pi_sq, V_pi_sq],
    #             'NH2':[V_pi_sq, V_sigma_sq],
    #             'O':[V_pi_sq, V_pi_sq],
    #             'OH':[V_sigma_sq, V_pi_sq],
    #             'proton':[V_sigma_sq, V_sigma_sq],
    # }

    # S_species = {'C':[S_pi, S_pi],
    #             'CO2':[S_pi, S_pi],
    #             'CH2':[S_pi, S_sigma],
    #             'CH3':[S_pi, S_sigma],
    #             'CO':[S_pi, S_sigma],
    #             'N':[S_pi, S_pi],
    #             'NH':[S_pi, S_pi],
    #             'NH2':[S_pi, S_sigma],
    #             'O':[S_pi, S_pi],
    #             'OH':[S_sigma, S_pi],
    #             'proton':[S_sigma, S_sigma]
    # }

    # f_lumo, f_homo = factor_mult[adsorbate]
    # V_lumo, V_homo = V_species[adsorbate]
    # S_lumo, S_homo = S_species[adsorbate]

    # E_hyb =  -f_lumo * ( f * V_lumo / ( eps_lumo - eps_d ) \
    #               + f * S_lumo * np.sqrt(V_lumo) ) \
    #          -f_homo * ( (1-f) * V_homo / ( eps_d - eps_homo ) \
    #             + (1+f) * S_homo * np.sqrt(V_sigma_sq) )

    # mu_hyb =  -f_lumo * (  2 * f * V_lumo / ( eps_lumo - eps_d ) \
    #               + f * S_lumo * np.sqrt(V_lumo) ) \
    #          -f_homo * ( 2 * (1-f) * V_homo / ( eps_d - eps_homo ) \
    #             + (1+f) * S_homo * np.sqrt(V_sigma_sq) )

    # E = f[metal] *  V_sd_sq[metal] /  eps_d #(V_species[adsorbate][-1] - eps_d) #+ S_lumo * np.sqrt(V_lumo)

    q = V_sd_sq[metal] / ( epsilon[adsorbate] - eps_d )
    # q = V_sd_sq[metal]

    return q

    # return E_hyb, mu_hyb

@click.command()
@click.option('--qexp_dict', default='output/explicit_charge.json')
@click.option('--dbname', default='../databases/transition_metal_vacuum.db')
def main(qexp_dict, dbname):

    get_plot_params()

    results = {}
    parsedb(results, connect(dbname))
    qdict = json.load(open(qexp_dict, 'r'))
    colors = ['tab:blue', 'tab:red']

    fig, ax = plt.subplots(1, 1, figsize=(8,6))

    all_points = []
    for i, facet in enumerate(qdict):
        if '_' in facet:
            continue
        ax.plot([],[], '-o', label='('+facet+')', color=colors[i])
        for metal in qdict[facet]:
            try:
                q = results[facet][metal]['d_width'] #hybridization_single(results[facet][metal]['max_hilbert'], metal , 'CO2') #
                ax.plot(q, -1 * qdict[facet][metal], 'o', color=colors[i])
                ax.errorbar(results[facet][metal]['d_width'], -1 * qdict[facet][metal], 0.025, color=colors[i])
                ax.annotate(metal, \
                    xy=(q, -1 * qdict[facet][metal]), \
                        )
                all_points.append([q, -1 *  qdict[facet][metal]])
            except KeyError:
                continue
    all_points = np.array(all_points)
    eu_all, q_all = all_points.transpose()
    fit = get_fit_from_points(eu_all, q_all, 1)
    ax.plot(eu_all, fit['p'](eu_all), '-', lw=2 , color='k')
    ax.set_ylabel(r'$q_{CO_2^*} $ / e')
    ax.set_xlabel(r'width of metal d-states')
    ax.legend(loc='best', frameon=False)
    fig.tight_layout()
    fig.savefig('output/estructure.pdf')


if __name__ == "__main__":
    main()


