

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

def get_filled_states(energies, pdos):
    # does the integral of -inf to Ef for the pdos
    from scipy.integrate import quad, simps
    energy_fill = []
    pdos_occ = []
    for index, energy in enumerate(energies):
        if energy < 0:
            energy_fill.append(energy)
            pdos_occ.append(pdos[index])
    def integrand(e, rho):
        return rho

    f = simps(integrand(energy_fill, pdos_occ), energy_fill)

    return f

def moment_generator(energies, all_dos, moment):

    # Determine the reference value
    # In this case the reference is the d-band center

    def integrand_numerator(e, rho):
        return e * rho
    def integrand_denominator(e,rho):
        return rho
    def integrand_moment_numerator(e, epsilon_d, moment, rho):
        return ( (e - epsilon_d ) ** moment ) * rho

    moments = []
    eps_d = []

    for list_dos in all_dos:
        dos = np.array(list_dos)
        epsilon_d_num = simps(integrand_numerator(energies,dos), energies)
        epsilon_d_den = simps(integrand_denominator(energies,dos), energies)
        epsilon_d = epsilon_d_num / epsilon_d_den

        moment_numerator = simps( integrand_moment_numerator(energies, epsilon_d, moment, \
                dos), energies)
        moment_denom = epsilon_d_den

        moment_spin = moment_numerator / moment_denom
        moments.append(moment_spin)
        eps_d.append(epsilon_d)
    return moments, eps_d

def adsorbate_states(delta, hilbert, eps_a, eps, V_sq=1):
    na =  1 / np.pi * ( delta ) / ( ( (eps - eps_a)/V_sq - hilbert )**2 + delta**2  ) 
    return na

def semi_elliptical_dos(a, b, c, energy):
    """Plotting semi-ellipse as Delta (state of the metal)
    Equation is of the form
    delta = a ( c - energy**2 / b ) ** 0.5
    """
    delta = []
    for E in energy:
        d = a * ( c - E**2 / b) #** 1/2
        if d < 0:
            delta.append(0)
        else:
            delta.append(d)
    area = np.trapz(energy, delta)
    # print(area)
    return delta

def laurentz_dos(a, b1,energy):
    delta = 1 / np.pi * ( a / ( (energy - b1)**2 + a**2  ) ) #+ w * 1 / np.pi * ( a / ( (energy - b2)**2 + a**2  ) )  
    return delta

def hybridization_energy(energy, e_a, Delta, Lambda, V_sq):
    function = (2 / np.pi) * ( np.arctan( Delta / ( (energy - e_a)/V_sq - Lambda ) ) ) 
    dE_hyb = np.trapz(function, energy) - e_a
    return dE_hyb

def schematic(ideal_params, eps_a=0):

    fig, ax = plt.subplots(len(ideal_params), len(ideal_params['100']), figsize=(16,6))
    energy_val = np.linspace(-6,6,500)

    V_metal = {'Cu':1, 'Ag':2.26, 'Au':3.35, 'Pt':3.90, 'Pd':2.78, 'Al':1}

    # j=0
    for j, facet in enumerate(ideal_params):
        for i, metal in enumerate(ideal_params[facet]):
            V = V_metal[metal]

            # delta = laurentz_dos(*values[:-2]+[energy_val])
            parameters = ideal_params[facet][metal]
            delta = np.array(fit_semi_ellipse(energy_val, *parameters))
            hilbert = np.imag(signal.hilbert(delta))
            hilbert_prime = np.diff(hilbert) / np.diff(energy_val)
            n_occupancy = 1 / ( 1 - hilbert_prime )

            # ## find the energy value at which epsilon = epsilon_a + lambda
            root_match_1 = energy_val - eps_a
            root_match_2 = hilbert
            na = adsorbate_states(np.array(delta), hilbert, eps_a, energy_val, V_sq=V )
            
            ax[j,i].plot(delta, energy_val, color='tab:red', lw=4)
            ax[j,i].plot(na, energy_val, color='k')
            ads_range = np.linspace(*ax[j,i].get_ybound())
            ax[j,i].plot((ads_range-eps_a)/V, ads_range, color='tab:grey', ls='--', alpha=0.25)
            ax[j,i].plot(hilbert, energy_val, color='tab:blue', alpha=0.25)

            ax[j,i].set_xlim([-1.5, 3])
            ax[j,i].axhline(eps_a, color='tab:green', ls='-', alpha=0.25)
            filled_indices = [a for a in range(len(energy_val)) if energy_val[a] < 0 ]
            filled_na = na[filled_indices]
            occupancy = np.trapz(filled_na, energy_val[filled_indices]) / np.trapz(na, energy_val)
            ax[j,i].fill_between(filled_na, energy_val[filled_indices], color='k')


            ## get the hybridization energy 
            E_hyb = hybridization_energy(energy_val[filled_indices], eps_a, delta[filled_indices], hilbert[filled_indices], V)

            ax[j,i].set_xticks([])
            ax[j,i].annotate(r'$E_{hyb} = %1.2f eV$'%E_hyb, xy=(0.4, 0.8), xycoords='axes fraction', fontsize=14)
            # ax[i].set_title(r'$n_a = %1.1f eV$'%(values[1]/2))
            ax[j,i].annotate(r'$\left < n_a \right >  = %1.2f e$'%occupancy, xy=(0.4, 0.7), xycoords='axes fraction', fontsize=14)
            # ax[i].annotate(r'$\Lambda$', xy=(0.6, 0.8), xycoords='axes fraction', color='tab:blue', alpha=0.5)
            # ax[i].annotate(r'$n_a$', color='k', xy=(0.7,0.35), xycoords='axes fraction')
            # ax[i].annotate(r'$\epsilon-\epsilon_a$', color='tab:grey',  xy=(0.05,0.3), xycoords='axes fraction')

            # if i == 0:
            #     ax[i].annotate(r'$\Delta$: metal', color='tab:red', xy=(0.05,0.8), xycoords='axes fraction', fontsize=14)
            #     ax[i].annotate(r'$\Lambda$: Transformation', color='tab:blue', xy=(0.05,0.7), xycoords='axes fraction', fontsize=14)
            # elif i == 1:
            #     ax[i].annotate(r'$n_a (\epsilon) = \frac{1}{\pi} \frac{\Delta}{\left ( \epsilon - \epsilon_a - \Lambda \right )^2 + \Delta^2}$', \
            #         xy=(0.02, 0.8), xycoords='axes fraction', fontsize=14)
        j+= 1
    for i in range(len(ax)):
        ax[i,0].set_ylabel(r'$\epsilon - \epsilon_{F}$')
    fig.tight_layout()
    fig.savefig('output/dos.pdf')

def parsedb(results, database, sac=False, dbconfig={}):

    for row in database.select(*dbconfig):
        # try: 
        #     d_centre = row.d_centre
        # except AttributeError:
        #     continue
        if sac: 
            facet = row.vacancy_number.replace('vacancy_','') + '_' + row.dopant_number.replace('dopant_','')
            metal = row.metal_dopant.replace('metal_dopant_','').replace('_nonorth','')
        else:
            facet = row.facets.replace('facet_','')
            metal = row.sampling.replace('sampling_','')
        # results.setdefault(facet,{}).setdefault(metal,{})['d_centre'] = d_centre
        # results.setdefault(facet,{}).setdefault(metal,{})['d_width'] = row.d_width
        # results.setdefault(facet,{}).setdefault(metal,{})['max_hilbert'] = row.max_hilbert
        results.setdefault(facet,{}).setdefault(metal,{})['dipole'] = row.dipole_field
        # try:
        #     results.setdefault(facet,{}).setdefault(metal,{})['magmom'] = row.magmom
        # except AttributeError:
        #     results.setdefault(facet,{}).setdefault(metal,{})['magmom'] = 0

def fit_semi_ellipse(energy, a, b, epsilon_d):
    # ipdos = np.max(norm_pdos_metal) / np.sqrt(second_moment_metal) * ( second_moment_metal*2 - (E-epsilon_d_metal)**2 )**0.5
    dos = a / np.sqrt(b) * ( b - (energy - epsilon_d)**2  )**0.5
    # dos_semi_ellipse = [max(x, 0) for x in dos]
    dos_semi_ellipse = []
    for x in dos:
        if np.isfinite(x):
            dos_semi_ellipse.append(x)
        else:
            dos_semi_ellipse.append(0)
    return dos_semi_ellipse

def plot_dos():

    tm_files = '../databases/TM_dos.json'
    sac_files = '../databases/SAC_dos.json'
    gas_files = '../databases/Gas_dos.json'

    tm_data = json.load(open(tm_files, 'r'))
    sac_data = json.load(open(sac_files, 'r'))
    gas_data = json.load(open(gas_files, 'r'))    

    fig, ax = plt.subplots(2, 6, figsize=(20,7))
    figA, axA = plt.subplots(2, 6, figsize=(20,7), sharex=True, sharey=True)
    figS, axS = plt.subplots(1, 1, figsize=(3,4))
    dipole_desc = {}
    
    ideal_params = {}

    for i, metal in enumerate(tm_data):
        for j, facet in enumerate(tm_data[metal]):

            pdos_slab = tm_data[metal][facet]['pdos']['%s(d)'%metal]['+'][0]
            pdos_ads = tm_data[metal][facet]['pdos']['CO2(sp)']['+'][0]
            energies_ads = tm_data[metal][facet]['energies_ads']
            energies_slab = tm_data[metal][facet]['energies_slab']

            ax[j,i].fill_between(pdos_slab,\
                    energies_slab,\
                    color=jmol_colors[atomic_numbers[metal]])

            axA[j,i].fill_between(pdos_slab,\
                    energies_slab,\
                    color=jmol_colors[atomic_numbers[metal]])

            axA[j,i].plot(pdos_ads,\
                    energies_ads,\
                    color='tab:green')
            axA[j,i].set_xticks([])

            second_moment, eps_d = moment_generator(energies_slab, [pdos_slab], 2)
            width = 4 * np.sqrt(second_moment)
            ## get occupancy of adsorbate
            # zero_moment_ads, eps_ads = moment_generator(energies_ads, [pdos_ads], 0)
            zero_moment_ads = get_filled_states(energies_ads, pdos_ads)

            popt, pcov = curve_fit(lambda energy, a, b: fit_semi_ellipse(energy, a, b, eps_d), \
                        energies_slab, pdos_slab, [np.max(pdos_slab), second_moment[0]])
            
            ideal_params.setdefault(facet,{})[metal] = popt.tolist() + eps_d

            ideal_pdos = fit_semi_ellipse(np.array(energies_slab), popt[0], popt[1], eps_d)#np.array(ideal_pdos)
            ax[j,i].plot(ideal_pdos, energies_slab, '-', color='k')
            axA[j,i].annotate(metal, xy=(0.8,0.8), color=jmol_colors[atomic_numbers[metal]])


            hilbert_t = np.imag(signal.hilbert(pdos_slab)) 
            hilbert_ideal = np.imag(signal.hilbert(ideal_pdos))
            hilbert_prime = np.diff(hilbert_t) / np.diff(energies_slab)
            ax[j,i].plot(hilbert_t, energies_slab, color='tab:green')
            ax[j,i].plot(hilbert_ideal, energies_slab, color='tab:red')
            # axA[j,i].plot(hilbert_ideal, energies_slab, color='tab:red')
                
            # ax[j,i].plot(hilbert_prime, \
            #             energies_slab[:-1],
            #         color='tab:red', alpha=0.5)

            ax[j,i].set_ylim([-10,5])
            axA[j,i].set_ylim([-10,5])
            f = get_filled_states(energies_slab, pdos_slab)
            dipole_desc.setdefault(facet,{})[metal] = eps_d[0]#energies_slab[np.argmax(hilbert_t)] #eps_d[0]##width[0] 

    for axis in ax:
        for a in axis:
            a.set_xlabel(r'Energies / eV')
            a.axhline(0, color='k', ls='-')
    ax[0,0].set_ylabel(r'PDOS')
    ax[1,0].set_ylabel(r'PDOS')
    
    fig.tight_layout()
    fig.savefig('output/pdos_slab.pdf')
    figA.tight_layout()
    figA.savefig('output/pdos_ads.pdf')
    figS.savefig('output/single_metal_pdos_schematic.pdf')



    for metal in sac_data:
        fig, ax = plt.subplots(2, 4, figsize=(14,7),sharex=True, sharey=True)
        figA, axA = plt.subplots(2, 4, figsize=(14,7), sharex=True, sharey=True)
        # axA[0,0].plot(gas_data['gas_molecules']['none']['pdos']['CO2_bent(sp)']['+'][0], \
        #             gas_data['gas_molecules']['none']['energies_ads'], color='tab:blue')
        # axA[0,0].plot(-1*np.array(gas_data['gas_molecules']['none']['pdos']['CO2_bent(sp)']['-'][0]), \
        #             gas_data['gas_molecules']['none']['energies_ads'], color='tab:blue')
        # axA[0,0].set_ylim([-10,5])
        for value in sac_data[metal]:
            vacancy, dopant = value.split('_')
            index_v = int(vacancy)-1
            index_d = int(dopant)-1
            pdos_slab = []
            pdos_ads = []
            for spin in sac_data[metal][value]['pdos']['%s(d)'%metal]:
                pdos_slab.append(sac_data[metal][value]['pdos']['%s(d)'%metal][spin][0])
                pdos_ads.append(sac_data[metal][value]['pdos']['CO2(sp)'][spin][0])

            energies_ads = sac_data[metal][value]['energies_ads']
            energies_slab = sac_data[metal][value]['energies_slab']

            all_hilbert = []
            all_f = []
            for i in range(len(pdos_slab)):
                f = get_filled_states(energies_slab, pdos_slab[i])
                all_f.append(f)
            summed_dos_slab = np.array(pdos_slab).sum(axis=0)
            summed_dos_ads = np.array(pdos_ads).sum(axis=0)
            # hilbert_t = np.imag(signal.hilbert(pdos_slab[i]))
            # ax[index_v,index_d].fill_between((-1)**i * np.array(pdos_slab[i]), energies_slab, color='tab:red')
            ax[index_v,index_d].fill_between(summed_dos_slab, energies_slab, color=jmol_colors[atomic_numbers[metal]])
            ax[index_v,index_d].set_title('vac:%s dop:%s'%(vacancy, dopant))
            # axA[index_v,index_d].fill_between((-1)**i * np.array(pdos_ads[i]), energies_ads, color='tab:blue')
            axA[index_v,index_d].plot( summed_dos_ads, energies_ads, color='tab:green')
            # axA[index_v,index_d].fill_between((-1)**i * np.array(pdos_slab[i]), energies_slab, color='tab:red')
            axA[index_v,index_d].fill_between(summed_dos_slab, energies_slab, color=jmol_colors[atomic_numbers[metal]])
            axA[index_v,index_d].set_title('vac:%s dop:%s'%(vacancy, dopant))
            axA[index_v,index_d].set_ylim([-10,5])
            ax[index_v,index_d].set_ylim([-10,5])
            axA[index_v,index_d].set_xticks([])
            filled_indices = [i for i in range(len(energies_slab)) if energies_slab[i] < 0]
            all_hilbert.append(energies_slab[np.argmax(hilbert_t[filled_indices])])
            # ax[index_v,index_d].plot(hilbert_t, energies, color='tab:green')
            # dipole_desc.setdefault('%s_%s'%(vacancy, dopant),{})[metal] = energies[np.argmax(hilbert_t)]
            # dipole_desc.setdefault('%s_%s'%(vacancy, dopant),{})[metal] = all_hilbert[0] * all_f[0] + all_hilbert[1] * all_f[1]

        fig.tight_layout()
        fig.savefig('output/%s_pdos_sac.pdf'%metal)
        figA.tight_layout()
        figA.savefig('output/ads_%s_pdos_sac.pdf'%metal)

    return dipole_desc, ideal_params

def estructure_model(dipole_desc):

    qexp_dict = '../databases/explicit_charge.json'
    E0_dict = '../databases/zero_charge_energies.json'
    tmdbname = '../databases/transition_metal_vacuum.db'
    sacdbname = '../databases/single_atom_vacuum.db' 
    
    results = {}
    parsedb(results, connect(tmdbname), dbconfig={'states':'state_CO2_dos'} )
    # parsedb(results, connect(sacdbname), sac=True, dbconfig={'states':'state_CO2_dos'})
    pprint(results)

    qdict = json.load(open(qexp_dict, 'r'))
    E0dict = json.load(open(E0_dict, 'r'))
    colors = plt.get_cmap('viridis', len(qdict))

    fig, ax = plt.subplots(1, 2, figsize=(12,6))
    # pprint(qdict)

    all_points = []
    all_E_points = []

    for i, facet in enumerate(qdict):
        for metal in qdict[facet]:
            # if '_3' in facet and metal == 'Fe':
                # continue
            if '_' in facet:
                continue
            try:
                q = dipole_desc[facet][metal]
                mu =  qdict[facet][metal]#results[facet][metal]['dipole']#qdict[facet][metal] #
                E0 = E0dict[facet][metal]
                # magmom = results[facet][metal]['magmom']
                ax[0].plot(q, mu, 'o', color=colors(i))
                ax[1].plot(q, E0, 'o', color=colors(i))
                # ax[0].errorbar(q, mu, 0.025, color=colors(i))
                ax[0].annotate(r'%s %s'%(metal,facet), xy=(q, mu))
                ax[1].annotate(r'%s %s'%(metal,facet), xy=(q, E0))
                all_points.append([q, mu])

                if metal != 'Ag':
                    all_E_points.append([q, E0])
            except KeyError:
                continue

    all_points = np.array(all_points)
    all_E_points = np.array(all_E_points)
    eu_all, q_all = all_points.transpose()
    eu_all_E, E_all = all_E_points.transpose()
    fit = get_fit_from_points(eu_all, q_all, 1)
    fit_E = get_fit_from_points(eu_all_E, E_all, 1)
    ax[0].plot(eu_all, fit['p'](eu_all), '-', lw=2 , color='k')
    ax[1].plot(eu_all_E, fit_E['p'](eu_all_E), '-', lw=2, color='k')
    ax[0].set_ylabel(r'$\mu \left ( CO_2^* \right ) $ / $e \AA$')
    # ax[0].set_ylabel(r'$q \left ( \textrm{CO}_2^* \right ) $ / $e \AA$')
    ax[1].set_ylabel(r'$\Delta E \left ( CO_2^* \right ) $ / $eV$')
    ax[0].set_xlabel(r"$W_{d}$ / eV")
    ax[1].set_xlabel(r"$\epsilon_{u}$ / eV")
    fig.tight_layout()
    fig.savefig('output/estructure.pdf')


if __name__ == "__main__":
    get_plot_params()
    create_output_directory()


    dipole_desc, ideal_params = plot_dos()
    schematic(ideal_params)
    print(dipole_desc)

    estructure_model(dipole_desc)
    

