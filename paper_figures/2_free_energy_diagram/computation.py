
# from cruchparser.forces.findiff_adsorbates import ForceExtrapolation
import sys
sys.path.append('../utilities')
from findiff import ForceExtrapolation
import numpy as np
from ase.db import connect 
from useful_functions import get_vasp_nelect0
from pprint import pprint
from useful_functions import get_reference_energies
from ase import units
from useful_functions import get_fit_from_points
import matplotlib.pyplot as plt
from ase.thermochemistry import HarmonicThermo, IdealGasThermo
from ase.data import atomic_numbers
from ase.data.colors import jmol_colors
from matplotlib import ticker
import json
from scipy.interpolate import griddata
import traceback
import matplotlib.tri as tri

def gas_phase_energies(references, referencedb):
    for row in referencedb.select():
        functional = row.functional
        pw = row.pw
        state = row.states.replace('state_','')
        references.setdefault(functional,{}).setdefault(pw,{}).setdefault(state,{})['energy'] = row.energy
        references.setdefault(functional,{}).setdefault(pw,{}).setdefault(state,{})['atoms'] = row.toatoms()

def parse_data(database, results):
    """Parse keys from the database

    :param database: database
    :type database: db
    """

    for row in database.select():

        try:
            metal = row.sampling.replace('sampling_','')
            type_of_calc = 'TM'
        except AttributeError:
            metal = row.metal_dopant.replace('metal_dopant_','').replace('_nonorth','')
            type_of_calc = 'SAC'
            # In the case of Fe take only calculations with a U
            if 'Fe' in row.metal_dopant:
                try:
                    hubbard_U = [a for a in row.data['ldau']['U'] if a > 0.0 ][0]
                except KeyError:
                    continue
        try:
            displacement = row.displacement
            findiff_calc = True
        except AttributeError:
            findiff_calc = False
        
        state = row.states.replace('state_','').replace('implicit_','')
        implicit = 'implicit' if row.implicit else 'vacuum'
        tot_charge = row.tot_charge
        charge = tot_charge - get_vasp_nelect0(row.toatoms())

        if type_of_calc == 'TM':
            facet = row.facets.replace('facet_','')
        elif type_of_calc == 'SAC':
            facet = row.vacancy_number.replace('vacancy_','') + '_' + row.dopant_number.replace('dopant_','')

        if not findiff_calc:
            try:
                results.setdefault(facet,{}).setdefault(metal,{}).setdefault(state,{})\
                    .setdefault(implicit,{}).setdefault(charge,{})['energy'] = row.energy
            except AttributeError:
                continue
            results.setdefault(facet,{}).setdefault(metal,{}).setdefault(state,{})\
                .setdefault(implicit,{}).setdefault(charge,{})['atoms'] = row.toatoms()
            try:
                magmom = row.toatoms().get_magnetic_moment()
            except:
                magmom = 0.0
            results.setdefault(facet,{}).setdefault(metal,{}).setdefault(state,{})\
                .setdefault(implicit,{}).setdefault(charge,{})['magmom']  = magmom

            if implicit == 'vacuum':
                results.setdefault(facet,{}).setdefault(metal,{}).setdefault(state,{})\
                    .setdefault(implicit,{}).setdefault(charge,{})['dipole'] = row.dipole_field
                results.setdefault(facet,{}).setdefault(metal,{}).setdefault(state,{})\
                    .setdefault(implicit,{}).setdefault(charge,{})['wf'] = row.wf

                # get the area by multiplying lattice vectors
                atoms = row.toatoms()
                cell = atoms.get_cell()
                area = atoms.get_volume() / cell[-1, -1] * 1e-16
                results.setdefault(facet,{}).setdefault(metal,{}).setdefault(state,{})\
                    .setdefault(implicit,{}).setdefault(charge,{})['area'] = area

        if findiff_calc:

            ## now add in details for the finite difference
            findiff = row.findiff
            indices = int(findiff[:-2])
            field = row.field
            direction = findiff[-2]
            displacement = float(displacement.split('_')[-2])
            results.setdefault(facet,{}).setdefault(metal,{}).setdefault(state,{})\
                .setdefault(implicit,{}).setdefault(charge,{}).setdefault('findiff',{})\
                    .setdefault(indices,{}).setdefault(direction,{})\
                    .setdefault(displacement,{}).setdefault(field,{})['dipole'] = row.dipole_field
            results.setdefault(facet,{}).setdefault(metal,{}).setdefault(state,{})\
                .setdefault(implicit,{}).setdefault(charge,{}).setdefault('findiff',{})\
                    .setdefault(indices,{}).setdefault(direction,{})\
                    .setdefault(displacement,{}).setdefault(field,{})['forces'] = row.toatoms().get_forces()[indices]


def get_frequencies():

    frequencies = {}
    frequencies['CO'] = [1832.373539, 481.555294, 467.482512, 425.061714, 73.09318, 68.233697]
    frequencies['COOH'] = [3579.438347, 1573.100374, 1241.519251, 932.338544, 624.650858, 619.94479, 468.883838, 246.531351, 245.478053, 75.740879, 73.602391, 15.378853]
    frequencies['CO2'] = [1846.027301, 1179.574916, 546.070088, 534.269426, 168.712357, 136.301896, 66.113421, 51.400916, 36.766228]

    frequencies['COg'] = [2102.17, 30.57, 30.53]
    frequencies['CO2g'] = [2337.10, 1303.43, 626.20, 626.20, 8.98, 8.98 ]
    frequencies['H2g'] = [4357.74, 101.87, 101.851]
    frequencies['H2Og'] = [3823.99, 3715.40, 1599.34, 84.90, 76.593, 10 ]

    return frequencies


def create_reference_dict(referencedb, frequencies):
    """
    Create reference dictionary with the input
    The references are made similar to how CatMAP requires them to be made
    Eg. for CO2 reduction,
    H2 = 0, CO2 = 0, H2O = 0
    and CO is referenced accoring to the water gas shift reaction

    :param referencedb: Reference database with all gas references
    :type referencedb: ASE database
    """

    references = {}
    gas_phase_energies(references, referencedb=referencedb)

    frequencies = get_frequencies()

    # CatMAP references
    functional = 'RP'
    CO2g_correction = 0.45

    # store internal energies as fixed variables
    COg_E = references[functional][500.]['CO']['energy']
    CO2g_E = references[functional][500.]['CO2']['energy'] + CO2g_correction
    H2g_E = references[functional][500.]['H2']['energy']
    H2Og_E = references[functional][500.]['H2O']['energy']

    
    ## Corrections to the free energy 
    cmtoeV = 0.00012 
    COg_G = IdealGasThermo(
        cmtoeV * np.array(frequencies['COg']),
        geometry='linear', 
        atoms=references[functional][500.]['CO']['atoms'],
        symmetrynumber=2,
        spin=0,
    ).get_gibbs_energy(298.15, 101325, verbose=False)

    CO2g_G = IdealGasThermo(
        cmtoeV * np.array(frequencies['CO2g']),
        geometry='linear', 
        atoms=references[functional][500.]['CO2']['atoms'],
        symmetrynumber=2,
        spin=0,
    ).get_gibbs_energy(298.15, 101325, verbose=False) 

    H2g_G = IdealGasThermo(
        cmtoeV * np.array(frequencies['H2g']),
        geometry='linear', 
        atoms=references[functional][500.]['H2']['atoms'],
        symmetrynumber=2,
        spin=0,
    ).get_gibbs_energy(298.15, 101325, verbose=False)

    H2Og_G = IdealGasThermo(
        cmtoeV * np.array(frequencies['H2Og']),
        geometry='linear', 
        atoms=references[functional][500.]['H2O']['atoms'],
        symmetrynumber=3,
        spin=0,
    ).get_gibbs_energy(298.15, 101325, verbose=False)

    # Internal energy references
    reference_energies_E = {}
    reference_energies_E['CO'] = CO2g_E + H2g_E - H2Og_E 
    reference_energies_E['COOH'] = CO2g_E + 0.5 * H2g_E 
    reference_energies_E['CO2'] =  CO2g_E
    reference_energies_E['CO(g)'] = COg_E + H2Og_E - H2g_E - CO2g_E 

    gas_dict = {'CO2':0.0,
                'CO':reference_energies_E['CO(g)'],
                'H2':0.0,
                'H2O':0.0}

    writeout_gas = []
    for gas_name, gas_E in gas_dict.items():
        writ = ['None', 'gas', gas_name, round(gas_E, 3),  frequencies[gas_name+'g'], 'sv_calc']
        writeout_gas.append(writ)
    
    # Free energy references
    reference_energies = {}
    reference_energies['CO'] = CO2g_G - H2Og_G + H2g_G
    reference_energies['COOH'] = CO2g_G + 0.5 * H2g_G 
    reference_energies['CO2'] = CO2g_G
    reference_energies['CO(g)'] = COg_G + H2Og_G - CO2g_G - H2g_G 
    

    return reference_energies, reference_energies_E, writeout_gas



def get_explicit_charge(vibresults, atomsIS, atomsFS, 
            fields_to_choose=[0.1, 0.2], displacement=0.01, \
            direction='p'):

    method = ForceExtrapolation(
            vibresults=vibresults,
            fields_to_choose=fields_to_choose,
            displacement=displacement,
            direction=direction,
            atomsIS=atomsIS,
            atomsFS=atomsFS,
            )
    
    method.get_dFdG()
    method.get_dmudR()
    method.get_q()

    return method.q[0]


def get_free_energy_diagram_data(databases, referencedb, potential, pH):
    results = {}
    for database in databases:
        parse_data(connect(database), results)
    
    ## Create reference dictionary
    references, references_E, writeout_gas = create_reference_dict(referencedb, get_frequencies())
    colors = {'CO':'tab:blue', 'COOH':'tab:green', 'CO2':'tab:red'}
    U_RHE = potential + 0.059 * pH
    CHE_correction = {'CO':2 * U_RHE, 'COOH':U_RHE}
    frequencies = get_frequencies()
    diagram = {}
    explicit_charge = {}
    E0 = {}

    writeout = []
    writeout_header = ['surface_name', 'site_name', 'species_name', 'formation_energy', 'frequencies', 'reference']
    writeout.append(writeout_header)
    writeout += writeout_gas
    ## Have another writeout for the case where there is no surface charge
    writeout_zero = []
    writeout_zero.append(writeout_header)
    writeout_zero += writeout_gas

    for facet in results:
        for metal in results[facet]:
            if metal in ['Ni', 'Al'] and facet == '111':
                continue

            ## verification plot
            fig, ax = plt.subplots(1, 1, figsize=(6,4))
            figs, axs = plt.subplots(1,1,figsize=(6,4))
            for state in results[facet][metal]:
                # No need to parse data for the slab
                if state == 'slab':
                    continue
                if 'sp' in state or 'gas' in state or 'dos' in state:
                    continue
                if state == 'CO2':
                    ## get the explicit charge
                    try:
                        vibresults = results[facet][metal][state]['vacuum'][0.0]['findiff']
                        q_eff = get_explicit_charge(
                            vibresults=vibresults,
                            atomsFS=results[facet][metal][state]['implicit'][2.00]['atoms'],
                            atomsIS=results[facet][metal]['CO2_gas']['vacuum'][0.0]['atoms'],
                                                )
                        explicit_charge.setdefault(facet,{})[metal] = q_eff
                    except KeyError as e:
                        print(metal, facet)
                        pprint(results[facet][metal][state])
                        traceback.print_exc()
                        continue
                else:
                    # Another non CO2 adsorbate
                    q_eff = 0
                
                Eq = []
                q = []
                sigma = []
                magmoms = []

                try:
                    for charge in results[facet][metal][state]['implicit']:
                        try:
                            dE = results[facet][metal][state]['implicit'][charge]['energy'] \
                                -results[facet][metal]['slab']['implicit'][charge]['energy'] \
                                -references_E[state]
                        except KeyError:
                            continue
                        Eq.append(dE)
                        q.append(charge - q_eff/2)
                        area = results[facet][metal]['slab']['vacuum'][0.0]['area'] 
                        sigma = q / area * units._e * 1e6 # mu C / cm-2 
                        magmoms.append(np.abs(results[facet][metal][state]['implicit'][charge]['magmom']))
                except KeyError:
                    continue
                
                ## get the fit of the energy vs surface charge
                fit = get_fit_from_points(sigma, Eq, 1) # linear fit of energy to surface charge
                ax.plot(sigma, Eq, 'o', color=colors[state])
                ax.plot(sigma, fit['p'](sigma), color=colors[state])
                
                axs.plot(sigma, magmoms, 'o-', color=colors[state])

                ## surface charge corresponding to the requested potential
                ## Assume pzc is the same as wf for now
                ## TODO: Validate the assumption - check sensitivity
                if metal not in ['Fe', 'Ni']:
                    pzc = results[facet][metal]['slab']['vacuum'][0.0]['wf'] - 4.4
                else:
                    ## doped metal changes the wf a lot 
                    pzc = -0.05 # V vs SHE
                C_gap = 25 # mu F cm-2
                sigma_for_pot = C_gap * ( potential - pzc )
                # Save the results needed for the diagram
                ## Get the vibrational energies 
                dG_correct = HarmonicThermo(0.00012 * np.array(frequencies[state])).get_helmholtz_energy(298.15, verbose=False)
                if 'CO2' in state:
                    # Then we have a relationship with surface charge ready
                    diagram.setdefault(facet,{}).setdefault(metal,{})[state] = fit['p'](-1*sigma_for_pot) + dG_correct - references[state] 
                    correction = + dG_correct - references[state]  
                    # print('CO2 correction:%1.3f'%correction)
                    writeout.append([metal, facet, state, round(fit['p'](-1*sigma_for_pot),2), frequencies[state], 'sv_calc'])
                    writeout_zero.append([metal, facet, state, round(fit['p'](0),2), frequencies[state], 'sv_calc'])
                    E0.setdefault(facet,{})[metal] = fit['p'](0) - references[state]
                else:
                    dE0 = fit['p'](0)
                    dE = fit['p'](-1*sigma_for_pot) #mu * 1e-6 * sigma_for_pot / units._eps0 / eps_r * 1e-6 
                    diagram.setdefault(facet,{}).setdefault(metal,{})[state] = dE + CHE_correction[state] + dG_correct - references[state] 
                    correction = + dG_correct - references[state]  
                    # print('%s correction:%1.3f'%(state,correction))
                    writeout.append([metal, facet, state, round(dE,2), frequencies[state], 'sv_calc'])
                    writeout_zero.append([metal, facet, state, round(dE0,2), frequencies[state], 'sv_calc'])

            ## Add in the results for the gas phase species
            diagram.setdefault(facet,{}).setdefault(metal,{})['CO2(g)'] = 0
            diagram.setdefault(facet,{}).setdefault(metal,{})['CO(g)'] = references_E['CO(g)'] + references['CO(g)'] + CHE_correction['CO'] 

            ax.set_ylabel(r'$\Delta E$ / eV')
            ax.set_xlabel(r'$\sigma$ / $\mu C cm^{-2}$')

            axs.set_ylabel(r'Magnetic Moment')
            axs.set_xlabel(r'$\sigma$ / $\mu C cm^{-2}$')
            # axs.set_xlabel(r'Potential / V vs. SHE')

            for i, j in colors.items():
                ax.plot([],[], color=j, label=r''+i.replace('2','$_{2}$'))

            for i, j in colors.items():
                axs.plot([],[], color=j, label=r''+i.replace('2','$_{2}$'))

            ax.legend(loc='best', frameon=False, fontsize=12)
            axs.legend(loc='best', frameon=False, fontsize=12)
            fig.tight_layout()
            figs.tight_layout()
            fig.savefig('output/SI_charging_curve_metal_%s_facet_%s.pdf'%(metal, facet))
            figs.savefig('output/SI_magmom_curve_metal_%s_facet_%s.pdf'%(metal, facet))
            plt.close(fig)
            plt.close(figs)
    

    return diagram, writeout, writeout_zero, explicit_charge, E0


def plot_computational_diagram(data, ax):
    ## Plot the CO2 to CO free energy diagram
    all_potentials = []
    ls = ['-', '--', '-.', ':']

    for i, potential in enumerate(data):
        all_potentials.append(float(potential))
        # first the gold plot
        plot_data = [
            data[potential]['100']['Au']['CO2(g)'],
            data[potential]['100']['Au']['CO2'],
            data[potential]['100']['Au']['COOH'],
            data[potential]['100']['Au']['CO'],
            data[potential]['100']['Au']['CO(g)'],
        ]
        fed = [val for val in plot_data for _ in (0, 1)]
        ax[0].plot(fed, color=jmol_colors[atomic_numbers['Au']], lw=3,
                label=r'$%1.1f\mathregular{V}_{\mathregular{SHE}}$'%potential,\
                ls=ls[i])
        
    ax[0].legend(loc='best', frameon=False, \
        fontsize=11,
    )
        # color=jmol_colors[atomic_numbers['Au']])
    
    ## chooses the last potential to plot the rest
    potential = np.max(all_potentials)
    for i, dop in enumerate(range(1,5)):
        try:
            plot_data = [
                data[potential]['2_%s'%dop]['Fe']['CO2(g)'],
                data[potential]['2_%s'%dop]['Fe']['CO2'],
                data[potential]['2_%s'%dop]['Fe']['COOH'],
                data[potential]['2_%s'%dop]['Fe']['CO'],
                data[potential]['2_%s'%dop]['Fe']['CO(g)'],
            ]
            fed = [val for val in plot_data for _ in (0, 1)]
            ax[1].plot(fed, color=jmol_colors[atomic_numbers['Fe']], lw=3,
                    alpha=(i+1) / 4)#, label=r'$%1.1fV_{SHE}$'%potential)
                    # ls=ls[i])

            ax[1].annotate(r'FeN$_{%d}$'%dop, xy=(0.8,0.87-0.15*i), alpha=(i+1)/4, \
                        color=jmol_colors[atomic_numbers['Fe']], xycoords='axes fraction', \
                            fontsize=12)
        except KeyError:
            # continue
            pass
        # and now the NiNC plot
        plot_data = [
            data[potential]['2_%s'%dop]['Ni']['CO2(g)'],
            data[potential]['2_%s'%dop]['Ni']['CO2'],
            data[potential]['2_%s'%dop]['Ni']['COOH'],
            data[potential]['2_%s'%dop]['Ni']['CO'],
            data[potential]['2_%s'%dop]['Ni']['CO(g)'],
        ]
        fed = [val for val in plot_data for _ in (0, 1)]
        ax[2].plot(fed, color=jmol_colors[atomic_numbers['Ni']], lw=3,
                alpha=(i+2) / 5,)# label=r'$%1.1fV_{SHE}$'%potential)
                # ls=ls[i])
        ax[2].annotate(r'NiN$_{%d}$'%dop, xy=(0.8,0.87-0.15*i), alpha=(i+1)/4, \
                    color=jmol_colors[atomic_numbers['Ni']], xycoords='axes fraction', \
                        fontsize=12)
        
        ax[1].set_ylabel(r'$\Delta G$ / eV')

        ax[0].set_xticks([])
        ax[1].set_xticks([])
        ax[2].set_xticks([])
        ax[0].set_ylim([-1.6,2.2])
        ax[1].set_ylim([-1.6,2.2])
        ax[2].set_ylim([-1.6,2.2])

        ax[2].annotate(r'CO$_2$(g)', xy=(0.05, -0.2), xycoords='axes fraction', fontsize=14)
        ax[2].annotate(r'CO$_2^*$', xy=(0.25, -0.2), xycoords='axes fraction', fontsize=14)
        ax[2].annotate(r'COOH$^*$', xy=(0.4, -0.2), xycoords='axes fraction', fontsize=14)
        ax[2].annotate(r'CO$^*$', xy=(0.65, -0.2), xycoords='axes fraction', fontsize=14)
        ax[2].annotate(r'CO(g)', xy=(0.85, -0.2), xycoords='axes fraction', fontsize=14)

        ax[0].annotate(r'$\Delta G_{\mathregular{CO}_2^*} > \Delta G_{\mathregular{COOH}^*}$', color='k', 
                    xy=(0.72, -1),fontsize=14 )
        ax[1].annotate(r'$\Delta G_{\mathregular{CO}_2^*} > \Delta G_{\mathregular{COOH}^*}$', color='k', 
                    xy=(0.72, -1),fontsize=14 )
        ax[2].annotate(r'$\Delta G_{\mathregular{CO}_2^*} < \Delta G_{\mathregular{COOH}^*}$', color='k', 
                    xy=(0.72, -1),fontsize=14)
        
def plot_variation_with_potential(explicit_charge):
    fig, ax = plt.subplots(1, 4, figsize=(20,5), sharex=True, sharey=True)
    potential_range = np.linspace(-0.5, 0.5)

    for facet in explicit_charge:
        for metal in explicit_charge[facet]:
            fit = np.poly1d([explicit_charge[facet][metal], 0])

            if facet == '100':
                i = 0
                alpha = 0.75
                label = metal + '(100)'
                annotation = 'Transition Metal (100)'
            elif facet == '211':
                i = 1
                alpha = 0.75
                label = metal + '(211)'
                annotation = 'Transition Metal (211)'
            elif '2_' in facet:
                i = 3
                alpha = 1 / int(facet.split('_')[-1])
                label = r'$'+metal + 'N_%s'%facet.split('_')[-1] +'$'
                annotation = 'Double Vacancy'
            
            elif '1_' in facet:
                i = 2
                alpha = 1 / int(facet.split('_')[-1])
                label = r'$'+metal + 'N_%s'%facet.split('_')[-1] +'$'
                annotation = 'Single Vacancy'

            ax[i].plot(potential_range, fit(-1*potential_range), '-', lw=3,\
                 color=jmol_colors[atomic_numbers[metal]], alpha=alpha, label=label)
            ax[i].annotate(annotation, xy=(0.4, 0.2), color='k', xycoords='axes fraction')
    
    for a in ax:
        a.set_xlabel(r'$\phi - \phi_{pzc}$ / V ')
        a.legend(loc='best', frameon=False, fontsize=14)

    ax[0].set_ylabel(r'$\Delta E_{CO_{2}^*} - \Delta E_{CO_{2}^*}^0 $ / eV')

    fig.tight_layout()

    return fig, ax



                
        



                
