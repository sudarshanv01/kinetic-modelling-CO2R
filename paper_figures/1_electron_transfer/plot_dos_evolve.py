


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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Circle
from ase.data import covalent_radii as radii
import string
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

def parsedb(database, results):
    for row in database.select():
        results.setdefault(row.states.replace('state_',''),{})['pdos'] = row.data.pdos.pdos
        results.setdefault(row.states.replace('state_',''),{})['energy'] = row.data.pdos.energies
        results.setdefault(row.states.replace('state_',''),{})['integrated_dos'] = row.data.pdos.integrated_dos
        results.setdefault(row.states.replace('state_',''),{})['total_dos'] = row.data.pdos.total_dos
        results.setdefault(row.states.replace('state_',''),{})['wf'] = row.wf
        results.setdefault(row.states.replace('state_',''),{})['atoms'] = row.toatoms()
        results.setdefault(row.states.replace('state_',''),{})['magmom'] = row.magmom

def lorentz_dos(a, b1,energy):
    delta = 1 / np.pi * ( a / ( (energy - b1)**2 + a**2  ) )
    return delta

@click.command()
@click.option('--mncdb', default='databases/single_atom_rls.db')
@click.option('--npdb', default='databases/nanoparticle_rls.db')
def main(mncdb, npdb):

    results = {}
    results['mnc'] = {} ; results['np']  = {}
    parsedb(connect(mncdb), results['mnc']) 
    # parsedb(connect(npdb), results['np'])

    hbar = 6.582 * 1e-16 # eV.s

    # fig, ax = plt.subplots(2, 5, figsize=(12,7), squeeze=False, constrained_layout=True)
    fig = plt.figure(constrained_layout=True, figsize=(12, 10))
    gs = fig.add_gridspec(2,10, wspace=0.05)
    ax = []
    for i in range(1):
        temp_ax = []
        for j in range(5):
            temp_ax.append(fig.add_subplot(gs[i+1,2*j:2*j+2]))
        ax.append(temp_ax)
    ax = np.array(ax)
    axl = fig.add_subplot(gs[0,5:])
    axr = fig.add_subplot(gs[0,0:5])
    # axf = fig.add_subplot(gs[2,3:])

    ## plot lorentzian data
    energy_range = np.linspace(-3,3,500)
    spec_lorentz = [[0.01,0], [0.1,0],[0.5,0], [1,0], [2,0]]
    color = ['tab:blue','tab:red', 'tab:green', 'tab:orange', 'tab:purple']
    for m, spec in enumerate(spec_lorentz):
        peak = lorentz_dos(*spec, energy_range)
        rate = 2 * np.pi / hbar * spec[0]
        axl.plot(peak, energy_range, '-', color=color[m])
        axl.fill_between(peak, energy_range, color=color[m], alpha=0.25)
        axr.plot(spec[0], rate, 'o', color=color[m], )
    axl.set_ylabel(r'Energy / eV')
    axl.set_xlim([0,1])
    axr.axhline(1e12, color='k', ls='--')
    axr.annotate('Diffusion rate', xy=(0.5,0.12), xycoords='axes fraction')
    axr.set_yscale('log')
    axl.set_xticks([])
    axr.set_ylabel(r'Rate / s$^{-1}$')
    axr.set_xlabel(r'Width / eV')
    axl.annotate(r'$\Delta = \Sigma_{k} V_{ak}^2 \delta \left ( \epsilon - \epsilon_{k} \right )$', \
                    xy=(0.4,0.8), xycoords='axes fraction')


    plot_states = {}
    plot_states['mnc'] = ['00','04','07', '08', '09' ]
    plot_states['np'] = ['00','01','02', '03', '04']
    width = [[-0.1,0.1], [-0.45,0.3], [-0.8,0.3]]

    # for i, conc in enumerate(sorted(results)):
    for i, species in enumerate(results):
        j=0
        n=0
        all_f = []
        all_E = []
        for state in results[species]:
            energy = np.array(results[species][state]['energy'])
            filled_indices = [a for a in range(len(energy)) if energy[a] <= 0.0]
            pdos_Fe = np.array(results[species][state]['pdos']['Fe']['+']) \
                    +np.array(results[species][state]['pdos']['Fe']['-'])
            pdos_co2 = np.array(results[species][state]['pdos']['co2']['+']) \
                    +np.array(results[species][state]['pdos']['co2']['-']) 
            if species == 'mnc':
                pdos_C = np.array(results[species][state]['pdos']['C']['+']) \
                        +np.array(results[species][state]['pdos']['C']['-'])
                pdos_N = np.array(results[species][state]['pdos']['N']['+']) \
                        +np.array(results[species][state]['pdos']['N']['-']) 
            co2 = pdos_co2[0].sum(axis=0)
            f = np.trapz(co2[filled_indices], energy[filled_indices]) / np.trapz(co2, energy)
            # if species == 'mnc':
            all_f.append(f)
            if state in plot_states[species]:
                ax[i,j].plot(co2, energy, color='tab:blue', alpha=0.5)
                ax[i,j].fill_between(co2, energy, color='tab:blue', alpha=0.1)
                ax[i,j].set_ylim([-8,5])
                ax[i,j].set_xlim([0.0, 0.1])
                axins = ax[i,j].inset_axes([0.5, 0.5, 0.47, 0.47])
                axins.plot(co2, energy, color='tab:blue')
                axins.fill_between(co2, energy, color='tab:blue', alpha=0.25)
                axins.set_xlim([0.0, 0.006])
                axins.set_ylim([-2,2])
                # ax[i,j].indicate_inset_zoom(axins,  ec='r', lw=3, fc='none')
                mark_inset(ax[i,j], axins, loc1=2, loc2=3, fc="none", lw=3, ec='k')

                axins.set_xticklabels('')
                axins.set_yticklabels('')
                if j== 0:
                    ax[i,j].set_ylabel(r'$\epsilon - \epsilon_{f}$ / eV')
                else:
                    ax[i,j].set_yticks([])
                ax[i,j].set_xticks([])
                ax[i,j].axhline(0, ls='--')
                if j>1 and species == 'mnc':
                    # ax[i,j].annotate(r'$\Delta > 0.1$ eV', xy=(0.05,0.40), color='tab:red', xycoords='axes fraction' )
                    ax[i,j].axhspan(*width[n], color='tab:red', alpha=0.5)
                    axins.axhspan(*width[n], color='tab:red', alpha=0.5)
                    n += 1

                ## plot the atoms as insets
                atoms = results[species][state]['atoms'] 
                axins4 = inset_axes(ax[i,j], width="40%", height="40%", loc=4, borderpad=1)
                for atom in atoms:
                    if species == 'mnc':
                        if atom.index not in [26, 32, 33, 27, 31, 30, 28, 29, 18, 24, 23, 11, 17]:
                            continue

                    color = jmol_colors[atom.number]
                    radius = radii[atom.number]
                    if species == 'mnc':
                        circle = Circle((atom.y, atom.z), radius, facecolor=color,
                                                edgecolor='k', linewidth=1)
                    else:
                        circle = Circle((atom.x, atom.y), radius, facecolor=color,
                                                edgecolor='k', linewidth=1)
                    axins4.add_patch(circle)

                # Enforce the circles are round (equal x and y scales) and turn off
                # tickmarks.
                axins4.axis('equal')
                axins4.set_xticks([])
                axins4.set_yticks([])
                axins4.axis('off')


                j+=1

            ax[0,2].annotate('TS', xy=(0.25,0.1), xycoords='axes fraction')
            # ax[1,2].annotate('TS', xy=(0.25,0.1), xycoords='axes fraction')
    
        # axf.plot(all_f, 'o--', color='tab:blue')
    fig.suptitle(r'Rate of e$^-$ transfer on MNC $ \approx 10^{14} \mathregular{s}^{-1}$')
    fig.text(0.6, 0.46, '(s,p) PDOS / arb. units.', ha='center')
    # fig.set_constrained_layout_pads(w_pad=2/72, h_pad=2/72, hspace=0.2, wspace=0.2)
    alphabet = list(string.ascii_lowercase)
    for i, a in enumerate([axr, axl] + ax.flatten().tolist()):
        a.annotate(alphabet[i]+')', xy=(0.1, 0.9), xycoords='axes fraction', fontsize=18)
    fig.savefig('output/dos.pdf')

if __name__ == "__main__":
    create_output_directory()
    get_plot_params()
    main()
