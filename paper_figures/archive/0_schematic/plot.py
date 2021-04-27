
import numpy as np
import matplotlib.pyplot as plt
from useful_functions import create_output_directory
from plot_params import get_plot_params

def main():
    co2 = [0, 0.75, 0.5, 0.30, 0.25]
    cooh = [0, 0.75, 1., 0.30, 0.25]
    co2_limiting = [val for val in co2 for _ in (0, 1)]
    cooh_limiting = [val for val in cooh for _ in (0,1)]
    
    fig, ax = plt.subplots(1, 1, figsize=(8,6))

    ax.plot(co2_limiting, '-', lw=4, color='tab:blue', label=r'$CO_2^*$ limiting')
    ax.plot(cooh_limiting, '-', lw=4, color='tab:red', label=r'$COOH^*$ limiting')

    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_ylabel(r'$\Delta G$ / eV')
    ax.annotate(r'$CO_2(g)$', xy=(0.05, -0.1), xycoords='axes fraction')
    ax.annotate(r'$CO_2^*$', xy=(0.25, -0.1), xycoords='axes fraction')
    ax.annotate(r'$COOH^*$', xy=(0.4, -0.1), xycoords='axes fraction')
    ax.annotate(r'$CO^*$', xy=(0.65, -0.1), xycoords='axes fraction')
    ax.annotate(r'$CO(g)$', xy=(0.85, -0.1), xycoords='axes fraction')

    ax.annotate(r'$A^{-1} = - \frac{1}{k_BT} \frac{C \mu_{CO_2^*}}{\epsilon}$', color='tab:red', \
         xy=(0.6, 0.9), xycoords='axes fraction', fontsize=12)
    ax.annotate(r'$A^{-1} = - \frac{1}{k_BT} \frac{C \left ( \mu_{COO-H^*} - \mu_{CO_2^*} \right )}{\epsilon}$',\
         color='tab:blue', xy=(0.6, 0.7), xycoords='axes fraction', fontsize=12)

    fig.savefig('output/schematic.pdf')



if __name__ == '__main__':
    create_output_directory()
    get_plot_params()
    main()