from os import path
from aiida import cmdline, engine
from aiida.plugins import DataFactory, CalculationFactory
import click
from aiida.orm import SinglefileData, List, Dict, Int, Float, Str

def run_calculation(facet, pH, energy_pk, label):
    """Run AiiDA CatMAP calculation

    :param catmap_code: Code object for CatMAP
    :type catmap_code: Code
    :param potential: SHE potential
    :type potential: float
    :param pH: pH value
    :type pH: float
    :param energy_pk: pk of the energy file
    :type energy_pk: int
    :param label: Label to denote the computation
    :type label: str
    """
    catmap_code = load_code('catmap-sigma@localhost')

    energies = load_node(energy_pk)
    potential = energies.get_extra('potential_SHE')

    species_definitions = {}
    species_definitions['CO2_g'] = {'pressure':0.2} 
    species_definitions['CO_g'] = {'pressure':0.1}
    species_definitions['H2_g'] = {'pressure':1}
    species_definitions['H2O_g'] = {'pressure':1}
    species_definitions['s'] = {'site_names': [facet], 'total':1}

    scaling_constraint_dict = {
                           'COOH_s':['+','+',None],
                           'CO2_s':['+','+',None],
                           'CO_s':['+', '+', None],
                           }

    ideal_gas_params = { 
                         'CO2_g' : [2,'linear', 0],
                         'CO_g' : [1,'linear', 0],
                         'H2_g': [2,'linear', 0],
                         'H2O_g': [2,'nonlinear', 0],
                        }

    surfaces = ['Pt', 'Pd', 'Cu', 'Ag', 'Au']

    inputs = {
        'code': catmap_code,
        'energies': energies,
        'rxn_expressions':List(list=[
                        'CO2_g + *_s <-> ^0.01eV_s <-> CO2_s',
                        'CO2_s + + H_g + ele_g <-> COOH_s',
                        'COOH_s + H_g + ele_g <-> CO_s + H2O_g', 
                        'CO_s <-> CO_g + *_s',
                        'H2_g <-> H2_g',
        ]), 
        'surface_names':List(list=surfaces), 
        'descriptor_names':List(list=['COOH_s','CO2_s']), 
        'descriptor_ranges':List(list=[[-2.5, 1.5 ], [-2.5, 1.5]]), 
        'resolution':Int(50), 
        'voltage':Float(potential),
        'pH': Float(pH),
        'temperature':Float(300), 
        'species_definitions':Dict(dict=species_definitions), 
        'gas_thermo_mode':Str('ideal_gas'), 
        'adsorbate_thermo_mode':Str('harmonic_adsorbate'),
        # Note: here we have the surface charge correction pre-included
        'electrochemical_thermo_mode':List(list=['simple_electrochemical']),
        'scaling_constraint_dict':Dict(dict=scaling_constraint_dict), 
        'decimal_precision':Int(150), 
        'tolerance':Float(1e-20), 
        'max_rootfinding_iterations':Int(100), 
        'max_bisections':Int(3), 
        'ideal_gas_params':Dict(dict=ideal_gas_params),
        'metadata': {
            'label':label%potential,
            'description': "CatMAP calculation for CO2R",
        },
    }

    future = engine.submit(CalculationFactory('catmap'), **inputs)
    grp = Group.get(label='kinetic_models/descriptors_CO2_COOH')
    grp.add_nodes(future)


def main():
    """Main function to run the calculation."""
    FACET = '211' # Use the 211 facet to define the scaling line
    PH = 2 # The pH value at which the experiments were done
    ENERGY_PK = 403 # pk of the energy file
    LABEL = 'CatMAP calculation at potential: %1.2f' 

    run_calculation(FACET, PH, ENERGY_PK, LABEL)


if __name__ == '__main__':
    main() 
