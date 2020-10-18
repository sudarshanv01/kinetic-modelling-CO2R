#!/usr/bin/env python
"""Run a test calculation on localhost.

Usage: ./example_01.py
"""
from os import path
from aiida_catmap import helpers
from aiida import cmdline, engine
from aiida.plugins import DataFactory, CalculationFactory
import click
from aiida.orm import SinglefileData, List, Dict, Int, Float, Str

INPUT_DIR = path.join(path.dirname(path.realpath(__file__)), 'input_files')


def test_run(catmap_code, potential, energy_pk, label):
    """Run a calculation on the localhost computer.

    Uses test helpers to create AiiDA Code on the fly.
    """
    if not catmap_code:
        # get code
        computer = helpers.get_computer()
        catmap_code = helpers.get_code(entry_point='catmap', computer=computer)

    # Prepare input parameters

    energies = load_node(energy_pk)

    species_definitions = {}
    species_definitions['CO2_g'] = {'pressure':0.2} 
    species_definitions['CO_g'] = {'pressure':0.1}
    species_definitions['H2_g'] = {'pressure':0.1}
    species_definitions['H2O_g'] = {'pressure':0.1}
    species_definitions['s'] = {'site_names': ['graphene'], 'total':1}

    scaling_constraint_dict = {
                           'CO2_s':[None,None,None],
                           'COOH_s':[None,None,None],
                           #'O-CO_s':'initial_state',
                           #'O-O_s':'final_state',
                           }

    # set up calculation
    # surfaces = ['CoNNC', 'CoNNCC', 'FeNNC', 'FeNNNCC', 'FeNNNNCC', 'NiNNC', 'NiNNCC', 'NiNNNCC', 'MnNCC', 'MnNNCC', 'MnNNNNCC', 'RhNCC', 'RhNNC', 'RhNNCC', 'RhNNNCC', 'RhNNNNCC', 'RuNNC', 'RuNNCC', 'RuNNNCC', 'TiNNC', 'TiNNCC', 'TiNNNC', 'TiNNNCC', 'TiNNNNCC', 'CrNCC', 'CrNNC', 'CrNNCC', 'CrNNNC', 'CrNNNCC', 'CrNNNNCC']
    # surfaces = ['CoNNC', 'CoNNCC', 'CoNNNCC', 'FeNNC', 'FeNNCC', 'FeNNNCC', 'FeCC', 'FeNNNNCC', 'NiNC', 'NiNCC', 'NiNNC', 'NiNNCC', 'NiNNNCC', 'MnNCC', 'MnNNC', 'MnNNCC', 'MnNNNCC', 'MnNNNNCC', 'RhNCC', 'RhNNC', 'RhNNCC', 'RhNNNCC', 'RhNNNNCC', 'RuNNC', 'RuNNCC', 'RuNNNCC', 'RuNNNNCC', 'TiNC', 'TiNNC', 'TiNNCC', 'TiNNNC', 'TiNNNCC', 'TiNNNNCC', 'CrNCC', 'CrNNC', 'CrNNCC', 'CrNNNC', 'CrNNNCC', 'CrNNNNCC']
    # surfaces = ['CoNNC', 'CoNNCC', 'CoNNNCC', 'FeNNC', 'FeNNCC', 'FeNNNCC', 'FeCC', 'FeNNNNCC', 'NiNC', 'NiNCC', 'NiNNC', 'NiNNCC', 'NiNNNCC', 'VNNC', 'MnNCC', 'MnNNC', 'MnNNCC', 'MnNNNCC', 'MnNNNNCC', 'RhNCC', 'RhNNC', 'RhNNCC', 'RhNNNCC', 'RhNNNNCC', 'RuNNC', 'RuNNCC', 'RuNNNCC', 'RuNNNNCC', 'TiNC', 'TiNNC', 'TiNNCC', 'TiNNNC', 'TiNNNCC', 'TiNNNNCC', 'CrNCC', 'CrNNCC', 'CrNNNC', 'CrNNNCC', 'CrNNNNCC']
    # surfaces = ['CoNC', 'CoNNC', 'CoNNCC', 'CoNNNCC', 'FeNNC', 'FeNNCC', 'FeNNNCC', 'FeCC', 'FeNNNNCC', 'NiNC', 'NiNCC', 'NiNNC', 'NiNNCC', 'NiNNNCC', 'VNCC', 'VNNC', 'VNNCC', 'VNNNC', 'VNNNNCC', 'MnNC', 'MnNCC', 'MnNNC', 'MnNNCC', 'MnNNNCC', 'MnNNNNCC', 'RhNC', 'RhNCC', 'RhNNC', 'RhNNCC', 'RhNNNC', 'RhNNNCC', 'RhNNNNCC', 'RuNNC', 'RuNNCC', 'RuNNNC', 'RuNNNCC', 'RuNNNNCC', 'TiNC', 'TiNNC', 'TiNNCC', 'TiNNNC', 'TiNNNCC', 'TiNNNNCC', 'CrNC', 'CrNCC', 'CrNNCC', 'CrNNNC', 'CrNNNCC', 'CrNNNNCC']
    # surfaces = ['CoNC', 'CoNNC', 'CoNNCC', 'FeNNC', 'FeNNCC', 'FeNNNCC', 'FeCC', 'FeNNNNCC', 'NiNC', 'NiNCC', 'NiNNC', 'NiNNCC', 'NiNNNCC', 'VNCC', 'VNNC', 'VNNCC', 'VNNNCC', 'VNNNNCC', 'MnNC', 'MnNCC', 'MnNNC', 'MnNNCC', 'MnNNNCC', 'MnNNNNCC', 'RhNC', 'RhNCC', 'RhNNC', 'RhNNCC', 'RhNNNC', 'RhNNNCC', 'RhNNNNCC', 'RuNNC', 'RuNNCC', 'RuNNNC', 'RuNNNCC', 'RuNNNNCC', 'TiNC', 'TiNNC', 'TiNNCC', 'TiNNNC', 'TiNNNCC', 'TiNNNNCC', 'CrNC', 'CrNCC', 'CrNNC', 'CrNNCC', 'CrNNNC', 'CrNNNCC', 'CrNNNNCC']
    ## Try low spin cases
#    surfaces = ['CoNC', 'CoNNC', 'CoNNCC', 'FeNNC', 'FeNNCC', 'FeCC', 'FeNNNNCC', 'NiNC', 'NiNCC', 'NiNNC', 'NiNNCC', 'NiNNNCC', 'MnNC', 'MnNNC', 'RhNC', 'RhNCC', 'RhNNC', 'RhNNCC', 'RhNNNC', 'RhNNNCC', 'RhNNNNCC', 'RuNNC', 'RuNNCC', 'RuNNNC', 'RuNNNCC', 'RuNNNNCC']
    surfaces = ['CoNC', 'CoNNC', 'CoNNCC', 'FeNNC', 'FeNNCC', 'FeNNNCC', 'FeCC', 'FeNNNNCC', 'NiNC', 'NiNCC', 'NiNNC', 'NiNNCC', 'NiNNNCC', 'VNCC', 'VNNC', 'VNNCC', 'VNNNCC', 'VNNNNCC', 'MnNC', 'MnNCC', 'MnNNC', 'MnNNCC', 'MnNNNCC', 'MnNNNNCC', 'RhNC', 'RhNCC', 'RhNNC', 'RhNNCC', 'RhNNNC', 'RhNNNCC', 'RhNNNNCC', 'RuNNC', 'RuNNCC', 'RuNNNC', 'RuNNNCC', 'RuNNNNCC', 'TiNC', 'TiNNC', 'TiNNCC', 'TiNNNC', 'TiNNNCC', 'TiNNNNCC', 'CrNC', 'CrNCC', 'CrNNC', 'CrNNCC', 'CrNNNC', 'CrNNNCC', 'CrNNNNCC']

    inputs = {
        'code': catmap_code,
        'energies': energies,
        'rxn_expressions':List(list=[
                        'CO2_g + *_s <-> CO2_s',
                        'CO2_s + + H_g + ele_g <-> ^0.01eV_s <-> COOH_s',
                        'COOH_s + H_g + ele_g <-> CO_s + H2O_g', 
                        'CO_s <-> CO_g + *_s',
                        'H2_g <-> H2_g',
        ]), 
        'surface_names':List(list=surfaces), 
        'descriptor_names':List(list=['CO2_s','COOH_s']), 
        'descriptor_ranges':List(list=[[-4, 1], [-4.5, 2]]), 
        'resolution':Int(50), 
        'voltage':Float(potential),
        'temperature':Float(300), 
        'species_definitions':Dict(dict=species_definitions), 
        'gas_thermo_mode':Str('ideal_gas'), 
        'adsorbate_thermo_mode':Str('harmonic_adsorbate'),
        'electrochemical_thermo_mode':List(list=['simple_electrochemical']),
        'scaling_constraint_dict':Dict(dict=scaling_constraint_dict), 
        'decimal_precision':Int(150), 
        'tolerance':Float(1e-20), 
        'max_rootfinding_iterations':Int(100), 
        'max_bisections':Int(3), 
        #'numerical_solver':Str('coverages'),
        'metadata': {
            'label':label, #'Calculation with CO2 assumed solvation and corrected dipole',
            'description': "Submissions with AiiDA CatMAP",
        },
    }

    # Note: in order to submit your calculation to the aiida daemon, do:
    # from aiida.engine import submit
    future = engine.submit(CalculationFactory('catmap'), **inputs)
    grp = Group.get(label='mkm_doped_graphene/test_runs')
    grp.add_nodes(future)
    # result = engine.run(CalculationFactory('catmap'), **inputs)

    # computed_result = result['log'].get_content()


@click.command()
@click.option('--potential')
@click.option('--energy_pk')
@click.option('--label')
@cmdline.utils.decorators.with_dbenv()
@cmdline.params.options.CODE()
def cli(code, potential, energy_pk, label):
    """Run example.

    Example usage: $ ./example_01.py --code diff@localhost

    Alternative (creates diff@localhost-test code): $ ./example_01.py

    Help: $ ./example_01.py --help
    """
    test_run(code, potential, energy_pk, label)


if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter
