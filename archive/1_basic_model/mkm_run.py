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


def test_run(catmap_code):
    """Run a calculation on the localhost computer.

    Uses test helpers to create AiiDA Code on the fly.
    """
    if not catmap_code:
        # get code
        computer = helpers.get_computer()
        catmap_code = helpers.get_code(entry_point='catmap', computer=computer)

    # Prepare input parameters

    energies = load_node(5909)

    species_definitions = {}
    species_definitions['CO2_g'] = {'pressure':0.2} 
    species_definitions['CO_g'] = {'pressure':0.1}
    species_definitions['H2_g'] = {'pressure':0.1}
    species_definitions['H2O_g'] = {'pressure':0.1}
    species_definitions['s'] = {'site_names': ['graphene'], 'total':1}

    scaling_constraint_dict = {
                           'CO_s':['+',None,None],
                           'COOH_s':[None,'+',None],
                           #'O-CO_s':'initial_state',
                           #'O-O_s':'final_state',
                           }

    # set up calculation
    inputs = {
        'code': catmap_code,
        'energies': energies,
        'rxn_expressions':List(list=[
                        'CO2_g + *_s + + H_g + ele_g <-> ^0.01eV_s <-> COOH_s',
                        'COOH_s + H_g + ele_g <-> CO_s + H2O_g', 
                        'CO_s <-> CO_g + *_s',
                        'H2_g <-> H2_g',
        ]), 
        'surface_names':List(list=['FeNCC', 'FeNNC', 'FeNNNCC', 'FeNNNNCC', 
                'NiNNC', 'NiNNCC', 'NiNNNC', 'NiNNNCC', 'NiNNNNCC', 
                'CoNNC', 'CoNNCC', 'CoNNNNCC', 'CrNCC', 'CrNNC', 'CrNNCC', 
                'CrNNNC', 'CrNNNCC', 'CrNNNNCC', 'MnNCC', 'MnNNCC', 'MnNNNNCC', 
                'RhNCC', 'RhNNC', 'RhNNCC', 'RhNNNCC', 'RhNNNNCC', 'RuNCC', 
                'RuNNC', 'RuNNCC', 'RuNNNCC', 'VNCC', 'VNNC', 'VNNCC', 'TiNCC', 
                'TiNNC', 'TiNNCC', 'TiNNNC', 'TiNNNCC', 'TiNNNNCC', 'TiNNNNCC', 
                'TiNNNNCC', 'TiNNNNCC', 'TiNNNNCC', 
                ]), 
        'descriptor_names':List(list=['CO_s','COOH_s']), 
        'descriptor_ranges':List(list=[[-2, 2], [-2, 2]]), 
        'resolution':Int(20), 
        'voltage':Float(-0.5),
        'temperature':Float(300), 
        'species_definitions':Dict(dict=species_definitions), 
        'gas_thermo_mode':Str('ideal_gas'), 
        'adsorbate_thermo_mode':Str('harmonic_adsorbate'),
        'electrochemical_thermo_mode':Str('simple_electrochemical'),
        'scaling_constraint_dict':Dict(dict=scaling_constraint_dict), 
        'decimal_precision':Int(150), 
        'tolerance':Float(1e-20), 
        'max_rootfinding_iterations':Int(100), 
        'max_bisections':Int(3), 
        #'numerical_solver':Str('coverages'),
        'metadata': {
            'label':'Calculation without CO2 and scaling dictionary provided',
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
@cmdline.utils.decorators.with_dbenv()
@cmdline.params.options.CODE()
def cli(code):
    """Run example.

    Example usage: $ ./example_01.py --code diff@localhost

    Alternative (creates diff@localhost-test code): $ ./example_01.py

    Help: $ ./example_01.py --help
    """
    test_run(code)


if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter
