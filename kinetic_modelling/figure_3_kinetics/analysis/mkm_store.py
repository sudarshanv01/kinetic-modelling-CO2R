
""" For a given CalcJobNode store the production map """

import sys
import os
import json
from pathlib import Path
from aiida.plugins import CalculationFactory

GROUPNAME = "kinetic_models/descriptors_CO2_COOH"
TYPE_OF_CALC = CalculationFactory('catmap') 

def main():
    """Get all the quantities from the calculation and put them into a json."""

    qb = QueryBuilder()
    qb.append(Group, filters={'label':GROUPNAME}, tag='Group')
    qb.append(TYPE_OF_CALC, with_group='Group', tag='calctype')

    data_tot = {}
    for node in qb.all(flat=True):
        print(node)

        descriptors = node.inputs.descriptor_names.get_list()
        species_definitions = node.inputs.species_definitions
        facet = species_definitions['s']['site_names']
        surfaces = node.inputs.surface_names.get_list()
        potential = node.inputs.voltage.value
        pH = node.inputs.pH.value
        coverage_map = node.outputs.coverage_map.get_list()
        production_rate = node.outputs.production_rate_map.get_list()
        energy_file = node.inputs.energies.get_content()
        pk = node.pk

        data = {}
        data['facet'] = facet
        data['surfaces'] = surfaces
        data['descriptors'] = descriptors
        data['potential'] = potential
        data['pH'] = pH
        data['coverage_map'] = coverage_map
        data['production_rate'] = production_rate 
        data['energy_file'] = energy_file
        data_tot[pk] = data

    with open('aiida_output/kinetic_model_data.json', 'w') as handle:
        json.dump(data_tot, handle, indent=4)


if __name__ == '__main__':
    Path('aiida_out').mkdir(parents=True, exist_ok=True)
    main()

