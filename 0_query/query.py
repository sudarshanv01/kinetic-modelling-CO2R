

from aiida.orm import QueryBuilder, Node
from useful_classes import bcolors
from pprint import pprint

qb = QueryBuilder()
qb.append(Group, filters={'id':{'==':22}}, tag='group')
qb.append(CalcJobNode, with_group='group')
results = qb.all()

for index, row in enumerate(results):
    node = row[0]
    print(f'{bcolors.OKBLUE} Calculation {index}: "{node.label}" {bcolors.ENDC}')
    print(f'\t Node: {node.pk}')
    print(f'\t Descriptors: {node.inputs.descriptor_names.get_list()}')
    print(f'\t Voltage: {node.inputs.voltage.value} V vs. {node.inputs.potential_reference_scale.value}')
    print(f'\t Reaction Expressions')
    pprint(f'{node.inputs.rxn_expressions.get_list()}')

