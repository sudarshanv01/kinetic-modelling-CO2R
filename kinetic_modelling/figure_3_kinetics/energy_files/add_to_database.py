
import sys
import os
import click

@click.command()
@click.option('--file')
@click.option('--potential')
def main(file, potential):

    SinglefileData = DataFactory('singlefile')
    energies_file = os.path.abspath(file)
    energies = SinglefileData(energies_file)
    node = energies.store_all()

    node.set_extra('potential_SHE', float(potential))

    node.label = 'Energy files needed for running CatMAP'
    grp = Group.get(label='energy_files')
    grp.add_nodes(node)

    print(f'Node of pk: {node.pk} created ')



if __name__ == '__main__':
    main()
