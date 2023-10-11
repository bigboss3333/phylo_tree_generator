import argparse
import pathlib

from Bio import AlignIO
from Bio import Phylo
import matplotlib
import matplotlib.pyplot as plt
from Bio.Cluster import Tree
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

from phylogenetics import ALIGNMENT_FORMAT


def generate_tree(alignment_path: pathlib.Path, rooted=True) -> Tree:
    """ This will take an alignment file and convert it to a tree object.
    The tree will represent the relationship between the clades, and provides
    data needed to map a taxonomic relationship.

    :param alignment_path:
    :param rooted:
    :return: Tree object
    """
    with open(alignment_path, encoding='utf-8') as alignment_file_stream:
        alignment = AlignIO.read(alignment_file_stream, ALIGNMENT_FORMAT)
    calculator = DistanceCalculator('identity')
    constructor = DistanceTreeConstructor(calculator)
    tree_object = constructor.build_tree(alignment)
    tree_object.rooted = rooted
    tree_object.id = alignment_path.name.replace('.aln', '')
    Phylo.write(tree_object, alignment_path.name.replace('.aln', '.xml'), "phyloxml")

    return tree_object


def graph_tree(tree_object: Tree):
    """ Generates a phylogenetic tree plot from the tree object

    :param tree_object: Object created from alignment file mapping the distance between species
    :return: Nothing
    """
    fig = plt.figure(figsize=(20, 5))
    matplotlib.rc('font', size=14)
    matplotlib.rc('ytick', labelsize=12)
    matplotlib.rc('xtick', labelsize=12)

    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree_object, axes=axes)


def parse_args() -> argparse.Namespace:
    """ Parse user args

    :return: namespace of arguments
    """
    parser = argparse.ArgumentParser(description='Generate Phylogenetic Tree')
    parser.add_argument('alignment', type=pathlib.Path,
                        help='an alignment file used to generate phylogenetic tree, must be in clustal format')
    parser.add_argument('-r', '--rooted',
                        action='store_true')

    return parser.parse_args()


def main():
    """ Entry point for script

    :return: Nothing
    """
    args = parse_args()
    tree = generate_tree(args.alignment, args.rooted)
    graph_tree(tree)


if __name__ == '__main__':
    main()
