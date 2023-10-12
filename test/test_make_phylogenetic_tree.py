import os
import pathlib
import tempfile
from unittest import TestCase

from Bio import Phylo

from phylogenetics.make_phylogenetic_tree import generate_tree, graph_tree_object
from test import TEST_DIR, DATA


class TestTreeGen(TestCase):
    def test_generate_tree(self):
        with tempfile.TemporaryDirectory as tmp:
            base_path = pathlib.Path(tmp)
            alignment_path = base_path / 'influenza-a.aln'
            with open(alignment_path, 'w', encoding='utf-8') as aln_writer:
                aln_writer.write('CLUSTAL' + os.linesep)
                aln_writer.write(
                    'Influenza-A_1918      NGCCTGCTTTTGCTATT-------------------------------------------' + os.linesep)
                aln_writer.write(
                    'Influenza-A_H1N1      AGCCTAATTATTGCTGCTAGGAACATAGTGAGAAGAGCTGCAGTGTCAGCAGATCCACTA' + os.linesep)
                aln_writer.write(
                    '                       ****                                                       ' + os.linesep)
            tree = generate_tree(alignment_path, True)
            Phylo.draw_ascii(tree)
            graph_tree_object(tree)

    def test_vars(self):
        self.assertIn('phylo_tree_generator/test', str(TEST_DIR))
        self.assertIn('influenza-a', str(DATA))
