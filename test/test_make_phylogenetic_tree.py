import os
import pathlib
import tempfile
from unittest import TestCase

from Bio import Phylo

from phylogenetics.make_phylogenetic_tree import generate_tree


class TestTreeGen(TestCase):
    def test_generate_tree(self):
        #with tempfile.TemporaryDirectory as tmp:
        tmp = tempfile.mkdtemp()
        base_path = pathlib.Path(tmp)
        alignment_path = base_path / 'influenza-a.aln'
        with open(alignment_path, 'w', encoding='utf-8') as aln_writer:
            aln_writer.write('CLUSTAL'+os.linesep)
            aln_writer.write('Influenza-A_1918      NGCCTGCTTTTGCTATT-------------------------------------------'+os.linesep)
            aln_writer.write('Influenza-A_H1N1      AGCCTAATTATTGCTGCTAGGAACATAGTGAGAAGAGCTGCAGTGTCAGCAGATCCACTA'+os.linesep)
            aln_writer.write('                       ****                                                       '+os.linesep)

        tree = generate_tree(alignment_path, True)
        Phylo.draw_ascii(tree)

