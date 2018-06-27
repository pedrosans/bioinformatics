"""
Copyright 2018 Pedro Santos <pedrosans@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from bio.parameters99ff import Parameters
from bio.pdb import Molecule
from bio.topology import Dihedral
import numpy as np

code_map = {
	'A': 'ALA',  # Alanine
	'B': 'ASX',  # Aspartic acid (D) or Asparagine (N)
	'C': 'CYS',  # Cysteine
	'D': 'ASP',  # Aspartic acid
	'E': 'GLU',  # Glutamic acid
	'F': 'PHE',  # Phenylalanine
	'G': 'GLY',  # Glycine
	'H': 'HIS',  # Histidine
	'I': 'ILE',  # Isoleucine
	'K': 'LYS',  # Lysine
	'L': 'LEU',  # Leucine
	'M': 'MET',  # Methionine
	'N': 'ASN',  # Asparagine
	'P': 'PRO',  # Proline
	'Q': 'GLN',  # Glutamine
	'R': 'ARG',  # Arginine
	'S': 'SER',  # Serine
	'T': 'THR',  # Threonine
	'V': 'VAL',  # Valine
	'W': 'TRP',  # Tryptophan
	'Y': 'TYR',  # Tyrosine
	'Z': 'GLX',  # Glutamic acid (E) or Glutamine (Q)
}


class Fasta:

	def __init__(self, sequence):
		self.sequence = sequence
		self.molecule = Molecule()

	def to_pdb(self):
		self._populate_molecule()
		self.molecule.add_sequential_ids()
		self.molecule.update_internal_state()
		parameters = Parameters()
		self.molecule.set_force_field_parameters(parameters)
		self._rotate_to_default_position()
		return self.molecule

	def _populate_molecule(self):
		pointer = np.array([0, 0, 0])
		c_n_separator = np.array([1.3350, 0, 0])
		i = 0
		for a in self.sequence:
			i += 1
			abbreviation = code_map[a]
			residue = Molecule(pdb_file_location='data/catalog/' + abbreviation + '.pdb')
			if i == 1:
				Fasta._add_n_terminal_hydrogen(residue)

			if i == len(self.sequence):
				Fasta._add_c_terminal_oxygen(residue)

			n_terminal = residue.amino_acids[0].atoms_map['N']
			c_terminal = residue.amino_acids[0].atoms_map['C']
			delta = np.array(n_terminal.point) - np.array(pointer)
			residue.translate(-delta)
			for a in residue.atoms:
				a.res_seq = i
				self.molecule.atoms.append(a.copy())
			pointer = np.array(c_terminal.point) + c_n_separator

	def _rotate_to_default_position(self):
		topology = self.molecule.get_topology()
		for i in range(len(self.molecule.amino_acids)):
			amino_acid = self.molecule.amino_acids[i]
			backbone_torsions_map = topology.backbone_torsions_map[amino_acid.sequence]

			if not amino_acid.is_last():
				topology.fix_nitrogen(amino_acid, self.molecule.amino_acids[i + 1])
				topology.fix_alpha_carbon(amino_acid, self.molecule.amino_acids[i + 1])

			if 'PSI' in backbone_torsions_map:
				topology.rotate_proper_dihedral(backbone_torsions_map['PSI'], 180.0)
			if 'PHI' in backbone_torsions_map:
				topology.rotate_proper_dihedral(backbone_torsions_map['PHI'], 180.0)
			if 'OMEGA' in backbone_torsions_map:
				topology.rotate_proper_dihedral(backbone_torsions_map['OMEGA'], 180.0)

			# if not amino_acid.is_first(): topology.set_side_chain_torsion(amino_acid, 0)

		for i in range(len(self.molecule.amino_acids)):
			amino_acid = self.molecule.amino_acids[i]

			if not amino_acid.is_last():
				topology.fix_oxygen(amino_acid, self.molecule.amino_acids[i + 1])
				topology.fix_amine_hydrogen_torsion(amino_acid, self.molecule.amino_acids[i + 1])
				topology.fix_amine_hydrogen(amino_acid, self.molecule.amino_acids[i + 1])
			topology.fix_amine_hydrogen_angle_to_ca(amino_acid)

	@staticmethod
	def _add_n_terminal_hydrogen(residue):
		n = residue.amino_acids[0].atoms_map['N']
		h1 = residue.amino_acids[0].atoms_map['H']
		h1.name = 'H1'
# ATOM      9  H1  ASN A   1      -8.330   3.957   0.261  1.00  0.00           H
# ATOM     10  H2  ASN A   1      -8.740   5.068  -0.889  1.00  0.00           H
# ATOM     11  H3  ASN A   1      -9.877   4.041  -0.293  1.00  0.00           H
		h2 = h1.copy().translate([-0.3, -0.5, -1.3])
		h3 = h1.copy().translate([0, -1.1, -0.3])
		h2.name = 'H2'
		h3.name = 'H3'
		residue.atoms.append(h2)
		residue.atoms.append(h3)
		residue.update_internal_state()

	@staticmethod
	def _add_c_terminal_oxygen(residue):
		o_atom = residue.amino_acids[0].atoms_map['O']
		oxt_atom = o_atom.copy().translate([-1.5,  +0.0, 1.7])
		oxt_atom.name = 'OXT'
		residue.atoms.append(oxt_atom)
		residue.update_internal_state()
