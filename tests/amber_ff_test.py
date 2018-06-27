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
import unittest
import bio
from unittest.mock import patch
from unittest.mock import MagicMock
from bio.parameters99ff import Parameters
from bio.pdb import Molecule


class AmberTestCase(unittest.TestCase):

	def setUp(self):
		self.parameters = Parameters()

	def test_read_amino_acids(self):
		self.assertEqual(len(self.parameters.amino_acids), 94)

	def test_read_angle_type(self):
		self.assertTrue(self.parameters.get_angle_type('CT', 'S', 'S') is not None)

	def test_bond_atom_type(self):
		for amino_acid in self.parameters.amino_acids.values():
			for bond in amino_acid.bonds:
				self.assertTrue(bond.atom_01_type is not None)
				self.assertTrue(bond.atom_02_type is not None)

	def test_read_molecule(self):
		molecule = Molecule(pdb_file_location='tests/1l2y-01.pdb')
		self.assertEqual(len(molecule.amino_acids), 20)


if __name__ == '__main__':
	unittest.main()
