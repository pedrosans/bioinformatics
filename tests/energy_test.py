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
from bio.parameters99ff import Parameters
from bio.pdb import Molecule
from bio.molecular_dynamics import ForceField


class EnergyTestCase(unittest.TestCase):

	def setUp(self):
		self.parameters = Parameters()
		self.ligand = Molecule(pdb_file_location='tests/data/ligand.pdb')
		self.ligand.set_force_field_parameters(self.parameters)
		self.protein = Molecule(pdb_file_location='tests/data/ref01.pdb')
		self.protein.set_force_field_parameters(self.parameters)
		self.force_field = ForceField(self.parameters)

	def test_non_standard_res_energy(self):
		self.force_field.calculate_energy(self.ligand, test_electrostatic_only=True)
		self.force_field.print_energy()
		self.assertEquals(self.force_field.energy, 553.28142232581)

	def test_standard_res_energy(self):
		self.force_field.calculate_energy(self.protein, test_electrostatic_only=False)
		self.force_field.print_energy()
		self.assertEquals(self.force_field.energy, 606.3322671759454)


if __name__ == '__main__':
	unittest.main()
