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
from inf.optimization import GradientDescent
from bio.molecular_dynamics import ForceField


class MinimizationTestCase(unittest.TestCase):

	def setUp(self):
		self.molecule = Molecule(pdb_file_location='tests/data/ref01.pdb')
		self.parameters = Parameters()
		self.molecule.set_force_field_parameters(self.parameters)
		self.optimization = GradientDescent(self.atom_energy_function)
		self.optimization.dimensions_number = 3
		self.optimization.iterations = 40
		self.optimization.precision = 0.01
		self.optimization.small_step = 0.00015
		self.optimization.smallest_step = self.optimization.small_step / 128
		self.force_field = ForceField(self.parameters)
		self.atom_index = 0

	def test_gradient_descent(self):
		self.atom_index = 0
		for j in range(5):
			self.optimization.run()
			atom = self.molecule.atoms[self.atom_index]
			atom.translate_starting_point(self.optimization.best_particle.dimensions)
			self.atom_index += 1

		energy = self.force_field.calculate_energy(self.molecule.get_topology())
		# self.force_field.print_energy()
		self.assertEqual(energy, 599.882043366674)

	def atom_energy_function(self, particle):
		a = self.molecule.atoms[self.atom_index]
		a.translate_starting_point(particle.dimensions)
		return self.force_field.calculate_energy(self.molecule.get_topology(), atom=a, protect_bonds=False)


if __name__ == '__main__':
	unittest.main()
