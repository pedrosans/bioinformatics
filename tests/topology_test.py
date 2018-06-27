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
from bio.topology import Topology
from bio.pdb import Molecule


class TopologyTestCase(unittest.TestCase):

	def setUp(self):
		self.parameters = Parameters()
		self.protein = Molecule(pdb_file_location='tests/data/ref01.pdb')
		self.protein.set_force_field_parameters(self.parameters)
		self.ligand = Molecule(pdb_file_location='tests/data/ligand.pdb')

	def test_read_standard_residue(self):
		topology = Topology(self.protein, self.parameters)
		topology.read_bonds_from_ff()
		topology.mount_topology(read_atom_types_from_bonds=False)
		topology.print_topology()
		self.assertEquals(len(topology.bonds), 415)

	def test_read_non_standard_residue(self):
		topology = Topology(self.ligand, self.parameters)
		topology.read_bonds_from_pdb_file()
		topology.mount_topology(read_atom_types_from_bonds=True)
		topology.print_topology()
		self.assertEquals(len(topology.bonds), 33)


if __name__ == '__main__':
	unittest.main()
