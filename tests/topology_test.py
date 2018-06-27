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
		self.ligand = Molecule(pdb_file_location='tests/data/ligand.pdb')
		self.topology = Topology(self.ligand, self.parameters)

	def test_read_molecule(self):
		self.topology.read_bonds_from_pdb_file()
		self.topology.mount_topology(read_atom_types_from_bonds=True)
		self.topology.print_topology()
		self.assertEquals(len(self.topology.bonds), 33)


if __name__ == '__main__':
	unittest.main()
