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

input_phi_psi = [
	['ASN', 360.00, -56.14],
	['LEU', -43.98, -51.31],
	['TYR', -66.47, -30.90],
	['ILE', -65.22, -45.94],
	['GLN', -64.75, -30.35],
	['TRP', -73.14, -43.42],
	['LEU', -64.88, -43.25],
	['LYS', -59.51, -25.70],
	['ASP', -77.99, -8.82],
	['GLY', 110.78, 8.08],
	['GLY', 55.24, -124.37],
	['PRO', -57.98, -28.77],
	['SER', -81.83, 19.13],
	['SER', -124.06, 13.40],
	['GLY', 67.93, 25.22],
	['ARG', -143.95, 131.30],
	['PRO', -70.10, 160.07],
	['PRO', -69.48, 145.67],
	['PRO', -77.26, 124.22],
	['SER', -78.10, 360.00]
]

class FastaTestCase(unittest.TestCase):

	def setUp(self):
		from bio.fasta import Fasta
		input = 'NLYIQWLKDGGPSSGRPPPS'
		fasta_sequence = Fasta(input)
		self.generated = fasta_sequence.to_pdb()

	def test_generate_molecule(self):
		expected = Molecule(pdb_file_location='tests/generated.pdb')
		rmsd = self.generated.rmsd(expected)
		self.assertTrue(rmsd < 0.0006)

	def test_rotation(self):
		parameters = Parameters()
		self.generated.set_force_field_parameters(parameters)
		topology = self.generated.get_topology()

		phi = []
		psi = []
		for i in range(len(input_phi_psi)):
			phi.append(input_phi_psi[i][1])
			psi.append(input_phi_psi[i][2])
		topology.rotate_backbone('PHI', phi)
		topology.rotate_backbone('PSI', psi)

		expected = Molecule(pdb_file_location='tests/rotated.pdb')
		rmsd = self.generated.rmsd(expected)
		self.assertTrue(rmsd < 0.022)

if __name__ == '__main__':
	unittest.main()
