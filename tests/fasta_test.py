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
from bio.fasta import Fasta

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
import numpy as np
import quaternion as quat
import inf.geometry


class FastaTestCase(unittest.TestCase):

	def setUp(self):
		import warnings
		warnings.filterwarnings("ignore", message="numpy.dtype size changed")
		warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
		self.input = 'NLYIQWLKDGGPSSGRPPPS'
		fasta_sequence = Fasta(self.input)
		self.generated = fasta_sequence.to_pdb()
		self.compute_matrix = 1
		self.compute_point = 10000000

	def rotate_quaternion(self, axis, v, theta):
		vector = np.array([0.] + v)
		rot_axis = np.array([0.] + axis)
		axis_angle = (theta * 0.5) * rot_axis / np.linalg.norm(rot_axis)

		vec = quat.quaternion(*v)
		qlog = quat.quaternion(*axis_angle)
		q = np.exp(qlog)

		for i in range(self.compute_point):
			result = q * vec * np.conjugate(q)
		return result

	def rotate_raw(self, axis, v, theta):
		M = inf.geometry.rotation_matrix(axis, theta)
		point = np.array(v)
		for i in range(self.compute_point):
			result = np.dot(M, point)
		return result

	def ignore_test_stress_01(self):
		for i in range(self.compute_matrix):
			self.rotate_raw([4, 4, 1], [3, 5, 0], 1.2)
		v = self.rotate_raw([4, 4, 1], [3, 5, 0], 1.2)
		self.assertEquals(v[0], 2.7491163796476545)
		self.assertEquals(v[1], 4.771809323355129)
		self.assertEquals(v[2], 1.9162971879888682)

	def ignore_test_stress_02(self):
		for i in range(self.compute_matrix):
			self.rotate_quaternion([4, 4, 1], [3, 5, 0], 1.2)
		v = self.rotate_quaternion([4, 4, 1], [3, 5, 0], 1.2)
		self.assertEquals(v.x, 2.7491163796476545)
		self.assertEquals(v.y, 4.771809323355129)
		self.assertEquals(v.z, 1.9162971879888677)

	def ignore_test_stress_03(self):
		for i in range(30):
			molecule = Fasta(self.input).to_pdb()
		self.assertEqual(len(molecule.atoms), 304)

	def test_generate_molecule(self):
		expected = Molecule(pdb_file_location='tests/data/generated.pdb')
		rmsd = self.generated.rmsd(expected)
		self.assertLess(rmsd, 0.0006)

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

		expected = Molecule(pdb_file_location='tests/data/rotated.pdb')
		rmsd = self.generated.rmsd(expected)
		self.assertLess(rmsd, 0.05)

if __name__ == '__main__':
	unittest.main()
