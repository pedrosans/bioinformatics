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

import os, math, re
import inf.geometry
import numpy as np


class Molecule:

	def __init__(self, pdb_file_content=None, pdb_file_location=None):
		self.amino_acids = []
		self.atoms = []
		self.backbone = []
		self.alphas = []
		self.backbone = []
		self.pdb_file_location = pdb_file_location
		self.force_field_parameters = None
		if pdb_file_location:
			f = open(pdb_file_location, 'r')
			pdb_file_content = f.read()
			f.close()
		self.pdb_file_content = pdb_file_content
		if pdb_file_content:
			for line in pdb_file_content.splitlines():
				if 'ATOM' in line:
					pdb_atom = Atom(line)
					self.atoms.append(pdb_atom)
			self._populate_subsets()
		self.topology = None

	def set_force_field_parameters(self, ffp):
		self.force_field_parameters = ffp
		for amino_acid in self.amino_acids:
			amino_acid.set_force_field_parameters(ffp)

	def get_topology(self):
		from bio.topology import Topology
		if not self.topology:
			self.topology = Topology(self, self.force_field_parameters)
			self.topology.mount_topology()
		return self.topology

	def _populate_subsets(self):
		last_amino_acid = None
		for a in self.atoms:
			if not last_amino_acid or last_amino_acid.sequence != a.res_seq:
				last_amino_acid = AminoAcid(a.res_name, a.res_seq, self)
				self.amino_acids.append(last_amino_acid)
			last_amino_acid.atoms_map[a.name] = a
			if a.is_alpha_carbon():
				self.alphas.append(a)
			if a.is_backbone():
				self.backbone.append(a)

	def shift_coordinate_origing(self):
		for a in self.atoms:
			a.shift_coordinate_origing()

	def dislocate(self, shift):
		for atom in self.atoms:
			atom.translate(shift)
		return self

	def copy(self):
		copied = Molecule()
		for atom in self.atoms:
			copied.atoms.append(atom.copy())
		copied._populate_subsets()
		return copied

	def rmsd(self, m2, comparison=None):
		if comparison is 'ALPHAS':
			m1_atoms = self.alphas
			m2_atoms = m2.alphas
		elif comparison is 'BACKBONE':
			m1_atoms = self.backbone
			m2_atoms = m2.backbone
		else:
			m1_atoms = self.atoms
			m2_atoms = m2.atoms

		delta_sum = 0
		for i in range(len(m1_atoms)):
			distance = m1_atoms[i].distance_of(m2_atoms[i])
			delta_sum += math.pow(distance, 2)
		return math.sqrt(delta_sum / len(self.atoms))

	def write_all(self, directory, file_name=None, format=None):
		f = self._open_file(directory, file_name)
		i = 0
		for a in self.atoms:
			i += 1
			line = 'ATOM  {:5} {:4} {:3} A{:4}    {:8.3f}{:8.3f}{:8.3f}  1.00  1.00          {:>2}  '.format(
				i, a.name, a.res_name, a.res_seq,
				a.x, a.y, a.z,
				a.symbol
			)
			f.write(line + '\n')
		f.close()

	def write(self, directory, file_name=None, format=None):
		if not file_name:
			segs = self.pdb_file_location.split('/')
			file_name = 'generated-' + segs[len(segs) - 1]
		f = self._open_file(directory, file_name)
		i = 0
		for line in self.pdb_file_content.splitlines():
			if 'ATOM' in line:
				atom = self.atoms[i]
				i += 1
				new_line = '{}{}{}\n'.format(line[0:30], '{:8.3f}{:8.3f}{:8.3f}'.format(atom.x, atom.y, atom.z), line[54:])
				if format == 'gmx':
					new_line = new_line.replace(atom.original_name, atom.name)
				f.write(new_line)
			else:
				f.write(line + '\n')
		f.close()

	def _open_file(self, directory, file_name):
		if not os.path.exists(directory):
			os.makedirs(directory)
		result_pdb_file = directory + '/' + file_name
		if os.path.exists(result_pdb_file):
			os.remove(result_pdb_file)
		return open(result_pdb_file, 'w')


class AminoAcid:

	def __init__(self, code, sequence, molecule):
		self.code = code
		self.ff_code = code
		self.sequence = sequence
		self.molecule = molecule
		self.atoms_map = {}
		self.force_field_atoms_map = {}
		self.amino_acid_parameters = None

	def set_force_field_parameters(self, parameters):
		if self.code == 'CYS' and 'HG' not in self.atoms_map:
			self.ff_code = 'CYX'

		self.amino_acid_parameters = parameters.get_amino_acid_parameters(self)

		if 'SG' in self.atoms_map and 'HG' not in self.atoms_map:  # disulfid brigde detected
			self.atoms_map['SG'].type = 'S'

		for a in self.atoms_map.values():
			a.set_force_field_parameters(parameters)
			if self.is_last():
				if a.name == 'O':
					a.ff_name = 'OC1'
					a.type = 'O2'
				if a.name == 'OXT':
					a.ff_name = 'OC2'
					a.type = 'O2'
			self.force_field_atoms_map[a.ff_name] = a

		for atom in self.atoms_map.values():
			if atom.ff_name in self.amino_acid_parameters.atoms_map:
				atom.charge = self.amino_acid_parameters.atoms_map[atom.ff_name].charge

	def is_first(self):
		return self.molecule.atoms[0].res_seq == self.sequence

	def is_last(self):
		last_index = len(self.molecule.atoms) - 1
		return self.molecule.atoms[last_index].res_seq == self.sequence


class Bond:
	def __init__(self, atom_01, atom_02, bond_type):
		self.atom_01 = atom_01
		self.atom_02 = atom_02
		self.atoms = (atom_01, atom_02)
		self.bond_type = bond_type


class Angle:

	cache = {}

	def __init__(self, atom_01, atom_02, atom_03, angle_type):
		self.atom_01 = atom_01
		self.atom_02 = atom_02
		self.atom_03 = atom_03
		self.atoms = (atom_01, atom_02, atom_03)
		self.angle_type = angle_type

	def get_angle(self):
		ba = np.array(self.atom_01.point) - np.array(self.atom_02.point)
		bc = np.array(self.atom_03.point) - np.array(self.atom_02.point)
		cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
		angle_radians = np.arccos(cosine_angle)
		return np.degrees(angle_radians)


class Dihedral:

	def __init__(self, atom_01, atom_02, atom_03, atom_04, dihedral_types):
		self.atom_01 = atom_01
		self.atom_02 = atom_02
		self.atom_03 = atom_03
		self.atom_04 = atom_04
		self.atoms = (atom_01, atom_02, atom_03, atom_04)
		self.dihedral_types = dihedral_types

	def get_torsion_angle(self):
		return inf.geometry.calculate_torsion(
			self.atom_01.point, self.atom_02.point, self.atom_03.point, self.atom_04.point)

	def get_out_of_plane_angle(self):
		return inf.geometry.calculate_out_of_plane_angle(
			self.atom_01.point, self.atom_02.point, self.atom_03.point, self.atom_04.point)


class Atom:

	# https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html
	def __init__(self, pdb_atom_line=None):
		if pdb_atom_line:
			# serial number
			self.serial = int(pdb_atom_line[6:11])
			self.original_name = pdb_atom_line[12:16].strip()
			self.name = self.original_name
			self.type = None
			self.charge = None
			# altLoc - Alternate location indicator.
			self.alt_loc = pdb_atom_line[16:17]
			# resName - Residue name
			self.res_name = pdb_atom_line[17:20]
			# resSeq Residue sequence number.
			self.res_seq = int(pdb_atom_line[22:26])
			self.x = self.original_x = float(pdb_atom_line[30:38])
			self.y = self.original_y = float(pdb_atom_line[38:46])
			self.z = self.original_z = float(pdb_atom_line[46:54])
			self.symbol = pdb_atom_line[76:77].strip()
			if not self.symbol:
				self.symbol = self.name[0:1]
			self.point = [self.x, self.y, self.z]

	def set_force_field_parameters(self, ffp):
		# self.name = ffp.get_force_field_atom_name(self.name, self.res_name)
		self.ff_name = self.name
		if self.res_name != 'ALA':
			self.ff_name = re.sub('^HB2$', 'HB1', self.ff_name)  # corrige SER/2
			self.ff_name = re.sub('^HB3$', 'HB2', self.ff_name)  # corrige SER/2
		self.ff_name = re.sub('^HG2$', 'HG1', self.ff_name)  # corrige GLU/4
		self.ff_name = re.sub('^HG3$', 'HG2', self.ff_name)  # corrige GLU/4
		if self.res_name != 'HIS':
			self.ff_name = re.sub('^HD2$', 'HD1', self.ff_name)  # corrige PRO/7
			self.ff_name = re.sub('^HD3$', 'HD2', self.ff_name)  # corrige PRO/7
		self.ff_name = re.sub('^HE2$', 'HE1', self.ff_name)  # corrige LYS/14
		self.ff_name = re.sub('^HE3$', 'HE2', self.ff_name)  # corrige LYS/14

		if self.res_name == 'ILE':
			self.ff_name = re.sub('^CD1$', 'CD', self.ff_name)  # corrige ILE/29
			self.ff_name = re.sub('^HD11$', 'HD1', self.ff_name)  # corrige ILE/29
			self.ff_name = re.sub('^HD12$', 'HD2', self.ff_name)  # corrige ILE/29
			self.ff_name = re.sub('^HD13$', 'HD3', self.ff_name)  # corrige ILE/29
			self.ff_name = re.sub('^HG12$', 'HG11', self.ff_name)  # corrige ILE/29
			self.ff_name = re.sub('^HG13$', 'HG12', self.ff_name)  # corrige ILE/29

	def shift_coordinate_origing(self):
		self.original_x = self.x
		self.original_y = self.y
		self.original_z = self.z

	def translate(self, position):
		self.x = self.original_x + position[0]
		self.y = self.original_y + position[1]
		self.z = self.original_z + position[2]
		self.point = [self.x, self.y, self.z]

	def is_in(self, *args):
		return self in args

	def is_alpha_carbon(self):
		return self.name.startswith('CA')

	def is_backbone(self):
		return self.name in ['N', 'CA', 'C']

	def distance_of(self, a2):
		dx = a2.x - self.x
		dy = a2.y - self.y
		dz = a2.z - self.z
		return math.sqrt(math.pow(dx, 2) + math.pow(dy, 2) + math.pow(dz, 2))

	def copy(self):
		copied = Atom()
		copied.serial = self.serial
		copied.name = self.name
		copied.type = self.type
		copied.charge = self.charge
		copied.alt_loc = self.alt_loc
		copied.res_name = self.res_name
		copied.res_seq = self.res_seq
		copied.x = copied.original_x = self.x
		copied.y = copied.original_y = self.y
		copied.z = copied.original_z = self.z
		copied.symbol = self.symbol
		copied.point = (self.x, self.y, self.z)
		return copied

