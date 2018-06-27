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
from pyquaternion import Quaternion


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
			self.update_internal_state()
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

	def add_sequential_ids(self):
		for i in range(len(self.atoms)):
			self.atoms[i].serial = i + 1

	def update_internal_state(self):
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

	def reset_starting_point(self):
		for a in self.atoms:
			a.reset_starting_point()

	def move_origin(self, point):
		import numpy
		delta = numpy.array(self.atoms[0].point) - numpy.array(point)
		for a in self.atoms:
			a.translate(-delta)
			a.reset_starting_point()

	def translate(self, shift):
		for atom in self.atoms:
			atom.translate(shift)
		return self

	def translate_starting_point(self, shift):
		for atom in self.atoms:
			atom.translate_starting_point(shift)
		return self

	def rotate_starting_angle(self, angles):
		for atom in self.atoms:
			atom.rotate_starting_angle(angles)
		return self

	def copy(self):
		copied = Molecule()
		for atom in self.atoms:
			copied.atoms.append(atom.copy())
		copied.update_internal_state()
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

		sum = 0
		for i in range(len(m1_atoms)):
			distance = m1_atoms[i].distance_of(m2_atoms[i])
			sum += math.pow(distance, 2)
		return math.sqrt(sum / len(m1_atoms))

	def write_all(self, directory, file_name=None, format=None):
		f = self._open_file(directory, file_name)
		i = 0
		for a in self.atoms:
			i += 1
			atom_name = a.name
			if len(atom_name) <= 3:
				# prepend space to keep compatibility to Tinker pdbtoxyz command
				atom_name = ' ' + atom_name
			line = 'ATOM  {:5} {:4} {:3} A{:4}    {:8.3f}{:8.3f}{:8.3f}  1.00  1.00          {:>2}  '.format(
				i, atom_name, a.res_name, a.res_seq,
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
		# TODO: rename to res_name to match Atom
		self.code = code
		self.ff_code = code
		# TODO: rename to res_seq to match Atom
		self.sequence = sequence
		self.molecule = molecule
		self.atoms_map = {}
		self.force_field_atoms_map = {}
		self.amino_acid_parameters = None

	def set_force_field_parameters(self, parameters):
		if self.code == 'CYS' and 'HG' not in self.atoms_map:
			self.ff_code = 'CYX'

		for a in self.atoms_map.values():
			a.ff_name = a.name

		self.amino_acid_parameters = parameters.get_amino_acid_parameters(self)

		if 'SG' in self.atoms_map and 'HG' not in self.atoms_map:  # disulfid brigde detected
			self.atoms_map['SG'].type = 'S'

		if 'HB1' not in self.atoms_map and 'HB3' in self.atoms_map:
			self.atoms_map['HB2'].ff_name = 'HB1'
			self.atoms_map['HB3'].ff_name = 'HB2'

		if 'HD1' not in self.atoms_map and 'HD3' in self.atoms_map:
			self.atoms_map['HD2'].ff_name = 'HD1'
			self.atoms_map['HD3'].ff_name = 'HD2'

		if 'HE1' not in self.atoms_map and 'HE3' in self.atoms_map:
			self.atoms_map['HE2'].ff_name = 'HE1'
			self.atoms_map['HE3'].ff_name = 'HE2'

		if 'HG1' not in self.atoms_map and 'HG3' in self.atoms_map:
			self.atoms_map['HG2'].ff_name = 'HG1'
			self.atoms_map['HG3'].ff_name = 'HG2'

		# from ILE fix
		if 'HG11' not in self.atoms_map and 'HG13' in self.atoms_map:
			self.atoms_map['HG12'].ff_name = 'HG11'
			self.atoms_map['HG13'].ff_name = 'HG12'

		# TODO: get name mapping from gromacs parameters
		if self.code == 'ILE':
			if 'CD1' in self.atoms_map:
				self.atoms_map['CD1'].ff_name = 'CD'
			if 'HD11' in self.atoms_map:
				self.atoms_map['HD11'].ff_name = 'HD1'
			if 'HD12' in self.atoms_map:
				self.atoms_map['HD12'].ff_name = 'HD2'
			if 'HD13' in self.atoms_map:
				self.atoms_map['HD13'].ff_name = 'HD3'

		if self.is_last():
			self.atoms_map['O'].ff_name = 'OC1'
			self.atoms_map['O'].type = 'O2'
			self.atoms_map['OXT'].ff_name = 'OC2'
			self.atoms_map['OXT'].type = 'O2'

		for a in self.atoms_map.values():
			self.force_field_atoms_map[a.ff_name] = a

		for atom in self.atoms_map.values():
			if atom.ff_name in self.amino_acid_parameters.atoms_map:
				atom.charge = self.amino_acid_parameters.atoms_map[atom.ff_name].charge

	def is_first(self):
		return self.molecule.atoms[0].res_seq == self.sequence

	def is_last(self):
		last_index = len(self.molecule.atoms) - 1
		return self.molecule.atoms[last_index].res_seq == self.sequence

# TODO move to topology
class Bond:
	def __init__(self, atom_01, atom_02, bond_type):
		self.atom_01 = atom_01
		self.atom_02 = atom_02
		self.atoms = (atom_01, atom_02)
		self.bond_type = bond_type


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
			self.ff_name = None

	def move(self, position):
		self.x = position[0]
		self.y = position[1]
		self.z = position[2]
		self.point = [self.x, self.y, self.z]

	def reset_starting_point(self):
		self.original_x = self.x
		self.original_y = self.y
		self.original_z = self.z
		return self

	def translate_starting_point(self, position):
		self.x = self.original_x + position[0]
		self.y = self.original_y + position[1]
		self.z = self.original_z + position[2]
		self.point = [self.x, self.y, self.z]
		return self

	def translate(self, position):
		self.x = self.x + position[0]
		self.y = self.y + position[1]
		self.z = self.z + position[2]
		self.point = [self.x, self.y, self.z]
		# TODO: reset_starting_point ?
		return self

	def rotate_starting_angle(self, angles):
		x_rotation = Quaternion(axis=[1, 0, 0], radians=angles[0])
		y_rotation = Quaternion(axis=[0, 1, 0], radians=angles[1])
		z_rotation = Quaternion(axis=[0, 0, 1], radians=angles[2])
		rotation = x_rotation * y_rotation * z_rotation
		rotated = rotation.rotate([self.original_x, self.original_y, self.original_z])
		self.x = rotated[0]
		self.y = rotated[1]
		self.z = rotated[2]
		self.point = [self.x, self.y, self.z]

	def is_in(self, *args):
		return self in args

	def is_alpha_carbon(self):
		return self.name.startswith('CA')

	def is_backbone(self):
		return self.name in ['N', 'CA', 'C', 'O', 'OXT', 'H', 'H1', 'H2', 'H3', 'HA']

	# http://kinemage.biochem.duke.edu/teaching/anatax/html/anatax.1b.html
	def is_forming_important_dihedral(self):
		return self.name in ['N', 'CA', 'C']

	def is_forming_side_chain_dihedral(self):
		return self.name[0] != 'H' and (not self.is_backbone() or self.name in ['CA'])

	def distance_of(self, a2):
		dx = a2.x - self.x
		dy = a2.y - self.y
		dz = a2.z - self.z
		return math.sqrt(math.pow(dx, 2) + math.pow(dy, 2) + math.pow(dz, 2))

	def copy(self):
		copied = Atom()
		copied.serial = self.serial
		copied.name = self.name
		copied.ff_name = self.ff_name
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

