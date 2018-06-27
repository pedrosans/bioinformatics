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

import re


class Parameters:

	def __init__(self):
		self.atom_types = {}
		self.force_field_atom_names = {}
		self.amino_acids = {}
		self.bond_types = {}
		self.angle_types = {}
		self.proper_dihedral_types = {}
		self.improper_dihedral_types = {}
		self.load_bond_types()
		self.load_amino_acid_parameters()
		self.load_atom_types()
		self._load_ff_atom_name()

	def get_amino_acid_parameters(self, amino_acid):
		amino_acid_code = amino_acid.ff_code
		if amino_acid_code == 'HIS':
			amino_acid_code = 'HID'
		if amino_acid.is_first():
			amino_acid_code = 'N' + amino_acid_code
		if amino_acid.is_last():
			amino_acid_code = 'C' + amino_acid_code
		if amino_acid_code in self.amino_acids:
			return self.amino_acids[amino_acid_code]
		else:
			return None

	def get_atom_type(self, type):
		return self.atom_types[type]

	def get_bond_type(self, atom_01_code, atom_02_code):
		if atom_01_code in self.bond_types and atom_02_code in self.bond_types[atom_01_code]:
			return self.bond_types[atom_01_code][atom_02_code]
		elif atom_02_code in self.bond_types and atom_01_code in self.bond_types[atom_02_code]:
			return self.bond_types[atom_02_code][atom_01_code]
		else:
			return None

	def get_angle_type(self, type_01, type_02, type_03):
		angle_type = self._get_matching_angle_type(type_01, type_02, type_03)
		if not angle_type:
			angle_type = self._get_matching_angle_type(type_03, type_02, type_01)
		return angle_type

	def _get_matching_angle_type(self, type_01, type_02, type_03):
		matches_01 = self._get_map_match(type_01, self.angle_types)
		for match_01 in matches_01:
			matches_02 = self._get_map_match(type_02, match_01)
			for match_02 in matches_02:
				matches_03 = self._get_map_match(type_03, match_02)
				if len(matches_03) > 0:
					return matches_03[0]

	def get_proper_dihedral_type(self, type_01, type_02, type_03, type_04):
		result = []
		dihedral_type = self._get_dihedral_type(type_01, type_02, type_03, type_04, self.proper_dihedral_types)
		if dihedral_type:
			result = list(dihedral_type)
		mirror_components = self._get_dihedral_type(type_04, type_03, type_02, type_01, self.proper_dihedral_types)
		if mirror_components:
			result += list(mirror_components)
		return result

	def get_improper_dihedral_type(self, type_01, type_02, type_03, type_04):
		return self._get_dihedral_type(type_01, type_02, type_03, type_04, self.improper_dihedral_types)

	def _get_dihedral_type(self, type_01, type_02, type_03, type_04, map):
		matches_01 = self._get_map_match(type_01, map)
		for match_01 in matches_01:
			matches_02 = self._get_map_match(type_02, match_01)
			for match_02 in matches_02:
				matches_03 = self._get_map_match(type_03, match_02)
				for match_03 in matches_03:
					matches_04 = self._get_map_match(type_04, match_03)
					if len(matches_04) > 0:
						return matches_04[0]
		return None

	@staticmethod
	def _get_map_match(type, map):
		matches = []
		if type in map:
			matches.append(map[type])
		if type.startswith('C') and 'C*' in map:
			matches.append(map['C*'])
		if type.startswith('N') and 'N*' in map:
			matches.append(map['N*'])
		if 'X' in map:
			matches.append(map['X'])
		return matches

	def load_amino_acid_parameters(self):
		f = open('data/amber99.ff/aminoacids.rtp', 'r')
		content = f.read()
		f.close()
		last_aminoacid = None
		last_property = None
		for line in content.splitlines():
			if not line or line.startswith(';'):
				continue
			if line.startswith('[ '):
				a = AminoAcidParameters(line[1:6].strip())
				self.amino_acids[a.code] = a
				last_aminoacid = a
			elif line.startswith(' ['):
				last_property = line[2:8].strip()
			else:
				if last_property == 'atoms':
					atom = AminoAcidAtomParameters(line)
					last_aminoacid.atoms_map[atom.code] = atom
				elif last_property == 'bonds':
					bond = AminoAcidBondParameters(line, last_aminoacid)
					if not bond.is_marker():
						last_aminoacid.bonds.append(bond)

	def load_bond_types(self):
		f = open('data/amber99.ff/ffbonded.itp', 'r')
		content = f.read()
		f.close()
		last_type = None
		for line in content.splitlines():
			if line.startswith('['):
				last_type = line
				continue
			if line.startswith(';') or not line:
				continue
			if 'bondtypes' in last_type:
				b = BondType(line)
				if b.atom_01_type not in self.bond_types:
					self.bond_types[b.atom_01_type] = {}
				self.bond_types[b.atom_01_type][b.atom_02_type] = b
			elif 'angletypes' in last_type:
				a = AngleType(line)
				if not a.atom_01_type in self.angle_types:
					self.angle_types[a.atom_01_type] = {}
				if not a.atom_02_type in self.angle_types[a.atom_01_type]:
					self.angle_types[a.atom_01_type][a.atom_02_type] = {}
				self.angle_types[a.atom_01_type][a.atom_02_type][a.atom_03_type] = a
			elif 'dihedraltypes' in last_type:
				d = DihedralType(line)
				dihedral_map = None
				if d.function_type == 4:
					dihedral_map = self.improper_dihedral_types
				elif d.function_type == 9:
					dihedral_map = self.proper_dihedral_types
				if d.atom_01_type not in dihedral_map:
					dihedral_map[d.atom_01_type] = {}
				if d.atom_02_type not in dihedral_map[d.atom_01_type]:
					dihedral_map[d.atom_01_type][d.atom_02_type] = {}
				if d.atom_03_type not in dihedral_map[d.atom_01_type][d.atom_02_type]:
					dihedral_map[d.atom_01_type][d.atom_02_type][d.atom_03_type] = {}
				if d.atom_04_type not in dihedral_map[d.atom_01_type][d.atom_02_type][d.atom_03_type]:
					dihedral_map[d.atom_01_type][d.atom_02_type][d.atom_03_type][d.atom_04_type] = []
				dihedral_map[d.atom_01_type][d.atom_02_type][d.atom_03_type][d.atom_04_type].append(d)


	def load_atom_types(self):
		f = open('data/amber99.ff/ffnonbonded.itp', 'r')
		content = f.read()
		f.close()
		last_type = None
		for line in content.splitlines():
			if not line.startswith('[') and not line.startswith(';'):
				a = AtomType(line)
				self.atom_types[a.type] = a

	def _load_ff_atom_name(self):
		f = open('data/amber99.ff/xlateat.dat', 'r')
		content = f.read()
		f.close()
		for line in content.splitlines():
			if not line.startswith('26'):
				a = ForceFieldAtomParameter(line)
				self.force_field_atom_names[a.pdb_name] = a

	def get_force_field_atom_name(self, atom_name, residue_name):
		if residue_name == 'ILE':
			print('asdf')
		if atom_name in self.force_field_atom_names:
			p = self.force_field_atom_names[atom_name]
			if p.scope in (residue_name, 'protein'):
				return self.force_field_atom_names[atom_name].ff_name
		return atom_name


class ForceFieldAtomParameter:

	def __init__(self, line):
		values = line.split()
		self.scope = values[0]
		self.pdb_name = values[2]
		self.ff_name = values[1]


class AtomType:

	def __init__(self, line):
		values = line.split()
		self.type = values[0]
		self.number = values[1]
		# m a.m.u.
		self.mass = float(values[2])
		# q electron
		self.charge = float(values[3])
		self.particle_type = values[4]
		# VDW radius nm
		self.sigma = float(values[5])
		# VDW attraction magnitude KJ/mol
		self.epsilon = float(values[6])

class DihedralType():

	def __init__(self, line):
		values = line.split()
		self.atom_01_type = values[0]
		self.atom_02_type = values[1]
		self.atom_03_type = values[2]
		self.atom_04_type = values[3]
		self.function_type = int(values[4])
		#degrees
		self.barrier_minimum_offset_angle = float(values[5])
		#kJ mol −1
		self.rotation_barrier_height = float(values[6])
		#as pn in ff99
		self.frequency_of_barrier = int(values[7])

class BondType():

	def __init__(self, line):
		self.atom_01_type = line[0:4].strip()
		self.atom_02_type = line[5:7].strip()
		self.equilibrium_bond_length = float(line[21:28].strip())
		#kJ mol −1 nm −2 from gromac manual
		self.bond_sprint_constant = float(line[31:39].strip())


class AngleType():

	def __init__(self, line):
		self.atom_01_type = line[0:4].strip()
		self.atom_02_type = line[4:8].strip()
		self.atom_03_type = line[8:12].strip()
		self.equilibrium_bond_angle = float(line[25:32].strip())
		#kJ mol −1 rad −2
		#k θ (kJ mol −1 rad −2 ) - from gromacs manual
		self.angle_sprint_constant = float(line[35:43].strip())


class AminoAcidParameters():
	def __init__(self, code):
		self.code = code
		self.atoms_map = {}
		self.bonds = []


class AminoAcidAtomParameters():
	def __init__(self, line):
		values = line.split()
		self.code = values[0]
		self.type = values[1]
		#symbol q unity electron
		self.charge = float(values[2])
		self.charge_group = int(values[3])


class AminoAcidBondParameters:

	def __init__(self, line, amino_acid_parameters):
		self.atom_01_parameter = line[0:6].strip()
		self.atom_02_parameter = line[6:12].strip()
		self.atom_01_code = self.atom_01_parameter.replace('-', '').replace('+', '')
		self.atom_02_code = self.atom_02_parameter.replace('-', '').replace('+', '')
		self.atom_01_type = self.atom_02_type = None
		if self.atom_01_code in amino_acid_parameters.atoms_map:
			self.atom_01_type = amino_acid_parameters.atoms_map[self.atom_01_code].type
		if self.atom_02_code in amino_acid_parameters.atoms_map:
			self.atom_02_type = amino_acid_parameters.atoms_map[self.atom_02_code].type

	def is_marker(self):
		if self.atom_01_parameter == '-C':
			return True
		if self.atom_02_parameter == '+N':
			return True
		return False

