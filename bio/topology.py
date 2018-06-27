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
from bio.parameters99ff import Parameters
from bio.pdb import Molecule
from bio.pdb import Bond
import inf.geometry
import numpy as np
import itertools


class Topology:

	def __init__(self, molecule, parameters):
		self.molecule = molecule
		self.parameters = parameters
		self.bonds_map = {}
		self.interacting = {}
		self.bonds = []
		self.angles = []
		self.improper_dihedrals = []
		self.proper_dihedrals = []
		self.phi_psi_torsions = []
		self.dihedral_map = {}
		# rotations around N - Cα(called Phi, φ) and Cα - C(called Psi, ψ)
		self.backbone_torsions_map = {}
		self.followers_cache = {}
		self.eletrostatic_map = []
		self.solid_side_chain_eletrostatic_map = []

	def _print_missing_types(self):
		for a in self.molecule.atoms:
			if not a.type:
				print('Atom type not found: {} {} {}'.format(a.name, a.res_name, a.res_seq))

	def get_following_atoms(self, atom, excluding_list):
		key = tuple([atom.serial] + list(map(lambda a: a.serial, excluding_list)))
		if key not in self.followers_cache:
			self.followers_cache[key] = self._get_following_atoms(atom, excluding_list)
		return self.followers_cache[key]

	def _get_following_atoms(self, atom, excluding_list):
		result = []
		bonds = []
		for bonded in self.bonds_map[atom.serial]:
			if bonded not in excluding_list:
				bonds.append(bonded)
		self._populate_with_bonds(result, excluding_list, bonds)
		return result

	def _populate_with_bonds(self, result, excluding_list, bonds):
		if not bonds:
			return
		deeper_bonds = []
		for bonded in bonds:
			if bonded not in result and bonded not in excluding_list:
				result.append(bonded)
		for bonded in bonds:
			for deeper in self.bonds_map[bonded.serial]:
				if deeper not in result and deeper not in excluding_list:
					deeper_bonds.append(deeper)
		if deeper_bonds:
			self._populate_with_bonds(result, excluding_list, deeper_bonds)

	def rotate_proper_dihedral(self, dihedral, target_angle, impacted_atoms=None):
		zero_atom = dihedral.atom_02
		axis_atom = dihedral.atom_03
		current_angle = dihedral.get_torsion_angle()

		delta = np.array(zero_atom.point) - np.array([0.0, 0.0, 0.0])
		self.molecule.translate(-delta)

		axis = np.array(axis_atom.point)
		current_psi = np.deg2rad(current_angle)
		wanted_psi = np.deg2rad(target_angle)
		delta_psi = wanted_psi - current_psi
		rotation_matrix = inf.geometry.rotation_matrix(axis, delta_psi)

		if not impacted_atoms:
			impacted_atoms = self.get_following_atoms(dihedral.atom_03, [dihedral.atom_01, dihedral.atom_02, dihedral.atom_03])

		for a in impacted_atoms:
			rotated = np.dot(rotation_matrix, a.point)
			a.move(rotated)

	def change_angle(self, a1, a2, a3, target_angle, impacted_atoms=None):
		delta = np.array(a2.point) - np.array([0.0, 0.0, 0.0])
		self.molecule.translate(-delta)

		current_angle = Angle(a1, a2, a3, None).get_angle()
		delta_angle = np.deg2rad(target_angle - current_angle)

		plane_parallel_i = inf.geometry.to_unity_vector(a2.point, a1.point)
		plane_parallel_j = inf.geometry.to_unity_vector(a2.point, a3.point)
		plane_orthogonal = inf.geometry.calculate_cross_product(plane_parallel_i, plane_parallel_j)
		rotation_matrix = inf.geometry.rotation_matrix(plane_orthogonal, delta_angle)

		if not impacted_atoms:
			impacted_atoms = self.get_following_atoms(a3, [a2])

		for a in impacted_atoms:
			rotated = np.dot(rotation_matrix, a.point)
			a.move(rotated)

	def rotate_backbone(self, dihedral_type, angles):
		for i in range(len(angles)):
			amino_acid = self.molecule.amino_acids[i]
			# TODO change backbone map index type
			if dihedral_type in self.backbone_torsions_map[amino_acid.sequence]:
				dihedral = self.backbone_torsions_map[amino_acid.sequence][dihedral_type]
				self.rotate_proper_dihedral(dihedral, angles[i])

	# CT  C   O            1   120.400    669.440 ;
	def fix_oxygen(self, r1, r2):
		a1 = r2.atoms_map['N']
		a2 = r1.atoms_map['C']
		a3 = r1.atoms_map['O']
		self.change_angle(a1, a2, a3, 122.9, impacted_atoms=[a3])

		a1 = r1.atoms_map['CA']
		a2 = r1.atoms_map['C']
		a3 = r1.atoms_map['O']
		self.change_angle(a1, a2, a3, 120.4, impacted_atoms=[a3])

		if r2.is_last():
			a1 = r2.atoms_map['CA']
			a2 = r2.atoms_map['C']
			a3 = r2.atoms_map['OXT']
			self.change_angle(a1, a2, a3, 120.4, impacted_atoms=[a3])

	# CT  C   N            1   116.600    585.760 ; AA general
	def fix_nitrogen(self, r1, r2):
		a1 = r1.atoms_map['CA']
		a2 = r1.atoms_map['C']
		a3 = r2.atoms_map['N']
		self.change_angle(a1, a2, a3, 116.6)

	# C   N   CT           1   121.900    418.400 ; AA general
	def fix_alpha_carbon(self, r1, r2):
		a1 = r1.atoms_map['C']
		a2 = r2.atoms_map['N']
		a3 = r2.atoms_map['CA']
		self.change_angle(a1, a2, a3, 121.9)

	def fix_amine_hydrogen(self, r1, r2):
		if 'H' not in r2.atoms_map:
			return
		a1 = r1.atoms_map['C']
		a2 = r2.atoms_map['N']
		a3 = r2.atoms_map['H']
		self.change_angle(a1, a2, a3, 120, impacted_atoms=[a3])

	# CT  N   H            1   118.040    418.400 ; new99 general,     changed based on NMA nmodes
	def fix_amine_hydrogen_angle_to_ca(self, r1):
		if 'H' not in r1.atoms_map and 'H1' not in r1.atoms_map:
			return
		a1 = r1.atoms_map['CA']
		a2 = r1.atoms_map['N']
		if r1.is_first():
			h2 = r1.atoms_map['H2']
			h3 = r1.atoms_map['H3']
			self.change_angle(a1, a2, h2, 118.04, impacted_atoms=[h2])
			self.change_angle(a1, a2, h3, 118.04, impacted_atoms=[h3])
		else:
			a3 = r1.atoms_map['H']
			self.change_angle(a1, a2, a3, 118.04, impacted_atoms=[a3])

	def fix_amine_hydrogen_torsion(self, r1, r2):
		if 'H' not in r2.atoms_map:
			return
		a1 = r1.atoms_map['O']
		a2 = r1.atoms_map['C']
		a3 = r2.atoms_map['N']
		a4 = r2.atoms_map['H']
		d = Dihedral(a1, a2, a3, a4, None)
		self.rotate_proper_dihedral(d, 180, impacted_atoms=[a4])

	def set_side_chain_torsion(self, residue, angle):
		i = self.molecule.amino_acids.index(residue)
		last_r = self.molecule.amino_acids[i - 1]
		if 'CB' in residue.atoms_map:
			a1 = last_r.atoms_map['C']
			a2 = residue.atoms_map['N']
			a3 = residue.atoms_map['CA']
			a4 = residue.atoms_map['CB']
			side_chain_torsion = Dihedral(last_r.atoms_map['C'], residue.atoms_map['N'], residue.atoms_map['CA'], residue.atoms_map['CB'], None)
			impacted = []
			impacted += self.get_following_atoms(a4, [a3])
			self.rotate_proper_dihedral(side_chain_torsion, angle, impacted_atoms=impacted)

	def read_bonds_from_pdb_file(self):
		already_mapped = []
		ser_map = self.molecule.serial_map
		for line in self.molecule.pdb_file_content.splitlines():
			if line.startswith('CONECT'):
				values = line.split()
				serial_01 = int(values[1])
				for i in range(2, len(values)):
					serial_02 = int(values[i])
					if serial_01 in ser_map and serial_02 in ser_map:
						key_01 = (serial_01, serial_02)
						key_02 = (serial_02, serial_01)
						if key_01 not in already_mapped and key_02 not in already_mapped:
							self.bonds.append(Bond(ser_map[serial_01], ser_map[serial_02], None))
							already_mapped.append(key_01)
							already_mapped.append(key_02)

	def read_bonds_from_ff(self):
		#mount bonds list
		unbond_sulfides = []
		bond_sulfides = []
		for i in range(len(self.molecule.amino_acids)):
			amino_acid = self.molecule.amino_acids[i]
			amino_acid_bonds = self.mount_amino_acid_bonds(amino_acid)
			if not amino_acid.is_last():
				next = self.molecule.amino_acids[i + 1]
				amino_acid_bonds += [Bond(amino_acid.atoms_map['C'], next.atoms_map['N'], self.parameters.get_bond_type('C', 'N'))]
			self.bonds += amino_acid_bonds
			if 'SG' in amino_acid.atoms_map and not 'HG' in amino_acid.atoms_map:
				unbond_sulfides.append(amino_acid.atoms_map['SG'])

		disulfide_bridge = self.parameters.get_bond_type('S', 'S')
		for s1 in unbond_sulfides:
			for s2 in unbond_sulfides:
				if s1.serial != s2.serial and s1.distance_of(s2) < 2.5 and s1 not in bond_sulfides and s2 not in bond_sulfides:
					self.bonds.append(Bond(s1, s2, disulfide_bridge))
					bond_sulfides.append(s1)
					bond_sulfides.append(s2)

		self._print_missing_types()

	def mount_bond_map(self):
		#mount bonds map
		for a in self.molecule.atoms:
			self.bonds_map[a.serial] = []
			self.interacting[a.serial] = {}
		for bond in self.bonds:
			self.bonds_map[bond.atom_01.serial].append(bond.atom_02)
			self.bonds_map[bond.atom_02.serial].append(bond.atom_01)
			self.interacting[bond.atom_01.serial][bond.atom_02.serial] = True
			self.interacting[bond.atom_02.serial][bond.atom_01.serial] = True

	def read_atom_types_from_bonds(self):
		for atom in self.molecule.atoms:
			atom.type = atom.symbol

	def mount_topology(self, read_atom_types_from_bonds=False):
		self.mount_bond_map()

		if read_atom_types_from_bonds:
			self.read_atom_types_from_bonds()

		#mount angles list
		for atom in self.molecule.atoms:
			for i in range(len(self.bonds_map[atom.serial])):
				for j in range(len(self.bonds_map[atom.serial])):
					if j <= i:
						continue
					left_atom = self.bonds_map[atom.serial][i]
					right_atom = self.bonds_map[atom.serial][j]
					angle_type = self.parameters.get_angle_type(left_atom.type, atom.type, right_atom.type)
					if angle_type:
						#print('angle: {} {} {} ({} {} {})'.format(left_atom.type, atom.type, right_atom.type, left_atom.name, atom.name, right_atom.name))
						self.angles.append(Angle(left_atom, atom, right_atom, angle_type))
						self.interacting[left_atom.serial][atom.serial] = True
						self.interacting[left_atom.serial][right_atom.serial] = True
						self.interacting[atom.serial][left_atom.serial] = True
						self.interacting[atom.serial][right_atom.serial] = True
						self.interacting[right_atom.serial][left_atom.serial] = True
						self.interacting[right_atom.serial][atom.serial] = True

		# mount proper dihedral list
		for atom_01 in self.molecule.atoms:
			for atom_02 in self.bonds_map[atom_01.serial]:
				for atom_03 in self.bonds_map[atom_02.serial]:
					if atom_03.serial == atom_01.serial:
						continue
					for atom_04 in self.bonds_map[atom_03.serial]:
						if atom_04.serial == atom_02.serial:
							continue
						#if atom_01.res_name == 'LEU' and atom_01.name == 'N' and atom_02.name == 'CA' and atom_03.name == 'C' and atom_04.name == 'N':
						#	print('asdf')
						dihedral_types = self.parameters.get_proper_dihedral_type(atom_01.type, atom_02.type, atom_03.type, atom_04.type)
						if dihedral_types:
							key = (atom_01, atom_02, atom_03, atom_04)
							yek = (atom_04, atom_03, atom_02, atom_01)
							if key in self.dihedral_map or yek in self.dihedral_map:
								continue
							d = Dihedral(atom_01, atom_02, atom_03, atom_04, dihedral_types)
							self.dihedral_map[key] = d
							self.proper_dihedrals.append(d)

		# mount backbone torsion angles map
		for d in self.proper_dihedrals:
			if d.atom_03.res_seq not in self.backbone_torsions_map:
				self.backbone_torsions_map[d.atom_03.res_seq] = {}
			a1 = d.atom_01.name
			a2 = d.atom_02.name
			a3 = d.atom_03.name
			a4 = d.atom_04.name
			if a1 == 'C' and a2 == 'N' and a3 == 'CA' and a4 == 'C':
				self.backbone_torsions_map[d.atom_03.res_seq]['PHI'] = d
				self.phi_psi_torsions.append(d)
			elif a1 == 'N' and a2 == 'CA' and a3 == 'C' and a4 == 'N':
				self.backbone_torsions_map[d.atom_03.res_seq]['PSI'] = d
				self.phi_psi_torsions.append(d)
			elif a1 == 'CA' and a2 == 'C' and a3 == 'N' and a4 == 'CA':
				self.backbone_torsions_map[d.atom_03.res_seq]['OMEGA'] = d

		#mount improper dihedrals list
		for atom_03 in self.molecule.atoms:
			if atom_03.serial not in self.bonds_map: continue
			bondeds = self.bonds_map[atom_03.serial]
			if len(bondeds) < 3: continue
			for permutation in itertools.permutations(bondeds):
				atom_01 = permutation[0]
				atom_02 = permutation[1]
				atom_04 = permutation[2]
				improper_dihedral_types = self.parameters.get_improper_dihedral_type(atom_01.type, atom_02.type, atom_03.type, atom_04.type)
				if improper_dihedral_types:
					self.improper_dihedrals.append(Dihedral(atom_01, atom_02, atom_03, atom_04, improper_dihedral_types))

		self.mount_electrostatic_map(self.molecule)

	def mount_electrostatic_map(self, m):
		for i in range(len(m.atoms)):
			for j in range(i + 1, len(m.atoms)):
				if m.atoms[j].serial in self.interacting[m.atoms[i].serial]:
					continue  # skip atoms with covalent bond
				self.eletrostatic_map.append((m.atoms[i], m.atoms[j]))
				if m.atoms[i].res_seq != m.atoms[j].res_seq:
					self.solid_side_chain_eletrostatic_map.append((m.atoms[i], m.atoms[j]))

	def add_participant(self, participant):
		self.eletrostatic_map += participant.eletrostatic_map
		self.eletrostatic_map = []
		for i in range(len(self.molecule.atoms)):
			for p_atom in participant.molecule.atoms:
				self.eletrostatic_map.append((self.molecule.atoms[i], p_atom))

	def mount_amino_acid_bonds(self, amino_acid):
		result = []
		atoms_map = amino_acid.force_field_atoms_map
		for ff_res_bond in amino_acid.amino_acid_parameters.bonds:
			if ff_res_bond.atom_01_code == 'C' and ff_res_bond.atom_02_code == 'N':
				# pepitidic bond is mapped in another place
				continue
			if ff_res_bond.atom_01_code in atoms_map and ff_res_bond.atom_02_code in atoms_map:
				ff_bond = self.parameters.get_bond_type(ff_res_bond.atom_01_type, ff_res_bond.atom_02_type)
				if not ff_bond:
					print(' cant find boundtype for: '+ ff_res_bond.atom_01_type + ' ' + ff_res_bond.atom_02_type)
				result.append(Bond(atoms_map[ff_res_bond.atom_01_code], atoms_map[ff_res_bond.atom_02_code], ff_bond))
		return result

	def print_topology(self):
		print('bonds: {} angles: {} proper dihedrals: {} improper dihedrals: {}'.format(len(self.bonds), len(
			self.angles), len(self.proper_dihedrals), len(self.improper_dihedrals)))

	def print_joint_angles(self):
		residues = self.molecule.amino_acids
		for i in range(len(residues)):
			r1 = residues[i]
			caco = Angle(r1.atoms_map['CA'], r1.atoms_map['C'], r1.atoms_map['O'], None)
			print('{}({:2}) CA C O  {:8.3f}'.format(r1.code, r1.sequence, caco.get_angle()))
			r2 = None if i + 1 >= len(residues) else residues[i + 1]
			if r2:
				ocn = Angle(r1.atoms_map['O'], r1.atoms_map['C'], r2.atoms_map['N'], None)
				print('{}({:2}) O C N  {:8.3f}'.format(r1.code, r1.sequence, ocn.get_angle()))
				cnca = Angle(r1.atoms_map['C'], r2.atoms_map['N'], r2.atoms_map['CA'], None)
				print('{}({:2}) C  N CA {:8.3f}'.format(r1.code, r1.sequence, cnca.get_angle()))
				if 'H' in r2.atoms_map:
					cnh = Angle(r1.atoms_map['C'], r2.atoms_map['N'], r2.atoms_map['H'], None)
					print('{}({:2}) C  N H  {:8.3f}'.format(r1.code, r1.sequence, cnh.get_angle()))

	def print_amino_acids(self):
		for amino_acid in self.molecule.amino_acids:
			print('{}({}) {}'.format(amino_acid.code, amino_acid.sequence, len(amino_acid.atoms_map.values()) ))

	def print_bonds(self):
		for b in self.bonds:
			print('{}({}) {}({})'.format(b.atom_01.name, b.atom_01.serial, b.atom_02.name, b.atom_02.serial))

	def print_bonds_per_residue(self):
		count_map = {}
		for amino_acid in self.molecule.amino_acids:
			count_map[amino_acid.sequence] = 0
			for a in amino_acid.atoms_map.values():
				bonds = ''
				for b in self.bonds_map[a.serial]:
					bonds += str(b.serial) + ' '
					count_map[amino_acid.sequence] += 1
				print('{}({}){} {}'.format(a.res_name, a.res_seq, a.name, bonds))
		for k in count_map:
			print('{} {}'.format(k, count_map[k]))

	def print_charges(self):
		for a in self.molecule.atoms:
			print('{} {:7.3f}'.format(a.serial, a.charge))

	def get_phi_angles(self):
		angles = []
		for amino_acid in self.molecule.amino_acids:
			map = self.backbone_torsions_map[amino_acid.sequence]
			angles.append(360 if 'PHI' not in map else map['PHI'].get_torsion_angle())
		return angles

	def get_psi_angles(self):
		angles = []
		for amino_acid in self.molecule.amino_acids:
			map = self.backbone_torsions_map[amino_acid.sequence]
			angles.append(360 if 'PSI' not in map else map['PSI'].get_torsion_angle())
		return angles

	def print_phi_psi_omega(self):
		for amino_acid in self.molecule.amino_acids:
			map = self.backbone_torsions_map[amino_acid.sequence]
			angles = []
			angles.append(None if 'PHI' not in map else map['PHI'].get_torsion_angle())
			angles.append(None if 'PSI' not in map else map['PSI'].get_torsion_angle())
			angles.append(None if 'OMEGA' not in map else map['OMEGA'].get_torsion_angle())
			# print('AMINO ACID — PHI — PSI — OMEGA')
			line = amino_acid.code + ' \t '
			for i in range(len(angles)):
				angle = angles[i]
				if angle:
					line += '{:8.3f}'.format(angle)
				else:
					line += '        '
				if i + 1 < len(angles):
					line += '\t'
			print(line)

	def print_o_c_n_h_dihedral(self):
		residues = self.molecule.amino_acids
		for i in range(len(residues)):
			r1 = residues[i]
			r2 = None if i + 1 >= len(residues) else residues[i + 1]
			if r2 and 'H' in r2.atoms_map:
				d = Dihedral(r1.atoms_map['O'], r1.atoms_map['C'], r2.atoms_map['N'], r2.atoms_map['H'], None)
				print('{}({:2}) O C N H  {:8.3f}'.format(r1.code, r1.sequence, d.get_torsion_angle()))

	def print_C_N_CA_CB(self):
		residues = self.molecule.amino_acids
		for i in range(len(residues)):
			r1 = residues[i]
			r2 = None if i + 1 >= len(residues) else residues[i + 1]
			if r2 and 'CB' in r2.atoms_map:
				d = Dihedral(r1.atoms_map['C'], r2.atoms_map['N'], r2.atoms_map['CA'], r2.atoms_map['CB'], None)
				print('{}({:2}) C N CA CB  {:8.3f}'.format(r2.code, r2.sequence, d.get_torsion_angle()))


	def print_proper_dihedrals(self, amino_acid=None, backbone=None):
		counter = 0
		for d in self.proper_dihedrals:
			if amino_acid and d.atom_01.res_seq != amino_acid.sequence:
				continue
			if backbone and (
					not d.atom_01.is_backbone()
					or not d.atom_02.is_backbone()
					or not d.atom_03.is_backbone()
					or not d.atom_04.is_backbone()):
				continue
			counter += 1
			print('{} - {}({})\t{}({})\t{}({})\t{}({})\t{}({})\t{:7.3f}'.format(
				counter,
				d.atom_01.name, d.atom_01.res_seq,
				d.atom_02.name, d.atom_02.res_seq,
				d.atom_03.name, d.atom_03.res_seq,
				d.atom_04.name, d.atom_04.res_seq,
				d.atom_03.res_name, d.atom_03.res_seq,
				d.get_torsion_angle()
			))


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
		# TODO make sure the cos is bettewn -1 and 1
		cosine_angle = round(cosine_angle, 3)
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


