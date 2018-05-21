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
from bio.pdb import Angle
from bio.pdb import Dihedral
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
		self.dihedral_map = {}

	def _print_missing_types(self):
		for a in self.molecule.atoms:
			if not a.type:
				print('Atom type not found: {} {} {}'.format(a.name, a.res_name, a.res_seq))

	def mount_topology(self):
		#mount bonds list
		for i in range(len(self.molecule.amino_acids)):
			amino_acid = self.molecule.amino_acids[i]
			amino_acid_bonds = self.mount_aminoacid_bonds(amino_acid)
			if not amino_acid.is_last():
				next = self.molecule.amino_acids[i + 1]
				amino_acid_bonds += [Bond(amino_acid.atoms_map['C'], next.atoms_map['N'], self.parameters.get_bond_type('C', 'N'))]
			self.bonds += amino_acid_bonds
		# disulfide_bridge = self.parameters.get_bond_type('S', 'S')
		# self.bonds.append(Bond(self.molecule.atoms[34], self.molecule.atoms[260], disulfide_bridge))
		# self.bonds.append(Bond(self.molecule.atoms[71], self.molecule.atoms[330], disulfide_bridge))
		# self.bonds.append(Bond(self.molecule.atoms[127], self.molecule.atoms[356], disulfide_bridge))

		#mount bonds map
		for a in self.molecule.atoms:
			self.bonds_map[a.serial] = []
			self.interacting[a.serial] = {}
		for bond in self.bonds:
			self.bonds_map[bond.atom_01.serial].append(bond.atom_02)
			self.bonds_map[bond.atom_02.serial].append(bond.atom_01)
			self.interacting[bond.atom_01.serial][bond.atom_02.serial] = True
			self.interacting[bond.atom_02.serial][bond.atom_01.serial] = True

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

		self.eletrostatic_map = []
		for i in range(len(self.molecule.atoms)):
			for j in range(i + 1, len(self.molecule.atoms)):
				if self.molecule.atoms[j].serial in self.interacting[self.molecule.atoms[i].serial]:
					continue  # skip atoms with covalent bond
				self.eletrostatic_map.append((self.molecule.atoms[i], self.molecule.atoms[j]))

	def mount_aminoacid_bonds(self, amino_acid):
		result = []
		for mapped_bond in amino_acid.amino_acid_parameters.bonds:
			atom_01_code = mapped_bond.atom_01_code
			atom_02_code = mapped_bond.atom_02_code
			atom_01_type = mapped_bond.atom_01_type
			atom_02_type = mapped_bond.atom_02_type
			if atom_01_code == 'C' and atom_02_code == 'N':
				# pepitidic bond is mapped in another place
				continue
			if atom_01_code in amino_acid.force_field_atoms_map and atom_02_code in amino_acid.force_field_atoms_map:
				bond_type = self.parameters.get_bond_type(atom_01_type, atom_02_type)
				if not bond_type:
					print(' cant find boundtype for: '+atom_01_type + ' ' + atom_02_type)
				atom_01 = amino_acid.force_field_atoms_map[atom_01_code]
				atom_02 = amino_acid.force_field_atoms_map[atom_02_code]
				atom_01.type = atom_01_type
				atom_02.type = atom_02_type
				result.append(Bond(atom_01, atom_02, bond_type))
		return result

	def print_topology(self):
		print('bonds: {} angles: {} proper dihedrals: {} improper dihedrals: {}'.format(len(self.bonds), len(
			self.angles), len(self.proper_dihedrals), len(self.improper_dihedrals)))

	def print_aminoacids(self):
		for amino_acid in self.molecule.amino_acids:
			print('{}({}) {}'.format(amino_acid.code, amino_acid.sequence, len(amino_acid.atoms_map.values()) ))

	def print_bonds(self):
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

	# rotations around N - Cα(called Phi, φ) and Cα - C(called Psi, ψ)
	def print_phi_psi_omega(self):
		backbone_map = {}
		for d in self.proper_dihedrals:
			if (not d.atom_01.is_backbone()
					or not d.atom_02.is_backbone()
					or not d.atom_03.is_backbone()
					or not d.atom_04.is_backbone()):
				continue
			if d.atom_02.res_seq not in backbone_map:
				backbone_map[d.atom_02.res_seq] = [None, None, None]
			if d.atom_02.name == 'N' and d.atom_03.name == 'CA':
				# phi
				backbone_map[d.atom_02.res_seq][0] = d.get_torsion_angle()
			elif d.atom_02.name == 'CA' and d.atom_03.name == 'C':
				# psi
				backbone_map[d.atom_02.res_seq][1] = d.get_torsion_angle()
			elif d.atom_02.name == 'C' and d.atom_03.name == 'N':
				# omega
				backbone_map[d.atom_02.res_seq][2] = d.get_torsion_angle()

		for amino_acid in self.molecule.amino_acids:
			angles = backbone_map[amino_acid.sequence]
			if not angles[0]:
				print('{} —          — {:8.3f} — {:8.3f}'.format(amino_acid.code, angles[1], angles[2]))
			elif not angles[2]:
				print('{} — {:8.3f} —         —         '.format(amino_acid.code, angles[0]))
			else:
				print('{} — {:8.3f} — {:8.3f} — {:8.3f}'.format(amino_acid.code, angles[0], angles[1], angles[2]))
			# print('AMINO ACID — PHI — PSI — OMEGA')

	def print_side_chain_torsions(self):
		residue_map = {}
		for d in self.proper_dihedrals:
			if (d.atom_01.is_backbone()
					or d.atom_02.is_backbone()
					or d.atom_03.is_backbone()
					or d.atom_04.is_backbone()):
				continue
			if d.atom_02.res_seq not in residue_map:
				residue_map[d.atom_02.res_seq] = []
			residue_map[d.atom_02.res_seq].append(d.get_torsion_angle())

		for amino_acid in self.molecule.amino_acids:
			line = amino_acid.code + ' — '
			if amino_acid.sequence in residue_map:
				angles = residue_map[amino_acid.sequence]
				for i in range(len(angles)):
					line += '{:8.3f}'.format(angles[i])
					if i + 1 < len(angles):
						line += ' — '
			print(line)
		# print('AMINO ACID — PHI — PSI — OMEGA')

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
