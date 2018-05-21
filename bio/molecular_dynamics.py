#!/usr/bin/env python3
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

import math
from bio.parameters99ff import Parameters
from bio.pdb import Molecule
from bio.pdb import Bond
from bio.pdb import Angle
from bio.pdb import Dihedral
from bio.topology import Topology
import itertools
import inf.geometry

# http://gaussian.com/constants/
# Coulomb's constant
# http://www.ks.uiuc.edu/Research///namd/1.5/ug/node67.html
# http://m.wolframalpha.com/input/?i=332+kcal*+Angstrom+%2F%28mol+q%5E2%29+to+kj*nm%2F%28mol+q%5E2%29
ELETRON_TO_KJ_PER_MOL_CONSTANT = 138.93541024

class ForceField():

	def __init__(self, parameters):
		self.parameters = parameters
		self._clean_energy_values()

	def _clean_energy_values(self):
		self.energy = 0
		self.angles_e = 0
		self.bonds_e = 0
		self.proper_dihedrals_e = 0
		self.improper_dihedrals_e = 0
		self.vdw_e = 0
		self.eletrostatic_e = 0

	def calculate_energy(self, molecule, atom=None):
		self._clean_energy_values()
		self.molecule = molecule
		self.topology = molecule.get_topology()
		# calculate bond energy
		for i in range(len(self.topology.bonds)):
			bond = self.topology.bonds[i]
			if atom and atom not in bond.atoms:
				continue
			k = bond.bond_type.bond_sprint_constant
			r = bond.atom_01.distance_of(bond.atom_02) / 10 #converting from Angstrom to nm
			r_eq = bond.bond_type.equilibrium_bond_length
			self.bonds_e += k * math.pow(r - r_eq, 2)

		# calculate angle energy
		for i in range(len(self.topology.angles)):
			angle = self.topology.angles[i]
			if atom and atom not in angle.atoms:
				continue
			k = angle.angle_type.angle_sprint_constant
			a = angle.get_angle()
			a_eq = angle.angle_type.equilibrium_bond_angle
			self.angles_e += k * math.pow(math.radians(a - a_eq), 2)

		# calculate proper dihedral energy
		# C   N   CT  C     9     180.0      3.55640     2  ;
		# C   N   CT  C     9       0.0      3.34720     1  ;
		# http://www.wolframalpha.com/input/?i=3.34720%2F2+*+(+1+%2B+cos(+x+*+1+-+0+pi))+%2B+3.55640%2F2+*+(+1+%2B+cos(+x+*+2+-+1+pi))
		for dihedral in self.topology.proper_dihedrals:
			if atom and atom not in dihedral.atoms:
				continue
			for dihedral_type in dihedral.dihedral_types:
				Vn = dihedral_type.rotation_barrier_height
				gamma = dihedral_type.barrier_minimum_offset_angle
				n = dihedral_type.frequency_of_barrier
				torsion = dihedral.get_torsion_angle()
				self.proper_dihedrals_e += Vn / 2 * (1 + math.cos(math.radians(n * torsion - gamma)))

		# calculate improper dihedral energy
		for dihedral in self.topology.improper_dihedrals:
			if atom and atom not in dihedral.atoms:
				continue
			for dihedral_type in dihedral.dihedral_types:
				Vn = dihedral_type.rotation_barrier_height
				gamma = dihedral_type.barrier_minimum_offset_angle
				n = dihedral_type.frequency_of_barrier
				oop_angle = dihedral.get_out_of_plane_angle()
				self.improper_dihedrals_e += Vn / 2 * (1 + math.cos(math.radians(n * oop_angle - gamma)))

		# calculate vdw and eletrostatical energy
		if atom:
			for other in self.molecule.atoms:
				if atom.serial == other.serial:
					continue
				if atom.serial in self.topology.interacting[other.serial]:
					continue
				self._calculate_non_bonded(atom, other)
		else:
			for interaction in self.topology.eletrostatic_map:
				self._calculate_non_bonded(interaction[0], interaction[1])

		self.energy = self.bonds_e + self.angles_e + self.proper_dihedrals_e + self.improper_dihedrals_e + self.vdw_e + self.eletrostatic_e
		return self.energy

	def _calculate_non_bonded(self, atom_01, atom_02):
		# electrostatic
		distance = atom_01.distance_of(atom_02)
		distance = distance / 10  # converting from Angstrom to nm
		if atom_01.charge and atom_02.charge:
			self.eletrostatic_e += ELETRON_TO_KJ_PER_MOL_CONSTANT * atom_01.charge * atom_02.charge / distance
		# VDW
		atom_01_type = self.parameters.atom_types[atom_01.type]
		atom_02_type = self.parameters.atom_types[atom_02.type]
		radio_ratio = ((atom_01_type.sigma + atom_02_type.sigma) / 2) / distance
		epsilon_ij = math.sqrt(atom_01_type.epsilon * atom_02_type.epsilon)
		self.vdw_e += epsilon_ij * (radio_ratio ** 12 - 2.0 * (radio_ratio ** 6))

	def print_energy(self):
		print('  bond              :\t{0:.6f}'.format(self.bonds_e))
		print('  angles            :\t{0:.6f}'.format(self.angles_e))
		print('  proper dihedrals  :\t{0:.6f}'.format(self.proper_dihedrals_e))
		print('  improper dihedrals:\t{0:.6f}'.format(self.improper_dihedrals_e))
		print('  vdw               :\t{0:.6f}'.format(self.vdw_e))
		print('  eletrostatic      :\t{0:.6f}'.format(self.eletrostatic_e))
		print('  total             :\t{0:.6f} KJ/mol'.format(self.energy))
