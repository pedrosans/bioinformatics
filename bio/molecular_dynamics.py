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
from bio.topology import Angle
from bio.topology import Dihedral
from bio.topology import Topology
import itertools
import inf.geometry

# http://gaussian.com/constants/
# Coulomb's constant
# http://www.ks.uiuc.edu/Research///namd/1.5/ug/node67.html
# http://m.wolframalpha.com/input/?i=332+kcal*+Angstrom+%2F%28mol+q%5E2%29+to+kj*nm%2F%28mol+q%5E2%29
ELECTRON_TO_KJ_PER_MOL_CONSTANT = 138.93541024


class ForceField:

	# TODO: remove parameters, the molecule one should be used
	def __init__(self, parameters):
		self.parameters = parameters
		self.topology = None
		self.energy = self.bonds_e = self.angles_e = self.proper_dihedrals_e = self.improper_dihedrals_e = self.vdw_e = None
		self._clean_energy_values()

	def _clean_energy_values(self):
		self.energy = 0
		self.angles_e = 0
		self.bonds_e = 0
		self.proper_dihedrals_e = 0
		self.improper_dihedrals_e = 0
		self.vdw_e = 0
		self.eletrostatic_e = 0

	def _calculate_bond_energy(self, atom):
		# calculate bond energy
		for i in range(len(self.topology.bonds)):
			bond = self.topology.bonds[i]
			if atom and atom not in bond.atoms:
				continue
			k = bond.bond_type.bond_sprint_constant
			r = bond.atom_01.distance_of(bond.atom_02) / 10   # converting from Angstrom to nm
			r_eq = bond.bond_type.equilibrium_bond_length
			self.bonds_e += k * math.pow(r - r_eq, 2)

	def _calculate_angle_energy(self, atom, warn_high_energy):
		# calculate angle energy
		for i in range(len(self.topology.angles)):
			angle = self.topology.angles[i]
			if atom and atom not in angle.atoms:
				continue
			k = angle.angle_type.angle_sprint_constant
			a = angle.get_angle()
			a_eq = angle.angle_type.equilibrium_bond_angle
			angle_e = k * math.pow(math.radians(a - a_eq), 2)
			if warn_high_energy and angle_e > 100:
				print('warn: angle {} {}({}-{}) {}'.format(angle.atom_01.name, angle.atom_02.name, angle.atom_02.res_name, angle.atom_02.res_seq, angle.atom_03.name))
			self.angles_e += angle_e

	def calculate_energy(self, topology,
						 atom=None,
						 protect_bonds=False,
						 test_phi_psi_torsions_only=False,
						 test_electrostatic_only=False,
						 warn_high_energy=False):
		self._clean_energy_values()
		self.topology = topology

		if not test_electrostatic_only:
			self._calculate_bond_energy(atom)
			if protect_bonds:
				self.bonds_e = math.pow(self.bonds_e, 2)

		if not test_electrostatic_only:
			self._calculate_angle_energy(atom, warn_high_energy)

		# calculate proper dihedral energy
		# ;i  j   k  l	 func      phase      kd      pn
		# C   N   CT  C     9     180.0      3.55640     2  ;
		# C   N   CT  C     9       0.0      3.34720     1  ;
		# http://www.wolframalpha.com/input/?i=3.34720%2F2+*+(+1+%2B+cos(+x+*+1+-+0+pi))+%2B+3.55640%2F2+*+(+1+%2B+cos(+x+*+2+-+1+pi))
		for dihedral in self.topology.proper_dihedrals:
			if test_electrostatic_only: break
			if atom and atom not in dihedral.atoms:
				continue
			if test_phi_psi_torsions_only and dihedral not in self.topology.phi_psi_torsions:
				continue
			for dihedral_type in dihedral.dihedral_types:
				Vn = dihedral_type.rotation_barrier_height
				gamma = dihedral_type.barrier_minimum_offset_angle
				n = dihedral_type.frequency_of_barrier
				torsion = dihedral.get_torsion_angle()
				self.proper_dihedrals_e += Vn / 2 * (1 + math.cos(math.radians(n * torsion - gamma)))

		# calculate improper dihedral energy
		for dihedral in self.topology.improper_dihedrals:
			if test_electrostatic_only: break
			if atom and atom not in dihedral.atoms:
				continue
			if test_phi_psi_torsions_only:
				continue
			for dihedral_type in dihedral.dihedral_types:
				Vn = dihedral_type.rotation_barrier_height
				gamma = dihedral_type.barrier_minimum_offset_angle
				n = dihedral_type.frequency_of_barrier
				oop_angle = dihedral.get_out_of_plane_angle()
				self.improper_dihedrals_e += Vn / 2 * (1 + math.cos(math.radians(n * oop_angle - gamma)))

		# calculate vdw and electrostatic energy
		if atom:
			for other in topology.molecule.atoms:
				if atom.serial == other.serial:
					continue
				if atom.serial in self.topology.interacting[other.serial]:
					continue
				self._calculate_non_bonded(atom, other, warn_high_energy)
		else:
			for interaction in self.topology.eletrostatic_map:
			# for interaction in self.topology.solid_side_chain_eletrostatic_map:
				self._calculate_non_bonded(interaction[0], interaction[1], warn_high_energy)

		self.energy = self.bonds_e + self.angles_e + self.proper_dihedrals_e + self.improper_dihedrals_e + self.vdw_e + self.eletrostatic_e
		return self.energy

	def _calculate_non_bonded(self, atom_01, atom_02, warn_high_energy):
		# electrostatic
		distance = atom_01.distance_of(atom_02)
		distance = distance / 10  # converting from Angstrom to nm
		if atom_01.charge and atom_02.charge:
			self.eletrostatic_e += ELECTRON_TO_KJ_PER_MOL_CONSTANT * atom_01.charge * atom_02.charge / distance
		# VDW
		atom_01_type = self.parameters.atom_types[atom_01.type]
		atom_02_type = self.parameters.atom_types[atom_02.type]
		radio_ratio = ((atom_01_type.sigma + atom_02_type.sigma) / 2) / distance
		epsilon_ij = math.sqrt(atom_01_type.epsilon * atom_02_type.epsilon)
		bond_vdw = epsilon_ij * (radio_ratio ** 12 - 2.0 * (radio_ratio ** 6))
		if warn_high_energy and bond_vdw > 8000:
			print('warn: vdw {}({}) {}({}) energy: {}'.format(atom_01.name, atom_01.serial, atom_01.res_name, atom_01.res_seq, bond_vdw))
		self.vdw_e += bond_vdw

	def print_energy(self):
		print('  bond              :\t{:15.6f}'.format(self.bonds_e))
		print('  angles            :\t{:15.6f}'.format(self.angles_e))
		print('  proper dihedral   :\t{:15.6f}'.format(self.proper_dihedrals_e))
		print('  improper dihedral :\t{:15.6f}'.format(self.improper_dihedrals_e))
		print('  vdw               :\t{:15.6f}'.format(self.vdw_e))
		print('  electrostatic     :\t{:15.6f}'.format(self.eletrostatic_e))
		print('  total             :\t{:15.6f} KJ/mol'.format(self.energy))
