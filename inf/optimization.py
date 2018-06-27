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
#!/usr/bin/env python3

import math, random, sys
from bio.pdb import Molecule
from inf import particle
from inf.particle import Particle
from inf.particle import Velocity
try:
	from bio import ui
	ui_imported = True
except ImportError:
	ui_imported = False

ITERATIONS = 40


class GradientDescent:

	small_value = 0.000000000001

	def __init__(self, fitness_function):
		self.fitness_function = fitness_function
		self.best_particle = None
		self.best_fitness = None
		self.gradient = None
		self.precision = None
		self.precision = 0.01
		self.small_step = 0.00015

	def _clean(self):
		self.small_step = 0.00015
		self.best_particle = Particle()
		self.best_fitness = self.fitness_function(self.best_particle)
		self.gradient = self.calculate_gradient(self.best_particle)

	def run(self):
		self._clean()
		# print('{:03}\t{:10.7f}'.format(-1, self.fitness))
		for j in range(ITERATIONS):
			# print('gradient: {:6.2f} {:6.2f} {:6.2f}'.format(gradient[0], gradient[1], gradient[2]))
			direction = -1
			displacement = [
				self.gradient[0] * self.small_step * direction,
				self.gradient[1] * self.small_step * direction,
				self.gradient[2] * self.small_step * direction
			]
			point = self.best_particle.copy().translate(displacement)
			current_fitness = self.fitness_function(point)
			delta_fitness = self.best_fitness - current_fitness
			if current_fitness < self.best_fitness:
				self.best_fitness = current_fitness
				self.best_particle = point
				self.gradient = self.calculate_gradient(self.best_particle)
			# print('{:03}\t{:10.7f}'.format(j, self.fitness))
			if delta_fitness < 0:
				# going in the wrong direction, possible overshoot
				self.small_step = self.small_step / 2
				continue
			if self.precision and delta_fitness < self.precision:
				break

	def calculate_gradient(self, p):
		# TODO make a constant
		delta_x = p.copy().translate([GradientDescent.small_value, 0, 0])
		delta_y = p.copy().translate([0, GradientDescent.small_value, 0])
		delta_z = p.copy().translate([0, 0, GradientDescent.small_value])
		particle_fitness = self.fitness_function(p)
		x_fitness = self.fitness_function(delta_x) - particle_fitness
		y_fitness = self.fitness_function(delta_y) - particle_fitness
		z_fitness = self.fitness_function(delta_z) - particle_fitness
		return [
				x_fitness / GradientDescent.small_value,
				y_fitness / GradientDescent.small_value,
				z_fitness / GradientDescent.small_value
			]

	def print_result(self):
		print('Fitness: {} '.format(self.best_fitness))


PLOT = False
if PLOT and ui_imported:
	ui.init()


class Pso:

	def __init__(self, fitness_function):
		self.fitness_function = fitness_function
		self.lower_bound = self.upper_bound = None
		self.personal_best_retention = 0.7
		self.global_best_retention = 0.3
		self.inertia = 0.8
		self.iterations = 40
		self.number_of_particles = 15

	def set_bounds(self, lower_bound, upper_bound):
		self.lower_bound = lower_bound
		self.upper_bound = upper_bound
		self.set_view_bounds(lower_bound, upper_bound)

	def set_view_bounds(self, lower_bound, upper_bound):
		self.view_lower_bound = lower_bound
		self.view_upper_bound = upper_bound

	@staticmethod
	def keep(value, proportion):
		return value * min(random.random(), proportion)

	def _clean(self):
		self.p_best_fitness = [None] * self.number_of_particles
		self.p_best_particle = [None] * self.number_of_particles
		self.best_particle = Particle()
		self.g_best_fitness = self.fitness_function(self.best_particle)
		self.particles = Particle.create_particles(self.number_of_particles)
		self.velocities = Velocity.create_velocities(self.number_of_particles)
		# print(self.g_best_fitness)

	def calculate_new_velocity(self, current, personal_best, global_best):
		for i in range(len(current)):
			current[i] = (current[i] * self.inertia
					+ self.keep(personal_best[i], self.personal_best_retention)
					+ self.keep(global_best[i], self.global_best_retention))

	def run(self):
		self._clean()
		for j in range(self.iterations):
			#Assess fitness
			iteration_best = None
			fitness_sum = 0
			for i in range(self.number_of_particles):
				fitness = self.fitness_function(self.particles[i])
				fitness_sum += fitness
				if not iteration_best or iteration_best > fitness:
					iteration_best = fitness
				if not self.p_best_fitness[i] or self.p_best_fitness[i] > fitness:
					self.p_best_fitness[i] = fitness
					self.p_best_particle[i] = self.particles[i].copy()
				if not self.g_best_fitness or self.g_best_fitness > fitness:
					self.g_best_fitness = fitness
					self.best_particle = self.p_best_particle[i].copy()
			m = fitness_sum / self.number_of_particles
			print('{}\t{}'.format(self.g_best_fitness, m))

			#Determine how to mutate
			for i in range(self.number_of_particles):
				v = self.velocities[i]
				p = self.particles[i]

				l_best_pos_delta = p.position_delta(self.p_best_particle[i])
				g_best_pos_delta = p.position_delta(self.best_particle)
				self.calculate_new_velocity(v.shift, l_best_pos_delta, g_best_pos_delta)

				l_best_rotation_delta = p.rotation_delta(self.p_best_particle[i])
				g_best_rotation_delta = p.rotation_delta(self.best_particle)
				self.calculate_new_velocity(v.rotation, l_best_rotation_delta, g_best_rotation_delta)

			#Mutate particles
			for i in range(self.number_of_particles):
				self.particles[i].move(self.velocities[i])
				if self.upper_bound:
					self.particles[i].truncate(self.lower_bound, self.upper_bound)

			if PLOT:
				ui.show(self.particles, self.view_lower_bound, self.view_upper_bound)
		if PLOT:
			ui.show_and_block([self.best_particle], self.view_lower_bound, self.view_upper_bound)

	def print_result(self):
		print(self.best_particle.position)
		print('Fitness: {} '.format(self.g_best_fitness))
		#print('RMSD: {} for molecure {} compared with {}'.format(g_best_rmsd, molecule.pdb_file_location, target.pdb_file_location))


