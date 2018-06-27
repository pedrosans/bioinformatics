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

import math, random, sys
from inf.particle import Particle
from inf.particle import Velocity
try:
	from bio import ui
	ui_imported = True
except ImportError:
	ui_imported = False


class GradientDescent:

	small_value = 0.000000000001

	def __init__(self, fitness_function, debug=False):
		self.debug = debug
		self.fitness_function = fitness_function
		self.iterations = 10
		self.dimensions_number = 3
		self.precision = 0.01
		self.current_step = self.small_step = 0.15
		self.smallest_step = self.small_step / 4
		self.best_particle = self.best_fitness = self.gradient = None

	def _clean(self):
		self.current_step = self.small_step
		self.best_particle = Particle(dimensions_number=self.dimensions_number)
		self.best_fitness = self.fitness_function(self.best_particle)
		self.gradient = self.calculate_gradient()

	def run(self):
		self._clean()
		# print('{:03}\t{:10.7f}'.format(-1, self.best_fitness))
		for j in range(self.iterations):
			# print('gradient: {:6.2f} {:6.2f} {:6.2f}'.format(self.gradient[0], self.gradient[1], self.gradient[2]))
			direction = -1
			displacement = []
			for i in range(self.dimensions_number):
				displacement.append(self.gradient[i] * self.current_step * direction)

			point = self.best_particle.copy().translate(displacement)
			current_fitness = self.fitness_function(point)
			delta_fitness = self.best_fitness - current_fitness
			if current_fitness < self.best_fitness:
				self.best_fitness = current_fitness
				self.best_particle = point
				self.gradient = self.calculate_gradient()
			if self.debug:
				print('{:03}\t{:10.7f}\tDelta: {:10.7f}\tStep: {:10.7f}\tSample disp: {:10.7f}'.format(j, self.best_fitness, delta_fitness, self.current_step, displacement[0]))
			if delta_fitness < 0 and self.current_step > self.smallest_step:
				# going in the wrong direction, possible overshoot
				self.current_step = self.current_step / 2
				continue
			if self.precision and delta_fitness < self.precision:
				break

	def calculate_gradient(self):
		deltas = []
		for i in range(self.dimensions_number):
			direction = [0] * self.dimensions_number
			direction[i] = GradientDescent.small_value
			deltas.append(self.best_particle.copy().translate(direction))
		gradient = []
		for i in range(self.dimensions_number):
			partial_slope = (self.fitness_function(deltas[i]) - self.best_fitness) / GradientDescent.small_value
			gradient.append(partial_slope)
		return gradient

	def print_result(self):
		print('Fitness: {} '.format(self.best_fitness))


PLOT = False
if PLOT and ui_imported:
	ui.init()


class Pso:

	def __init__(self, fitness_function, dimensions_number=3):
		self.fitness_function = fitness_function
		self.dimensions_number = dimensions_number
		self.lower_bound = self.upper_bound = None
		self.personal_best_retention = 0.7
		self.global_best_retention = 0.3
		self.inertia = 0.8
		self.iterations = 40
		self.number_of_particles = 15
		self.best_particle = self.g_best_fitness = None

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
		self.best_particle = self.g_best_fitness = None
		self.particles = Particle.create_particles(self.number_of_particles, self.lower_bound, self.upper_bound, dimensions_number=self.dimensions_number)
		self.velocities = Velocity.create_velocities(self.number_of_particles, dimensions_number=self.dimensions_number)

	def calculate_new_velocity(self, v, personal_best_delta, global_best_delta):
		for j in range(self.dimensions_number):
			v.dimensions[j] = (v.dimensions[j] * self.inertia
						+ self.keep(personal_best_delta[j], self.personal_best_retention)
						+ self.keep(global_best_delta[j], self.global_best_retention))

	def run(self):
		self._clean()
		for j in range(self.iterations):
			# Assess fitness
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
			print('{:15.3f}\t{:25.3f}\t'.format(self.g_best_fitness, m), end='')
			print(["{0:7.2f}".format(i) for i in self.particles[0].dimensions[:int(self.dimensions_number/2)]])
			# print(["{0:7.2f}".format(i) for i in self.particles[0].dimensions[int(self.dimensions_number/2):]])

			# Determine how to mutate
			for i in range(self.number_of_particles):
				v = self.velocities[i]
				p = self.particles[i]

				l_best_pos_delta = p.delta(self.p_best_particle[i])
				# if i == 0: print(["{0:7.2f}".format(i) for i in l_best_pos_delta])
				g_best_pos_delta = p.delta(self.best_particle)
				self.calculate_new_velocity(v, l_best_pos_delta, g_best_pos_delta)
			# print(["{0:7.2f}".format(i) for i in self.velocities[0].dimensions])

			# Mutate particles
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


