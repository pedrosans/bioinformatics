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
from bio import ui

ITERATIONS = 10


class GradientDescent:

	small_value = 0.000000000001

	def __init__(self, fitness_function):
		self.fitness_function = fitness_function
		self.best_particle = self.point = None
		self.best_fitness = self.fitness = None
		self.precision = None
		self.precision = 0.01
		self.small_step = 0.00015

	def _clean(self):
		self.small_step = 0.00015
		self.point = Particle()
		self.best_particle = self.point.copy()
		self.best_fitness = self.fitness = self.fitness_function(self.point)

	def run(self):
		self._clean()
		# print('{:03}\t{:10.7f}'.format(-1, self.fitness))
		for j in range(ITERATIONS):
			gradient = self.calculate_gradient(self.best_particle)
			# print('gradient: {:6.2f} {:6.2f} {:6.2f}'.format(gradient[0], gradient[1], gradient[2]))
			direction = -1
			displacement = [
				gradient[0] * self.small_step * direction,
				gradient[1] * self.small_step * direction,
				gradient[2] * self.small_step * direction
			]
			self.point = self.best_particle.copy().translate(displacement)
			current_fitness = self.fitness_function(self.point)
			delta_fitness = self.fitness - current_fitness
			self.fitness = current_fitness
			if self.fitness < self.best_fitness:
				self.best_fitness = self.fitness
				self.best_particle = self.point
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
if PLOT:
	ui.init()
NUMBER_OF_PARTICLE = 10
PB = 0.5
GB = 0.001
UPPER_BOUND=5
LOWER_BOUND=-5


class Pso:

	def __init__(self, fitness_function):
		self.fitness_function = fitness_function

	@staticmethod
	def keep(value, proportion):
		return value * min(random.random(), proportion)

	def _clean(self):
		self.p_best_fitness = [None] * NUMBER_OF_PARTICLE
		self.p_best_particle = [None] * NUMBER_OF_PARTICLE
		self.best_particle = Particle()
		self.g_best_fitness = self.fitness_function(self.best_particle)
		self.particles = particle.create_particles(NUMBER_OF_PARTICLE)
		self.velocities = particle.create_velocities(NUMBER_OF_PARTICLE)
		# print(self.g_best_fitness)

	def run(self):
		self._clean()
		for j in range(ITERATIONS):
			#Assess fitness
			iteration_best = None
			fitness_sum = 0
			for i in range(NUMBER_OF_PARTICLE):
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
			m = fitness_sum / NUMBER_OF_PARTICLE
			# print('{}\t{}'.format(self.g_best_fitness, m))

			#Determine how to mutate
			for i in range(NUMBER_OF_PARTICLE):
				v = self.velocities[i]
				p = self.particles[i]
				l_best_pos_delta = p.position_delta(self.p_best_particle[i])
				g_best_pos_delta = p.position_delta(self.best_particle)

				v.shift[0] = (v.shift[0] * 1) + self.keep(l_best_pos_delta[0], PB) + self.keep(g_best_pos_delta[0], GB)
				v.shift[1] = (v.shift[1] * 1) + self.keep(l_best_pos_delta[1], PB) + self.keep(g_best_pos_delta[1], GB)
				v.shift[2] = (v.shift[2] * 1) + self.keep(l_best_pos_delta[2], PB) + self.keep(g_best_pos_delta[2], GB)

			#Mutate particles
			for i in range(NUMBER_OF_PARTICLE):
				self.particles[i].move(self.velocities[i])
				if UPPER_BOUND:
					self.particles[i].truncate(LOWER_BOUND, UPPER_BOUND)

			if PLOT:
				ui.show(PLOT, LOWER_BOUND, UPPER_BOUND)

	def print_result(self):
		print(self.best_particle.position)
		print('Fitness: {} '.format(self.g_best_fitness))
		#print('RMSD: {} for molecure {} compared with {}'.format(g_best_rmsd, molecule.pdb_file_location, target.pdb_file_location))
		if PLOT:
			ui.show_and_block([self.best_particle], LOWER_BOUND, UPPER_BOUND)


