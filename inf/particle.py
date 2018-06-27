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
import random, math


class Particle:
	sequence = 1

	def __init__(self, dimensions_number=3):
		self.dimensions_number = dimensions_number
		self.dimensions = [0] * dimensions_number
		self.id = Particle.sequence
		Particle.sequence += 1

	def move(self, velocity):
		for i in range(self.dimensions_number):
			self.dimensions[i] = self.dimensions[i] + velocity.dimensions[i]

	def translate(self, delta):
		for i in range(self.dimensions_number):
			self.dimensions[i] = self.dimensions[i] + delta[i]
		return self

	def truncate(self, lbound, ubound):
		for i in range(self.dimensions_number):
			self.dimensions[i] = max(self.dimensions[i], lbound)
			self.dimensions[i] = min(self.dimensions[i], ubound)

	def copy(self):
		copied = Particle(dimensions_number=self.dimensions_number)
		copied.dimensions = list(self.dimensions)
		return copied

	def delta(self, reference_particle):
		delta = [None] * self.dimensions_number
		for i in range(self.dimensions_number):
			delta[i] = reference_particle.dimensions[i] - self.dimensions[i]
		return delta

	@staticmethod
	def create_particles(number, lbound, ubound, dimensions_number=3):
		particles = []
		for i in range(number):
			particles.append(Particle(dimensions_number=dimensions_number).move_to_random_place(lbound, ubound))
		return particles

	def move_to_random_place(self, lbound, ubound):
		for i in range(self.dimensions_number):
			self.dimensions[i] = random.randrange(lbound, ubound)
		return self


class Velocity:

	def __init__(self, dimensions_number=3):
		self.dimensions_number = dimensions_number
		self.dimensions = [None] * dimensions_number

	def set_random_speed(self):
		for i in range(self.dimensions_number):
			self.dimensions[i] = random.random() * 2 - 1
		return self

	@staticmethod
	def create_velocities(number, dimensions_number=3):
		velocities = []
		for i in range(number):
			velocities.append(Velocity(dimensions_number=dimensions_number).set_random_speed())
		return velocities
