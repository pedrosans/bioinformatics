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
#from pyquaternion import Quaternion

class Particle():

	def __init__(self):
		self.position = [0, 0, 0]
		#self.direction = [0, 0, 0]

	def move(self, velocity):
		self.position[0] = self.position[0] + velocity.shift[0]
		self.position[1] = self.position[1] + velocity.shift[1]
		self.position[2] = self.position[2] + velocity.shift[2]
		#self.direction[0] = self.direction[0] + velocity.rotation[0]
		#self.direction[1] = self.direction[1] + velocity.rotation[1]
		#self.direction[2] = self.direction[2] + velocity.rotation[2]

	def translate(self, delta):
		self.position[0] = self.position[0] + delta[0]
		self.position[1] = self.position[1] + delta[1]
		self.position[2] = self.position[2] + delta[2]
		return self

	def truncate(self, lbound, ubound):
		self.position[0] = max(self.position[0] , lbound)
		self.position[0] = min(self.position[0] , ubound)
		self.position[1] = max(self.position[1] , lbound)
		self.position[1] = min(self.position[1] , ubound)
		self.position[2] = max(self.position[2] , lbound)
		self.position[2] = min(self.position[2] , ubound)

	def copy(self):
		copied = Particle()
		copied.position = list(self.position)
		#copied.direction = list(self.direction)
		return copied

	def position_delta(self, reference_particle):
		delta_x = reference_particle.position[0] - self.position[0]
		delta_y = reference_particle.position[1] - self.position[1]
		delta_z = reference_particle.position[2] - self.position[2]
		return [delta_x, delta_y, delta_z]

	def rotation_delta(self, reference_particle):
		delta_x = reference_particle.direction[0] - self.direction[0]
		delta_y = reference_particle.direction[1] - self.direction[1]
		delta_z = reference_particle.direction[2] - self.direction[2]
		return [delta_x, delta_y, delta_z]

	def move_to_random_place(self):
		self.position = [random_position(), random_position(), random_position()]
		#self.direction = [random.random() * 2, random.random() * 2, random.random() * 2]
		return self

	def to_string(self):
		cp = self.position
		#cd = self.direction
		#return '{:+10.5f} {:+10.5f} {:+10.5f} | {:+10.5f} {:+10.5f} {:+10.5f}'.format(cp[0], cp[1], cp[2], cd[0], cd[1], cd[2])
		return '{:+10.5f} {:+10.5f} {:+10.5f}'.format(cp[0], cp[1], cp[2])

class Velocity():

	def __init__(self):
		self.shift = []
		self.rotation = []

	def set_random_speed(self):
		self.shift = [random_velocity(), random_velocity(), random_velocity()]
		self.rotation = [random_rotation(), random_rotation(), random_rotation()]
		return self


def create_particles(number):
	particles = []
	for i in range(number):
		particles.append(Particle().move_to_random_place())
	return particles

def create_velocities(number):
	velocities = []
	for i in range(number):
		velocities.append(Velocity().set_random_speed())
	return velocities

def random_position():
	return random.random() / 100

def random_velocity():
	return 0.01

def random_rotation():
	return random.random() * 0.2