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
import time, math
import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
import numpy as np

fig = None
ax = None

def init():
	plt.ion()
	global fig
	global ax
	fig = plt.figure()
	#ax = fig.add_subplot(111, projection='3d')
	ax = fig.gca(projection='3d')
	# ax = Axes3D(fig)

def show_and_block(P, lbound, ubound):
	plot_ax(ax, P, lbound, ubound)
	fig.canvas.draw()
	input("Please enter something: ")

def show(P, lbound, ubound):
	plot_ax(ax, P, lbound, ubound)
	fig.canvas.draw()
	plt.cla()

def plot_ax(ax, P, lbound, ubound):
	x = []
	y = []
	z = []
	u = []
	v = []
	w = []
	shift = ubound
	for i in range(len(P)):
		p = P[i]
		x.append(p.position[0] + shift)
		y.append(p.position[1] + shift)
		z.append(p.position[2] + shift)
		# u.append(math.sqrt(math.pow(np.cos(p.direction[1]), 2) + math.pow(np.cos(p.direction[2]), 2) ))
		# v.append(math.sqrt(math.pow(np.cos(p.direction[0]), 2) + math.pow(np.cos(p.direction[2]), 2) ))
		# w.append(math.sqrt(math.pow(np.cos(p.direction[0]), 2) + math.pow(np.cos(p.direction[1]), 2) ))
		u.append(np.cos(p.direction[1]) + np.cos(p.direction[2]))
		v.append(np.cos(p.direction[0]) + np.cos(p.direction[2]))
		w.append(np.cos(p.direction[0]) + np.cos(p.direction[1]))

	x.append(lbound + shift)
	x.append(ubound + shift)
	y.append(lbound + shift)
	y.append(ubound + shift)
	z.append(lbound + shift)
	z.append(ubound + shift)
	u.append(1)
	v.append(1)
	w.append(1)
	u.append(1)
	v.append(1)
	w.append(1)

	# ax.scatter(xs, ys, zs)
	ax.quiver(x, y, z, u, v, w)
	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')

if __name__ == '__main__':
	import inf.particle
	p1 = inf.particle.Particle()
	p2 = inf.particle.Particle()
	noventa = 0.5 * np.pi
	p2.position = [1, 1, 1]
	p2.direction = [noventa, noventa, noventa]
	p2.direction = [0, 0, 0]
	init()
	show_and_block([p1, p2], -2, 2)

