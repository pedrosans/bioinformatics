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
import time
import matplotlib.pyplot as plt

fig = None
ax = None

def init():
	plt.ion()
	global fig
	global ax
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

def show_and_block(P, lbound, ubound):
	plot_ax(ax, P, lbound, ubound)
	fig.canvas.draw()
	raw_input("Please enter something: ")

def show(P, lbound, ubound):
	plot_ax(ax, P, lbound, ubound)
	fig.canvas.draw()
	plt.cla()

def plot_ax(ax, P, lbound, ubound):
	xs = []
	ys = []
	zs = []
	shift = ubound
	for i in range(len(P)):
		p = P[i]
		xs.append(p.position[0] + shift)
		ys.append(p.position[1] + shift)
		zs.append(p.position[2] + shift)
	xs.append(lbound + shift)
	ys.append(lbound + shift)
	zs.append(lbound + shift)
	xs.append(ubound + shift)
	ys.append(ubound + shift)
	zs.append(ubound + shift)
	ax.scatter(xs, ys, zs)
	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')
