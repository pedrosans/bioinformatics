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
import numpy
import numpy as np
import math
from numpy import cross, eye, dot
from scipy.linalg import expm, norm

POINT_IDENTITY_MATRIX = eye(3)


def _rotation_matrix(axis, theta):
	return np.array(expm(cross(POINT_IDENTITY_MATRIX, axis / numpy.linalg.norm(axis) * theta)))


def rotation_matrix(axis, theta):
	"""
	Return the rotation matrix associated with counterclockwise rotation about
	the given axis by theta radians.
	Using the Euler-Rodrigues formula:
	from: https://stackoverflow.com/questions/6802577/rotation-of-3d-vector
	"""
	axis = np.asarray(axis)
	axis = axis/math.sqrt(np.dot(axis, axis))
	a = math.cos(theta/2.0)
	b, c, d = -axis*math.sin(theta/2.0)
	aa, bb, cc, dd = a*a, b*b, c*c, d*d
	bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
	return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
					 [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
					 [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


def calculate_torsion(point_i, point_j, point_k, point_l):
	plane_i_parallel = to_unity_vector(point_j, point_i)
	plane_ij_parallel = to_unity_vector(point_j, point_k)
	plane_j_parallel = to_unity_vector(point_k, point_l)

	plane_i_orthogonal = calculate_cross_product(plane_i_parallel, plane_ij_parallel)
	plane_j_orthogonal = -calculate_cross_product(plane_ij_parallel, plane_j_parallel)
	planes_angle_cos = numpy.dot(plane_i_orthogonal, plane_j_orthogonal)
	# TODO make sure the cos is between -1 and 1
	planes_angle_cos = round(planes_angle_cos, 3)
	sign = 1.0 if numpy.dot(plane_i_orthogonal, plane_j_parallel) <= 0.0 else -1.0
	return numpy.degrees(sign * math.acos(planes_angle_cos))


def calculate_out_of_plane_angle(point_i, point_j, point_k, point_l):
	plane_i_parallel_01 = to_unity_vector(point_k, point_i)
	plane_i_parallel_02 = to_unity_vector(point_k, point_j)
	plane_i_perpendicular = to_unity_vector(point_k, point_l)
	plane_i_orthogonal = calculate_cross_product(plane_i_parallel_01, plane_i_parallel_02)
	perpendicular_orthogonal_vectors_angle_cos = numpy.dot(plane_i_orthogonal, plane_i_perpendicular)
	return numpy.degrees(math.asin(perpendicular_orthogonal_vectors_angle_cos))


def to_unity_vector(point_01, point_02):
	delta = numpy.array(point_02) - numpy.array(point_01)
	delta_magnitude = numpy.linalg.norm(delta)
	return delta / delta_magnitude


def calculate_cross_product(point_01, point_02):
	cross = numpy.cross(point_01, point_02)
	return to_unity_vector([0, 0, 0], cross)

