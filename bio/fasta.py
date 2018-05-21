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

from bio.parameters99ff import Parameters
from bio.pdb import Molecule
import numpy as np

code_map = {
	'A': 'ALA',  # Alanine
	'B': 'ASX',  # Aspartic acid (D) or Asparagine (N)
	'C': 'CYS',  # Cysteine
	'D': 'ASP',  # Aspartic acid
	'E': 'GLU',  # Glutamic acid
	'F': 'PHE',  # Phenylalanine
	'G': 'GLY',  # Glycine
	'H': 'HIS',  # Histidine
	'I': 'ILE',  # Isoleucine
	'K': 'LYS',  # Lysine
	'L': 'LEU',  # Leucine
	'M': 'MET',  # Methionine
	'N': 'ASN',  # Asparagine
	'P': 'PRO',  # Proline
	'Q': 'GLN',  # Glutamine
	'R': 'ARG',  # Arginine
	'S': 'SER',  # Serine
	'T': 'THR',  # Threonine
	'V': 'VAL',  # Valine
	'W': 'TRP',  # Tryptophan
	'Y': 'TYR',  # Tyrosine
	'Z': 'GLX',  # Glutamic acid (E) or Glutamine (Q)
}

input = 'NLYIQWLKDGGPSSGRPPPS'

output = Molecule()
pointer = np.array([0, 0, 0])
separator = np.array([1, 0, 0])
i = 0
for a in input:
	i += 1
	abreviation = code_map[a]
	# residue = Molecule(pdb_file_location='/home/pedro/dev/src/bioinf/data/catalog/' + abreviation + '.pdb')
	residue = Molecule(pdb_file_location='data/catalog/' + abreviation + '.pdb')
	n_terminal = residue.amino_acids[0].atoms_map['N']
	c_terminal = residue.amino_acids[0].atoms_map['C']
	delta = np.array(n_terminal.point) - np.array(pointer)
	residue.dislocate(-delta)
	for a in residue.atoms:
		a.res_seq = i
		output.atoms.append(a.copy())
	print(abreviation)
	pointer = np.array(c_terminal.point) + separator
output.write_all('/home/pedro/tmp', 'generated.pdb')
