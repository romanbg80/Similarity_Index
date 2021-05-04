#! /usr/bin/python

"""
A program to compute the index to compare ADPs (Whitten and Spackman 06).

takes two lists of Uijs along with the cell parameters and
performs an orthogonalization, followed by the computation of

100 * (1 - 2^(2/3.)*(det(U1^(-1)*U2^(-1)))^(1/4.) / (det(U1^(-1) + U2^(-1)))^(1/2.))

Where U1 and U2 are the orthogonalized mean square displacement matrices.

"""

from Crystallographic_2 import ADPfrac2cart
from scipy.linalg import inv, det
from numpy import *
import numpy as np


# atoi converts to integer atof to float
# from string import atof

class SimilarityIndex:

    def __init__(self, file_U1, file_U2, cell_1, cell_2):
        number_atoms_file_1 = self._get_number_atoms(file_U1)
        number_atoms_file_2 = self._get_number_atoms(file_U2)
        self.cell_1 = cell_1
        self.cell_2 = cell_2

        self.result_dicts = []

        if number_atoms_file_1 != number_atoms_file_2:
            raise ValueError("The number of ADPs in your input files is different")


        for iatom in range(0, number_atoms_file_1):
            # get names of atoms from input files
            atomname_1 = self._get_atom_name(filename=file_U1, atomnumber=iatom)
            atomname_2 = self._get_atom_name(filename=file_U2, atomnumber=iatom)
            # read matrix from file
            U1 = self._get_Umatrix_from_file_1(filename=file_U1, atomnumber=iatom)
            U2 = self._get_Umatrix_from_file_2(filename=file_U2, atomnumber=iatom)
            index = self._get_overlap(U1=self.matrix_1, U2=self.matrix_2, cell_1=self.cell_1, cell_2=self.cell_2)

            # save result for printing or other use
            self.result_dicts.append({"atomname": atomname_1, "overlapindex": self.overlapindex})

    def _get_number_atoms(self, filename):
        return sum(1 for line in open(filename) if line != '\n')

    def _get_atom_name(self, filename, atomnumber):
        with open(filename, 'r') as f:
            for iline, line in enumerate(f):
                if iline == atomnumber:
                    splitline = line.split()
                    self.atom = splitline[0]
        return self.atom

    def _get_Umatrix_from_file_1(self, filename, atomnumber):
        with open(filename, 'r') as f:
            for iline, line in enumerate(f):
                if iline == atomnumber:
                    splitline = line.split()
                    newarray = splitline[1:]
                    self.matrix_1 = np.array(
                        [[float(newarray[0]), float(newarray[3]), float(newarray[4])],
                         [float(newarray[3]), float(newarray[1]), float(newarray[5])],
                         [float(newarray[4]), float(newarray[5]), float(newarray[2])]])
        #print(self.matrix_1)
        return self.matrix_1

    def _get_Umatrix_from_file_2(self, filename, atomnumber):
        with open(filename, 'r') as f:
            for iline, line in enumerate(f):
                if iline == atomnumber:
                    splitline = line.split()
                    newarray = splitline[1:]
                    self.matrix_2 = np.array(
                        [[float(newarray[0]), float(newarray[3]), float(newarray[4])],
                         [float(newarray[3]), float(newarray[1]), float(newarray[5])],
                         [float(newarray[4]), float(newarray[5]), float(newarray[2])]])
        #print(self.matrix_2)
        return self.matrix_2

    def _get_overlap(self, U1, U2, cell_1, cell_2):

        # orthogonalize:
        U1c = ADPfrac2cart(cell_1, U1)
        U2c = ADPfrac2cart(cell_2, U2)

        # invert:
        U1ci = inv(U1c)
        U2ci = inv(U2c)

        # compute the overlap index:
        overlapindex = 100 * (1 - 2 ** (3 / 2.) * (det(dot(U1ci, U2ci))) ** (1 / 4.) / (det(U1ci + U2ci)) ** (1 / 2.))
        self.overlapindex = "{:.2f}".format(overlapindex)
        print(self.overlapindex)
        return self.overlapindex

    def print_outputfile(self, filename="SI_20_40.txt"):

        with open(filename, 'w') as f:
            for dictionary in self.result_dicts:
                print(dictionary["atomname"], dictionary["overlapindex"], file=f)

# ambient pressure
#cell_1 = (9.9480, 9.9480, 9.9480, 90.00, 90.00, 90.00)
# 1 GPa
#cell_1 = (9.8971, 9.8971, 9.8971, 90.00, 90.00, 90.00)
# 5 GPa
#cell_1 = (9.72296, 9.72296, 9.72296, 90.00, 90.00, 90.00)
# 10 GPa
#cell_1 = (9.5449, 9.5449, 9.5449, 90.00, 90.00, 90.00)
# 20 GPa
cell_1 = (9.2623, 9.2623, 9.2623, 90.00, 90.00, 90.00)
# 40 GPa
cell_2 = (8.90018, 8.90018, 8.90018, 90.00, 90.00, 90.00)

file_U1 = 'ADPs_20.txt'
file_U2 = 'ADPs_40.txt'

s = SimilarityIndex(file_U1, file_U2, cell_1, cell_2)
s.print_outputfile()


