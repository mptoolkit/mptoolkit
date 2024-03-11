#!/usr/bin/env python
# Matrix Product Toolkit http://mptoolkit.qusim.net/
#
# python/mp-matrix-example.py
#
# Copyright (C) 2020 Ian McCulloch <ian@qusim.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Reseach publications making use of this software should include
# appropriate citations and acknowledgements as described in
# the file CITATIONS in the main source directory.
#----------------------------------------------------------------------------
# ENDHEADER

# example program for the mp-matrix command
# this assumes that the wavefunction is named psi.py in the current directory

# import the wavefunction output from mp-matrix
import psi

import numpy as np

# the MPS matrices
print("MPS has {} sites".format(len(psi.MPS)))

for i in range(len(psi.MPS)):
    print("site {} has a local basis dimension of {}".format(i,len(psi.MPS[i])))
    for j in range(len(psi.MPS[i])):
        print("site {} component {} is a matrix of shape {}".format(i,j, np.shape(psi.MPS[i][j])))

print()

# diagonal components of the density matrix
print("the diagonal components of the density matrix are ", np.diag(psi.RHO))

print()

# test matrix orthogonality
for i in range(len(psi.MPS)):
    print("testing left-orthogonality condition of the A-matrices at site {}".format(i))

    # contract over the A-matrices in left-orthogonalization order,
    # E = \sum_s A^{s \dagger} A^s
    E = np.einsum("ski,skj->ij", np.conj(psi.MPS[i]), psi.MPS[i])

    Ortho = np.allclose(E, np.eye(np.shape(E)[0]))
    print("A-matrix is orthogonal: {}".format(Ortho))
    print()
    
