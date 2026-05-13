// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/time-dependent-mpo.h
//
// Copyright (C) 2026 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#if !defined(MPTOOLKIT_MP_ALGORITHMS_TIME_DEPENDENT_MPO_H)
#define MPTOOLKIT_MP_ALGORITHMS_TIME_DEPENDENT_MPO_H

#include "lattice/infinitelattice.h"
#include "mpo/basic_triangular_mpo.h"
#include <complex>
#include <string>
#include <vector>

struct MagnusQuadrature
{
   std::vector<double> Nodes;
   std::vector<double> Weights;
   std::vector<std::vector<double>> CommutatorWeights;
};

void ValidateMagnusOrder(int MagnusOrder);

int DefaultMagnusQuadratureOrder(int MagnusOrder);

int ResolveMagnusQuadratureOrder(int MagnusOrder, int MagnusQuadratureOrder);

MagnusQuadrature GaussLegendreMagnusQuadrature(int Order);

BasicTriangularMPO
TimeDependentHamiltonianMPO(InfiniteLattice const& Lattice,
                            std::string const& HamOperator,
                            std::string const& TimeVar,
                            std::complex<double> Time,
                            std::complex<double> Timestep,
                            int MagnusOrder,
                            int MagnusQuadratureOrder = 0);

#endif
