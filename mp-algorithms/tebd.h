// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/tebd.h
//
// Copyright (C) 2020-2022 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#if !defined(MTOOLKIT_MP_ALGORITHMS_TEBD_H)
#define MTOOLKIT_MP_ALGORITHMS_TEBD_H

#include <complex>
#include <string>
#include <vector>
#include <map>
#include "lattice/infinitelattice.h"
#include "mpo/basic_triangular_mpo.h"
#include "mps/state_component.h"
#include "mps/truncation.h"

// generic Lie-Trotter-Suzuki decomposition with two slices, A and B.
class LTSDecomposition
{
   public:
      LTSDecomposition() : Order_(0) {}
      LTSDecomposition(int Order, std::string Description,
                       std::vector<double> a,
                       std::vector<double> b)
         : Order_(Order), Description_(Description), a_(a), b_(b)
      {
         CHECK(a.size() == b.size() || a.size() == b.size()+1);
         CHECK_CLOSE(std::accumulate(a.begin(), a.end(), 0.0), 1.0);
         CHECK_CLOSE(std::accumulate(b.begin(), b.end(), 0.0), 1.0);
      }

      int order() const { return Order_; }
      std::string description() const { return Description_; }
      std::vector<double> const& a() const { return a_; }
      std::vector<double> const& b() const { return b_; }

   private:
      int Order_;
      std::string Description_;
      std::vector<double> a_;
      std::vector<double> b_;
};

// The list of all available decompositions.  See tebd.cpp for the complete list
extern std::map<std::string, LTSDecomposition> Decompositions;

struct TEBDHamiltonianGates
{
   std::vector<std::vector<SimpleOperator>> EvenU;
   std::vector<std::vector<SimpleOperator>> OddU;
   std::vector<SimpleOperator> EvenContinuation;
};

SimpleOperator
TEBDBondIdentity(OperatorComponent const& Left, OperatorComponent const& Right);

SimpleOperator
ExponentiateTEBDBond(std::complex<double> Factor,
                     SimpleOperator const& BondTerm,
                     OperatorComponent const& Left,
                     OperatorComponent const& Right);

BasicTriangularMPO
PrepareFiniteTEBDHamiltonian(BasicTriangularMPO HamMPO, int PsiSize);

std::vector<SimpleOperator>
FiniteTEBDBondHamiltonian(BasicTriangularMPO const& HamMPO);

std::vector<SimpleOperator>
AssembleFiniteTEBDEvenSlice(BasicTriangularMPO const& HamMPO,
                            std::complex<double> SliceTimestep);

std::vector<SimpleOperator>
AssembleFiniteTEBDOddSlice(BasicTriangularMPO const& HamMPO,
                           std::complex<double> SliceTimestep);

TEBDHamiltonianGates
AssembleFiniteTEBDHamiltonian(BasicTriangularMPO const& HamMPO,
                              std::complex<double> Timestep,
                              LTSDecomposition const& decomp);

TEBDHamiltonianGates
AssembleFiniteTimeDependentTEBDHamiltonian(InfiniteLattice const& Lattice,
                                           std::string const& HamOperator,
                                           std::string const& TimeVar,
                                           std::complex<double> StepStart,
                                           std::complex<double> Timestep,
                                           LTSDecomposition const& decomp,
                                           int MagnusOrder,
                                           int MagnusQuadrature,
                                           int PsiSize);

BasicTriangularMPO
PreparePeriodicTEBDHamiltonianUnitCell(BasicTriangularMPO HamMPO,
                                       int Coarsegrain,
                                       int PsiSize);

std::vector<SimpleOperator>
PeriodicTEBDBondHamiltonian(BasicTriangularMPO const& HamMPO);

std::vector<SimpleOperator>
AssemblePeriodicTEBDEvenSlice(BasicTriangularMPO const& HamMPO,
                              std::complex<double> SliceTimestep);

std::vector<SimpleOperator>
AssemblePeriodicTEBDOddSlice(BasicTriangularMPO const& HamMPO,
                             std::complex<double> SliceTimestep);

TEBDHamiltonianGates
AssemblePeriodicTEBDHamiltonian(BasicTriangularMPO const& HamMPO,
                                std::complex<double> Timestep,
                                LTSDecomposition const& decomp);

TEBDHamiltonianGates
AssemblePeriodicTimeDependentTEBDHamiltonian(InfiniteLattice const& Lattice,
                                             std::string const& HamOperator,
                                             std::string const& TimeVar,
                                             std::complex<double> StepStart,
                                             std::complex<double> Timestep,
                                             LTSDecomposition const& decomp,
                                             int MagnusOrder,
                                             int MagnusQuadrature,
                                             int Coarsegrain,
                                             int PsiSize);

// Do a TEBD iteration.
//
// Input is A, B, Lambda
// A,B are in left-canonical form
// Lambda01 Gamma1 Lambda12 Gamma2 Lambda23
// A = Lambda01 Gamma1
// B = Lambda12 Gamma2
// Lambda = Lambda23
//
// On exit, A and B are in left canonical form,
// final Lambda' = Lambda12
//
TruncationInfo
DoTEBD(StateComponent& A, StateComponent& B, RealDiagonalOperator& Lambda,
       double& LogAmplitude, SimpleOperator const& U, StatesInfo const& SInfo);

#endif
