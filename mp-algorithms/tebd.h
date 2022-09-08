// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/tebd.h
//
// Copyright (C) 2020 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Reseach publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#if !defined(MTOOLKIT_MP_ALGORITHMS_TEBD_H)
#define MTOOLKIT_MP_ALGORITHMS_TEBD_H

#include <string>
#include <vector>
#include <map>
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

      int Order_;
      std::string Description_;
      std::vector<double> a_;
      std::vector<double> b_;
};

extern std::map<std::string, LTSDecomposition> Decompositions;

TruncationInfo
DoTEBD(StateComponent& A, StateComponent& B, RealDiagonalOperator& Lambda,
       double& LogAmplitude, SimpleOperator const& U, StatesInfo const& SInfo);

#endif
