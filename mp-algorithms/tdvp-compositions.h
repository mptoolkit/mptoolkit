// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/tdvp-compositions.h
//
// Copyright (C) 2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#if !defined(MPTOOLKIT_MP_ALGORITHMS_TDVP_COMPOSITIONS_H)
#define MPTOOLKIT_MP_ALGORITHMS_TDVP_COMPOSITIONS_H

#include "tebd.h"

LTSDecomposition
SymmetricDecomposition(int Order, std::string Description,
                       std::vector<double> a,
                       std::vector<double> b);

LTSDecomposition
LeapfrogDecompositionOdd(int Order, std::string Description, std::vector<double> w);

class Composition
{
   public:
      Composition() : Order(0) {}
      Composition(int Order_, std::string Description_, std::vector<double> Gamma_)
         : Order(Order_), Description(Description_), Gamma(Gamma_)
      {
         CHECK(Gamma.size() % 2 == 0);
         CHECK_CLOSE(std::accumulate(Gamma.begin(), Gamma.end(), 0.0), 1.0);
      }

      int Order;
      std::string Description;
      std::vector<double> Gamma;
};

Composition
ToComposition(std::string Description, LTSDecomposition d);

#endif
