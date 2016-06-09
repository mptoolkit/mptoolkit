// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mpo/test/test_operator_decompose.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "mpo/operator_component.h"
#include "models/spin.h"
#include "models/spin-u1.h"
#include "models/spin-su2.h"
#include "tensor/tensor_exponential.h"

int main()
{
   {
      SiteBlock Site = CreateSpinSite(0.5);
      SimpleOperator H = 0.5 * (tensor_prod(Site["Sp"], Site["Sm"]) + tensor_prod(Site["Sm"], Site["Sp"])) + tensor_prod(Site["Sz"], Site["Sz"]);
      H = Exponentiate(H * std::complex<double>(0.0, 0.1));
      std::pair<OperatorComponent, OperatorComponent> M = decompose_tensor_prod(H, Site.Basis1().Basis(), Site.Basis1().Basis());
      CHECK_EQUAL(M.first.Basis1().size(), 1);
      CHECK_EQUAL(M.second.Basis2().size(), 1);
      SimpleOperator HTest = local_tensor_prod(M.first, M.second)(0,0).scalar();
      CHECK(norm_frob(HTest - H) < std::numeric_limits<double>::epsilon()*100);
      //TRACE(H)(M.first)(M.second);
   }

   {
      SiteBlock Site = CreateU1SpinSite(0.5);
      SimpleOperator H = 0.5 * (tensor_prod(Site["Sp"], Site["Sm"]) + tensor_prod(Site["Sm"], Site["Sp"])) + tensor_prod(Site["Sz"], Site["Sz"]);
      H = Exponentiate(H * std::complex<double>(0.0, 0.1));
      std::pair<OperatorComponent, OperatorComponent> M = decompose_tensor_prod(H, Site.Basis1().Basis(), Site.Basis1().Basis());
      CHECK_EQUAL(M.first.Basis1().size(), 1);
      CHECK_EQUAL(M.second.Basis2().size(), 1);
      SimpleOperator HTest = local_tensor_prod(M.first, M.second)(0,0).scalar();
      CHECK(norm_frob(HTest - H) < std::numeric_limits<double>::epsilon()*100);
   }

   {
      SiteBlock Site = CreateSU2SpinSite(0.5);
      SimpleOperator H = -sqrt(3.0) * tensor_prod(Site["S"], Site["S"], QuantumNumber(Site.GetSymmetryList()));
      H = Exponentiate(H * std::complex<double>(0.0, 0.1));
      std::pair<OperatorComponent, OperatorComponent> M = decompose_tensor_prod(H, Site.Basis1().Basis(), Site.Basis1().Basis());
      CHECK_EQUAL(M.first.Basis1().size(), 1);
      CHECK_EQUAL(M.second.Basis2().size(), 1);
      SimpleOperator HTest = local_tensor_prod(M.first, M.second)(0,0).scalar();
      CHECK(norm_frob(HTest - H) < std::numeric_limits<double>::epsilon()*100);
   }
}
