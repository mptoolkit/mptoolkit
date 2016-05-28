// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-iorthogonalize.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "matrixproduct/infinitewavefunction.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include "common/environment.h"
#include "interface/inittemp.h"
#include "tensor/tensor_eigen.h"

int main(int argc, char** argv)
{
   if (argc != 2)
   {
      std::cout << "usage: mp-iorthogonalize <wavefunction>\n";
      return 1;
   }

   std::string FName = argv[1];

   pvalue_ptr<InfiniteWavefunction> Psi = pheap::OpenPersistent(FName, mp_pheap::CacheSize());

   orthogonalize(*Psi.mutate());

   std::cout.precision(14);
   std::cout << "Orthogonality fidelity = " << (1.0 - orthogonality_fidelity(*Psi)) << '\n';

   MatrixOperator I = MatrixOperator::make_identity(Psi->Psi.Basis1());
   MatrixOperator J = transfer_from_left(I, Psi->Psi);
   J = delta_shift(J, Psi->shift());
   TRACE(norm_frob(I-J));

   LinearWavefunction PsiR = Psi->Psi;
   MatrixOperator Lambda = Psi->C_right;
   MatrixOperator LambdaInv = delta_shift(InvertDiagonal(Lambda, InverseTol), Psi->shift());
   PsiR.set_front(prod(LambdaInv, PsiR.get_front()));
   PsiR.set_back(prod(PsiR.get_back(), Lambda));

   I = MatrixOperator::make_identity(PsiR.Basis2());
   J = transfer_from_right(I, PsiR);
   J = delta_shift(J, adjoint(Psi->shift()));
   TRACE(norm_frob(I-J));

   pheap::ShutdownPersistent(Psi);
}
