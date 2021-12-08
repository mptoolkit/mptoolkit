// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/transfer.h
//
// Copyright (C) 2021 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "transfer.h"
#include "mps/packunpack.h"
#include "mp-algorithms/arnoldi.h"
#include "common/statistics.h"
#include "wavefunction/operator_actions.h"
#include "linearalgebra/arpack_wrapper.h"

template <typename Func>
struct PackApplyFunc
{
   PackApplyFunc(PackStateComponent const& Pack_, Func f_) : Pack(Pack_), f(f_) {}

   void operator()(std::complex<double> const* In, std::complex<double>* Out) const
   {
      StateComponent x = Pack.unpack(In);
      x = f(x);
      Pack.pack(x, Out);
   }
   PackStateComponent const& Pack;
   Func f;
};

template <typename Func>
PackApplyFunc<Func>
MakePackApplyFunc(PackStateComponent const& Pack_, Func f_)
{
   return PackApplyFunc<Func>(Pack_, f_);
}

std::tuple<std::complex<double>, int, StateComponent>
get_right_transfer_eigenvector(LinearWavefunction const& Psi1, QuantumNumber const& QShift1,
                      LinearWavefunction const& Psi2, QuantumNumber const& QShift2,
                      ProductMPO const& StringOp,
                      double tol, int Verbose)
{
   int ncv = 0;
   int Length = statistics::lcm(Psi1.size(), Psi2.size(), StringOp.size());
   PackStateComponent Pack(StringOp.Basis1(), Psi1.Basis2(), Psi2.Basis2());
   int n = Pack.size();
   //   double tolsave = tol;
   //   int ncvsave = ncv;
   int NumEigen = 1;

   std::vector<std::complex<double> > OutVec;
      LinearAlgebra::Vector<std::complex<double> > LeftEigen =
         LinearAlgebra::DiagonalizeARPACK(MakePackApplyFunc(Pack,
                                                            RightMultiplyOperator(Psi1, QShift1,
                                                                                 StringOp,
                                                                                 Psi2, QShift2, Length, Verbose-1)),
                                          n, NumEigen, tol, &OutVec, ncv, false, Verbose);

   StateComponent LeftVector = Pack.unpack(&(OutVec[0]));

   return std::make_tuple(LeftEigen[0], Length, LeftVector);
}

std::tuple<std::complex<double>, int, StateComponent>
get_left_transfer_eigenvector(LinearWavefunction const& Psi1, QuantumNumber const& QShift1,
                     LinearWavefunction const& Psi2, QuantumNumber const& QShift2,
                     ProductMPO const& StringOp,
                     double tol, int Verbose)
{
   int ncv = 0;
   int Length = statistics::lcm(Psi1.size(), Psi2.size(), StringOp.size());
   PackStateComponent Pack(StringOp.Basis1(), Psi1.Basis2(), Psi2.Basis2());
   int n = Pack.size();
   //   double tolsave = tol;
   //   int ncvsave = ncv;
   int NumEigen = 1;

   std::vector<std::complex<double> > OutVec;
      LinearAlgebra::Vector<std::complex<double> > LeftEigen =
         LinearAlgebra::DiagonalizeARPACK(MakePackApplyFunc(Pack,
                                                            LeftMultiplyOperator(Psi1, QShift1,
                                                                                 StringOp,
                                                                                 Psi2, QShift2, Length, Verbose-1)),
                                          n, NumEigen, tol, &OutVec, ncv, false, Verbose);

   StateComponent LeftVector = Pack.unpack(&(OutVec[0]));

   return std::make_tuple(LeftEigen[0], Length, LeftVector);
}
