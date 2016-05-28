// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// misc/spin2.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "models/spin-su2.h"
#include "tensor/tensor_eigen.h"
#include "tensor/regularize.h"
#include "mps/state_component.h"

// The projector onto a singlet
//   MatrixOperator P0 = (5.0/12.0)*MS + (-1.0/6.0) * MQ + (5.0/9.0) * I;
//   TRACE(triple_prod(U, P0, herm(U)));

LinearAlgebra::Matrix<std::complex<double> > 
ToMatrix(MatrixOperator const& M)
{
   return LinearAlgebra::transform(M.data(), 
       LinearAlgebra::Trace<LinearAlgebra::Matrix<std::complex<double> > >());
}

int main()
{
   SiteBlock Site = CreateSU2SpinSite(2);

   QuantumNumbers::QNConstructor<QuantumNumbers::SU2> QN(Site["I"].GetSymmetryList());

   SiteOperator S = Site["S"];
   SiteOperator Q = Site["Q"];
   SiteOperator T = Site["T"];
   SiteOperator F = Site["F"];

   MatrixOperator I = ToMatrixOperator(tensor_prod(Site["I"], Site["I"], QN(0)));

   MatrixOperator MS = ToMatrixOperator(-sqrt(3.0) * tensor_prod(S, S, QN(0)));
   MatrixOperator MOld = MS;
   MatrixOperator U = DiagonalizeHermitian(MOld);
   TRACE(triple_prod(U, MS, herm(U)));

   MatrixOperator MQ = ToMatrixOperator(-sqrt(5.0) * tensor_prod(Q, Q, QN(0)));
   TRACE(triple_prod(U, MQ, herm(U)));

   MatrixOperator MT = ToMatrixOperator(-sqrt(7.0) * tensor_prod(T, T, QN(0)));
   TRACE(triple_prod(U, MT, herm(U)));

   MatrixOperator MF = ToMatrixOperator(-sqrt(9.0) * tensor_prod(F, F, QN(0)));
   TRACE(triple_prod(U, MF, herm(U)));

   LinearAlgebra::Matrix<double> MSd = LinearAlgebra::real(ToMatrix(MS));
   LinearAlgebra::Matrix<double> MQd = LinearAlgebra::real(ToMatrix(MQ));
   LinearAlgebra::Matrix<double> MTd = LinearAlgebra::real(ToMatrix(MT));
   LinearAlgebra::Matrix<double> MFd = LinearAlgebra::real(ToMatrix(MF));

   // Assemble the coefficient matrix
   LinearAlgebra::Matrix<double> C(5, 5, 0.0);
   for (int i = 0; i < 5; ++i)
   {
      C(i, 0) = MSd(i,i);
      C(i, 1) = MQd(i,i);
      C(i, 2) = MTd(i,i);
      C(i, 3) = MFd(i,i);
      C(i, 4) = 1;
   }

   TRACE(C);

   LinearAlgebra::Matrix<double> RHS(5,5,0.0);
   RHS(0,0) = 1;
   RHS(1,1) = 1;
   RHS(2,2) = 1;
   RHS(3,3) = 1;
   RHS(4,4) = 1;

   LinearAlgebra::Matrix<double> x = LinearAlgebra::LinearSolve(C, RHS);

   TRACE(x);
   TRACE(C*x);

   double p0 = 1;
   double p1 = 1;
   double p2 = 1;
   double p3 = 1;
   double p4 = 1;

   double Dipole     = (-1/50.0)  * p0 + (-1/20.0) * p1 + (-1/20.0)  * p2 + (0.0)     * p3 + (3/25.0)    * p4;
   double Quadrapole = (-1/105.0) * p0 + (-1/70.0) * p1 + (1/98.0)   * p2 + (4/105.0) * p3 + (-6/245.0)  * p4;
   double Hexapole   = (-1/180.0) * p0 + (0.0)     * p1 + (1/63.0)   * p2 + (-1/72.0) * p3 + (1/280.0)   * p4;
   double Octapole   = (-1/180.0) * p0 + (1/90.0)  * p1 + (-1/126.0) * p2 + (1/360.0) * p3 + (-1/2520.0) * p4;
   double c          = (1/25.0)   * p0 + (3/25.0)  * p1 + (1/5.0)    * p2 + (7/25.0)  * p3 + (9/25.0)    * p4;

   TRACE(Dipole)(Quadrapole)(Hexapole)(Octapole)(c);
   
   LinearAlgebra::Matrix<double> Mat = Dipole*MSd + Quadrapole*MQd + Hexapole*MTd + Octapole*MFd + c*RHS;

   TRACE(Mat);

}
