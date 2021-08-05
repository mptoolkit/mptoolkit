// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/tdvp.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
// Copyright (C) 2021 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "tdvp.h"
//#include "lanczos-exponential.h"
#include "tensor/tensor_eigen.h"
#include "linearalgebra/eigen.h"
#include "linearalgebra/exponential.h"

template <typename VectorType, typename MultiplyFunctor>
VectorType LanczosExponential(VectorType const& x,
                              MultiplyFunctor MatVecMultiply, int& Iterations,
                              std::complex<double> const& Theta, double& ErrTol)
{
   std::vector<VectorType> v(Iterations+1);
   std::vector<double> Alpha(Iterations+1);
   std::vector<double> Beta(Iterations);
   double const BetaTol = 1e-14;

   double xNorm = norm_frob(x);

   v[0] = x * (1.0/xNorm);

   VectorType w = MatVecMultiply(v[0]);
   Alpha[0] = real(inner_prod(w, v[0]));
   w = w - Alpha[0] * v[0];

   int i;
   for (i = 1; i <= Iterations; ++i)
   {
      Beta[i-1] = norm_frob(w);

      if (Beta[i-1] < BetaTol) {
         DEBUG_TRACE("Beta hit tolerance")(Beta[i-1]);
         break;
      }

      v[i] = w * (1.0/Beta[i-1]);

      w = MatVecMultiply(v[i]);
      Alpha[i] = real(inner_prod(w, v[i]));
      w = w - Alpha[i] * v[i] - Beta[i-1] * v[i-1];
   }

   Iterations = i-1;

   // Construct the tridiagonal matrix T.
   LinearAlgebra::Matrix<double> T(Iterations+1, Iterations+1, 0.0);
   for (int j = 0; j < Iterations; ++j)
   {
      T(j, j) = Alpha[j];
      T(j, j+1) = T(j+1, j) = Beta[j];
   }
   T(Iterations, Iterations) = Alpha[Iterations];
   
   LinearAlgebra::Matrix<std::complex<double>> P = Theta * T;
   LinearAlgebra::Matrix<std::complex<double>> X = LinearAlgebra::Exponentiate(1.0, P);

   // Calculate the time-evolved vector.
   VectorType Out = v[0] * X(0, 0);
   for (int j = 1; j <= Iterations; ++j)
   {
      Out += v[j] * X(j, 0);
   }
   
   // Normalise Out.
   Out *= xNorm;

   return Out;
}

struct HEff1
{
   HEff1(StateComponent const& E_, OperatorComponent const& H_, StateComponent const& F_)
      : E(E_), H(H_), F(F_)
   {
   }

   StateComponent operator()(StateComponent const& x) const
   {
      StateComponent R = operator_prod_inner(H, E, x, herm(F));
      return R;
   }

   StateComponent const& E;
   OperatorComponent const& H;
   StateComponent const& F;
};

struct HEff2
{
   HEff2(StateComponent const& E_, StateComponent const& F_)
      : E(E_), F(F_)
   {
   }

   MatrixOperator operator()(MatrixOperator const& x) const
   {
      MatrixOperator R = operator_prod(E, x, herm(F));
      return R;
   }

   StateComponent const& E;
   StateComponent const& F;
};

TDVP::TDVP(FiniteWavefunctionLeft const& Psi_, BasicTriangularMPO const& Ham_,
           std::complex<double> Timestep_, int MaxIter_, double ErrTol_, int Verbose_)
   : Hamiltonian(Ham_), Timestep(Timestep_), MaxIter(MaxIter_), ErrTol(ErrTol_),
     Verbose(Verbose_)
{
   // Construct the HamMatrix elements.
   if (Verbose > 0)
      std::cout << "Constructing Hamiltonian block operators..." << std::endl;
   H = Hamiltonian.begin();
   HamMatrices.push_left(Initial_E(Ham_, Psi_.Basis1()));
   for (FiniteWavefunctionLeft::const_mps_iterator I = Psi_.begin(); I != Psi_.end(); ++I)
   {
      if (Verbose > 1)
         std::cout << "site " << (HamMatrices.size_left()) << std::endl;
      HamMatrices.push_left(contract_from_left(*H, herm(*I), HamMatrices.left(), *I));
      Psi.push_back(*I);
      ++H;
   }
   HamMatrices.push_right(Initial_F(Ham_, Psi_.Basis2()));
   // Probably no need to incorporate lambda_r(), this is just a normalization
   Psi.set_back(prod(Psi.get_back(), Psi_.lambda_r()));

   // Initialize to the right-most site.
   HamMatrices.pop_left();
   C = Psi.end();
   --C;
   --H;
   Site = Psi.size() - 1;

   LeftStop = 0;
   RightStop = Psi.size() - 1;
}

FiniteWavefunctionLeft
TDVP::Wavefunction() const
{
   return FiniteWavefunctionLeft::Construct(Psi);
}

std::complex<double>
TDVP::Energy() const
{
   // FIXME?
   return inner_prod(HamMatrices.right(),
		     contract_from_left(*H, herm(HamMatrices.left()), *C, HamMatrices.right()));
}

void TDVP::IterateLeft()
{
   // Evolve current site.
   int Iter = MaxIter;
   double Err = ErrTol;

   *C = LanczosExponential(*C, HEff1(HamMatrices.left(), *H, HamMatrices.right()), Iter, 0.5*Timestep, Err);

   // Perform SVD to right-orthogonalise current site.
   MatrixOperator M = ExpandBasis1(*C);
   MatrixOperator U, Vh;
   RealDiagonalOperator D;

   //SingularValueDecompositionKeepBasis2(M, U, D, Vh);
   SingularValueDecomposition(M, U, D, Vh);

   *C = prod(Vh, *C);

   // Update the effective Hamiltonian.
   HamMatrices.push_right(contract_from_right(herm(*H), *C, HamMatrices.right(), herm(*C)));
   
   // Evolve the UD term backwards in time.
   Iter = MaxIter;
   Err = ErrTol;
   MatrixOperator UD = U*D;

   UD = LanczosExponential(UD, HEff2(HamMatrices.left(), HamMatrices.right()), Iter, -0.5*Timestep, Err);

   // Move to the next site.
   --Site;
   --H;
   --C;

   *C = prod(*C, UD);
   *C *= 1.0 / norm_frob(*C);

   HamMatrices.pop_left();
}

void TDVP::EvolveLeftmostSite()
{
   int Iter = MaxIter;
   double Err = ErrTol;

   *C = LanczosExponential(*C, HEff1(HamMatrices.left(), *H, HamMatrices.right()), Iter, Timestep, Err);
}

void TDVP::IterateRight()
{
   // Perform SVD to left-orthogonalise current site.
   MatrixOperator M = ExpandBasis2(*C);
   MatrixOperator U, Vh;
   RealDiagonalOperator D;

   //SingularValueDecompositionKeepBasis1(M, U, D, Vh);
   SingularValueDecomposition(M, U, D, Vh);

   *C = prod(*C, U);

   // Update the effective Hamiltonian.
   HamMatrices.push_left(contract_from_left(*H, herm(*C), HamMatrices.left(), *C));
   
   // Evolve the DVh term backwards in time.
   int Iter = MaxIter;
   double Err = ErrTol;
   MatrixOperator DVh = D*Vh;

   DVh = LanczosExponential(DVh, HEff2(HamMatrices.left(), HamMatrices.right()), Iter, -0.5*Timestep, Err);

   // Move to the next site.
   ++Site;
   ++H;
   ++C;

   *C = prod(DVh, *C);

   *C *= 1.0 / norm_frob(*C);

   HamMatrices.pop_right();

   // Evolve next site.
   Iter = MaxIter;
   Err = ErrTol;

   *C = LanczosExponential(*C, HEff1(HamMatrices.left(), *H, HamMatrices.right()), Iter, 0.5*Timestep, Err);
}

void TDVP::Evolve()
{
   while (Site > LeftStop)
   {
      this->IterateLeft();
      if (Verbose > 1)
      {
         std::cout << "Site=" << Site << std::endl;
      }
   }

   this->EvolveLeftmostSite();

   while (Site < RightStop)
   {
      this->IterateRight();
      if (Verbose > 1)
      {
         std::cout << "Site=" << Site << std::endl;
      }
   }
}
