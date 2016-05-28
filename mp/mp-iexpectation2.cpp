// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-iexpectation2.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "mpo/triangular_operator.h"
#include "mps/infinitewavefunction.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include "common/environment.h"
#include "interface/inittemp.h"
#include "mp-algorithms/gmres.h"
#include "mp-algorithms/arnoldi.h"
#include "common/polynomial.h"
#include "tensor/tensor_eigen.h"

#include "models/spin.h"
#include "models/spin-u1.h"
#include "models/spin-su2.h"
#include "models/tj-u1su2.h"
#include "models/tj-u1.h"
#include "models/spinlessfermion-u1.h"
#include "models/kondo-u1su2.h"
#include "models/kondo-u1.h"
#include "models/hubbard-so4.h"

long Binomial(int n, int k)
{
   if (k > n/2)
      k = n-k;     // take advantage of symmetry
   double r = 1.0;
   for (int i = 1; i <= k; ++i)
   {
      r *= double(n-k+i) / double(i);
   }
   return long(r+0.5); // round to nearest
}

struct OneMinusTransfer
{
   OneMinusTransfer(std::complex<double> ScaleFactor, LinearWavefunction const& Psi, MatrixOperator const& Rho,
                    MatrixOperator const& Identity)
      : Scale_(ScaleFactor), Psi_(Psi), Rho_(Rho), Identity_(Identity) {}

   MatrixOperator operator()(MatrixOperator const& x) const
   {
      MatrixOperator r = x-Scale_*transfer_from_left(x, Psi_);
      //      MatrixOperator r = x-Scale_*0.5*(operator_prod(herm(Psi_), x, Psi_) 
      // + operator_prod(Psi_, x, herm(Psi_)));
      r -= inner_prod(r, Rho_) * Identity_; // orthogonalize to the identity
      return r;
   }

   std::complex<double> Scale_;
   LinearWavefunction const& Psi_;
   MatrixOperator const& Rho_;
   MatrixOperator const& Identity_;
};

template <typename Func>
MatrixOperator
LinearSolve(Func F, MatrixOperator Rhs)
{
   MatrixOperator Guess = Rhs;
   int m = 30;
   int max_iter = 10000;
   double tol = 1e-15;
   GmRes(Guess, F, Rhs, m, max_iter, tol, LinearAlgebra::Identity<MatrixOperator>());
   return Guess;
}

template <typename T>
struct EigenPair
{
   T LeftVector;
   T RightVector;
   std::complex<double> Eigenvalue;
};

typedef Polynomial<MatrixOperator> PolynomialType;

// move to linearwavefunction.h
StateComponent 
inject_left(StateComponent const& In, 
            LinearWavefunction const& Psi1, 
            MPOperator const& Op,
            LinearWavefunction const& Psi2)
{
   PRECONDITION_EQUAL(Psi1.size(), Op.size());
   PRECONDITION_EQUAL(Psi1.size(), Psi2.size());
   StateComponent Result = In;
   LinearWavefunction::const_iterator I1 = Psi1.begin();
   LinearWavefunction::const_iterator I2 = Psi2.begin();
   MPOperator::const_iterator OpIter = Op.begin();

   while (OpIter != Op.end())
   {
      Result = operator_prod(herm(*OpIter), herm(*I1), Result, *I2);
      ++I1; ++I2; ++OpIter;
   }
   return Result;
}

StateComponent 
inject_right(LinearWavefunction const& Psi1, 
             MPOperator const& Op,
             LinearWavefunction const& Psi2,
             StateComponent const& In);

// define here.
// Or better as vector<Polynomial<MatrixOperator> > ?
// Probably even better, to do finite momentum at the same time.
std::vector<PolynomialType>
inject_left(std::vector<PolynomialType> const& In, 
            LinearWavefunction const& Psi1, 
            MPOperator const& Op,
            LinearWavefunction const& Psi2)
{
   PRECONDITION_EQUAL(Psi1.size(), Op.size());
   PRECONDITION_EQUAL(Psi1.size(), Psi2.size());

   std::vector<PolynomialType> Result(Op.Basis2().size());
   int MaxDegree = 0;
   for (unsigned i = 0; i < In.size(); ++i)
      MaxDegree = std::max(In[i].degree(), MaxDegree);

   for (int Degree = 0; Degree <= MaxDegree; ++Degree)
   {
      StateComponent E(Op.Basis1(), Psi1.Basis1(), Psi2.Basis1());
      for (unsigned k = 0; k < E.size(); ++k)
      {
         E[k] = In[k][Degree];
      }

      E = inject_left(E, Psi1, Op, Psi2);

      CHECK_EQUAL(E.size(), Result.size());
      for (unsigned i = 0; i < Result.size(); ++i)
      {
         Result[i][Degree] += E[i];
      }
   }
   return Result;
}

std::vector<PolynomialType>
inject_right(LinearWavefunction const& Psi1, 
             MPOperator const& Op,
             LinearWavefunction const& Psi2,
             std::vector<PolynomialType> const& In);

// Calculates result' = C[column] = sum_{j > Column} E[j] * T_Op(j, Column)
// Assumes that E[j] is defined, for j > Column
PolynomialType
MultiplyLeft(std::vector<PolynomialType> const& E, 
             TriangularOperator const& Op, 
             LinearWavefunction const& Psi, 
             QuantumNumber const& QShift, int Column)
{
   CHECK_EQUAL(Op.size(), Psi.size());
   PolynomialType Result;

   // replace this stuff with the inject_left implementation, and extract_column() in TriangularOperator

   MPOperator OpCol = extract_column(Op, Column);
   std::vector<PolynomialType> C = inject_left(E, Psi, OpCol, Psi);
   return C[Column];
}

// Solve an MPO in the left-handed sense, as x_L * Op = lambda * x_L
// We currently assume there is only one eigenvalue 1 of the transfer operator
PolynomialType
SolveMPO_Left(LinearWavefunction const& Psi, QuantumNumber const& QShift,
              TriangularOperator const& Op, MatrixOperator const& Rho,
              MatrixOperator const& Identity, bool Verbose = false)
{
   typedef Polynomial<MatrixOperator> PolyType;
   typedef Polynomial<std::complex<double> > ComplexPolyType;
   int Dim = Op.Basis1().size();  // dimension of the MPO
   std::vector<PolyType> EMat(Dim);  // The E-matrices

   // Initialize the first E matrix
   EMat[Dim-1] = PolyType(Identity);

   // solve recursively
   int Col = Dim-2;
   while (Col >= 0)
   {
      if (Verbose)
      {
         std::cerr << "Solving column " << Col << '\n';
      }

      // Generate the next C matrices, C(n) = sum_{j>Col} Op(j,Col) E_j(n)
      PolyType C;
      MPOperator M = extract_lower_column(Op, Col);
      C = inject_left(EMat, Psi, M, Psi)[0];

#if 0
      for (int j = Col+1; j < Dim; ++j)
      {
         MPOperator const M = Op(j,Col);
         if (!M.is_null())
         {
            C += inject_left(EMat, Psi, M, Psi)[0];
         }
      }
#endif

      // Now do the classification, based on the properties of the diagonal operator
      MPOperator Diag = Op(Col, Col);
      MPOperatorClassification Classification = classify(Diag);
      if (Classification.is_null())
      {
         // the operator is zero.  In this case the solution is simply
         // E(n) = C(n-1)
         DEBUG_TRACE("Zero diagonal element")(Col)(Diag);
         PolyType E;
         int MaxDegree = C.degree();
         for (int i = MaxDegree; i >= 0; --i)
         {
            E[i] = C[i];
            for (int j = i+1; j <= MaxDegree; ++j)
            {
               E[i] -= double(Binomial(j,i)) * E[j];
            }
         }
         EMat[Col] = E;
      }
      else if (Classification.is_prop_identity())
      {
         // Operator is proportional to the identity
         DEBUG_TRACE("diagonal element proportional to identity")(Col)(Diag)(Classification.factor());

         ComplexPolyType EParallel;  // components parallel to the identity, may be zero

         if (Classification.is_identity())
         {
            // diagonal element is the identity
            DEBUG_TRACE("Unit diagonal element")(Col);
         
            // decompose C into components parallel and perpendicular to the identity
            ComplexPolyType CParallel;
            for (PolyType::iterator I = C.begin(); I != C.end(); ++I)
            {
               std::complex<double> Overlap = inner_prod(I->second, Rho);
               I->second -= Overlap*Identity;
               if (norm_frob(Overlap) > 1E-16)
                  CParallel[I->first] = Overlap;
            }

            // components parallel to the identity satisfy equation (23) of the notes
            for (int m = CParallel.degree(); m >= 0; --m)
            {
               EParallel[m+1] = CParallel[m];
               for (int k = m+2; k <= CParallel.degree()+1; ++k)
               {
                  EParallel[m+1] -= double(Binomial(k,m)) * EParallel[k];
               }
               EParallel[m+1] *= 1.0 / (1.0 + m);
            }
         } 

         // Now the remaining components
      
         // Components perpendicular to the identity satisfy equation (24)
         PolyType E;
         for (int m = C.degree(); m >= 0; --m)
         {
            MatrixOperator Rhs = C[m];
            for (int k = m+1; k <= C.degree(); ++k)
            {
               Rhs -= double(Binomial(k,m)) * E[k];
            }
            
            // orthogonalize Rhs against the identity again, which is a kind of
            // DGKS correction
            Rhs -= inner_prod(Rho, Rhs) * Identity;

            double RhsNorm2 = norm_frob_sq(Rhs);
            RhsNorm2 = RhsNorm2 / (Rhs.Basis1().total_degree()*Rhs.Basis2().total_degree());
            //TRACE(RhsNorm2);
            if (RhsNorm2 > 1E-22)
            {
               E[m] = LinearSolve(OneMinusTransfer(Classification.factor(), Psi, Rho, Identity), Rhs);
               // do another orthogonalization
               E[m] -= inner_prod(Rho, E[m]) * Identity;
            }
         }

         // Reinsert the components parallel to the identity
         for (ComplexPolyType::const_iterator I = EParallel.begin(); I != EParallel.end(); ++I)
         {
            E[I->first] += I->second * Identity;
         }

         // finished this column
         EMat[Col] = E;
      }

      --Col;
   }

   return EMat[0];
}

Polynomial<std::complex<double> >
ExtractOverlap(Polynomial<MatrixOperator> const& E, MatrixOperator const& Rho)
{
   Polynomial<std::complex<double> > Overlap;
   for (Polynomial<MatrixOperator>::const_iterator I = E.begin(); I != E.end(); ++I)
   {
      Overlap[I->first] = inner_prod(I->second, Rho);
   }
   return Overlap;
}

int main(int argc, char** argv)
{
   if (argc < 2 || argc > 2)
   {
      std::cout << "usage: mp-iexpectation <psi>\n";
      return 1;
   }

   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pvalue_ptr<InfiniteWavefunction> PsiPtr = pheap::OpenPersistent(argv[1], CacheSize, true);
   InfiniteWavefunction Psi = *PsiPtr;

#if 0
   SiteBlock Site = CreateSpinSite(0.5);
   double J = 1.0;
   double Lambda = 1.0;
   TriangularOperator Op;
   Op = J * 4.0 * TriangularTwoSite(Site["Sz"], Site["Sz"])
      + Lambda * 2.0 * TriangularOneSite(Site["Sx"]);
#endif

#if 0
   SiteBlock Site = CreateSpinSite(0.5);
   double J = -1.0;
   double Jz = 1;
   TriangularOperator Op;
   Op = J * 0.5 * (TriangularTwoSite(Site["Sp"], Site["Sm"])
                   + TriangularTwoSite(Site["Sm"], Site["Sp"]))
      + Jz * TriangularTwoSite(Site["Sz"], Site["Sz"]);
#endif

#if 1
   SiteBlock Site = CreateSpinSite(1);
   TriangularOperator Op;
   double const k = 2*math_const::pi*0.5;
   Op = TriangularOneSite(Site["Sp"], k) * TriangularOneSite(Site["Sm"], -k);

   double J =  1.0;
   double Jz = 1;
   TriangularOperator H = J * 0.5 * (TriangularTwoSite(Site["Sp"], Site["Sm"])
                                     + TriangularTwoSite(Site["Sm"], Site["Sp"]))
      + Jz * TriangularTwoSite(Site["Sz"], Site["Sz"]);

   TriangularOperator HB = TriangularOneSite(Site["Sp"], k) * H * TriangularOneSite(Site["Sm"], -k);

   //Op = HB;
#endif

   TriangularOperator Ham = Op;

   //Op = Op + 1.2732394790856 * TriangularOneSite(Site["I"]);

   TriangularOperator Const = TriangularOneSite(Site["I"]);
   TriangularOperator Const2 = Const*Const;
   TriangularOperator Const4 = Const2*Const2;

   TRACE(Op);
   
   //Op = Op*Const;
   
   //Op = Op*Op;

   //Op = Op - 1.2883390008989e-07 * TriangularOneSite(Site["I"]);
   //Op = Op - 1.6211387711023 * Const2;
   //Op=Op*Op;

   //Op = Op - 6.1807994119478e-06 * TriangularOneSite(Site["I"]);
   //   Op = Op - 2.6280909151706540665585451556 * Const4;
   //Op = Op*Op; // H^8

   //   Op = Const4;

   Op = extend(Op, 2);

   MatrixOperator Identity = MatrixOperator::make_identity(Psi.Psi.Basis1());
   MatrixOperator Rho = scalar_prod(Psi.C_right, herm(Psi.C_right));
   LinearWavefunction Phi = Psi.Psi; // no need to bugger around with C_old,C_right
 
   MatrixOperator LambdaSqrt = SqrtDiagonal(Psi.C_old);
   MatrixOperator LambdaInvSqrt = InvertDiagonal(LambdaSqrt, InverseTol);

   bool Verbose = true;

   Phi.set_front(prod(LambdaInvSqrt, Phi.get_front()));
   Phi.set_back(prod(Phi.get_back(), LambdaSqrt));
   Rho = Psi.C_old;
   Identity = Rho;

   //   TRACE(norm_frob(operator_prod(herm(Phi), Identity, Phi) - Identity));
   //   TRACE(norm_frob(operator_prod(Phi, Rho, herm(Phi)) - Rho));

   std::cout.precision(14);

   Polynomial<MatrixOperator> E = SolveMPO_Left(Phi, Psi.QShift, Op, Rho, Identity, Verbose);
   Polynomial<std::complex<double> > a = ExtractOverlap(E, Rho);
   std::cout << a;

   pheap::Shutdown();
}
