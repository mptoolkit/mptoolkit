// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-iexpectation-sz2.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
#include "models/bosehubbard-spinless.h"
#include "models/bosehubbard-spinless-u1.h"

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
   OneMinusTransfer(std::complex<double> ScaleFactor, LinearWavefunction const& Psi, QuantumNumber const& QShift,
                    MatrixOperator const& Rho,
                    MatrixOperator const& Identity)
      : Scale_(ScaleFactor), Psi_(Psi), QShift_(QShift), Rho_(Rho), Identity_(Identity) {}

   MatrixOperator operator()(MatrixOperator const& x) const
   {
      MatrixOperator r = x-Scale_*delta_shift(transfer_from_left(x, Psi_), QShift_);
      //      r = delta_shift(r, QShift_);
      if (is_scalar(r.TransformsAs()))
          r -= inner_prod(r, Rho_) * Identity_; // orthogonalize to the identity
      return r;
   }

   std::complex<double> Scale_;
   LinearWavefunction const& Psi_;
   QuantumNumber const& QShift_;
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

PolynomialType
delta_shift(PolynomialType const& In, QuantumNumber const& QShift)
{
   PolynomialType Result(In);
   for (PolynomialType::iterator I = Result.begin(); I != Result.end(); ++I)
   {
      I->second = delta_shift(I->second, QShift);
   }
   return Result;
}

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

PolynomialType
inject_left(PolynomialType const& In,
            LinearWavefunction const& Psi1,
            MPOperator const& Op,
            LinearWavefunction const& Psi2)
{
   std::vector<PolynomialType> Vec(1, In);
   Vec = inject_left(Vec, Psi1, Op, Psi2);
   CHECK_EQUAL(Vec.size(), 1);
   return Vec[0];
}

std::vector<PolynomialType>
inject_right(LinearWavefunction const& Psi1,
             MPOperator const& Op,
             LinearWavefunction const& Psi2,
             std::vector<PolynomialType> const& In);

#if 0
// at finite momentum
typedef std::map<std::complex<double>, PolynomialType> MomentumPolynomialType;

std::vector<MomentumPolynomialType>
inject_left(std::vector<MomentumPolynomialType> const& In,
            LinearWavefunction const& Psi1,
            MPOperator const& Op,
            LinearWavefunction const& Psi2)
{

}
#endif

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
   return delta_shift(C[Column], QShift);
}

struct CompareComplex
{
   typedef std::complex<double> first_argument_type;
   typedef std::complex<double> second_argument_type;
   typedef bool result_type;
   bool operator()(std::complex<double> const& x, std::complex<double> const& y) const
   {
      return (x.real() < y.real()) || (x.real() == y.real() && x.imag() < y.imag());
   }
};

typedef std::map<std::complex<double>, PolynomialType, CompareComplex> MomentumPolynomialType;

MomentumPolynomialType
delta_shift(MomentumPolynomialType const& In, QuantumNumber const& QShift)
{
   MomentumPolynomialType Result(In);
   for (MomentumPolynomialType::iterator I = Result.begin(); I != Result.end(); ++I)
   {
      I->second = delta_shift(I->second, QShift);
   }
   return Result;
}

std::vector<MomentumPolynomialType>
delta_shift(std::vector<MomentumPolynomialType> const& In, QuantumNumber const& QShift)
{
   std::vector<MomentumPolynomialType> Result(In);
   for (std::vector<MomentumPolynomialType>::iterator I = Result.begin(); I != Result.end(); ++I)
   {
      *I = delta_shift(*I, QShift);
   }
   return Result;
}

void
add_triple_prod(MomentumPolynomialType& Result, std::complex<double> Factor,
                HermitianProxy<MatrixOperator> const& x,
                MomentumPolynomialType const& E,
                MatrixOperator const& y,
                QuantumNumber const& qxy,
                QuantumNumber const& qEp)
{
   // loop over momenta
   for (MomentumPolynomialType::const_iterator K = E.begin(); K != E.end(); ++K)
   {
      // loop over degrees of the polynomial
      for (PolynomialType::const_iterator D = K->second.begin(); D != K->second.end(); ++D)
      {
         Result[K->first][D->first] += Factor * triple_prod(x, D->second, y, qxy, qEp);
      }
   }
}

std::vector<MomentumPolynomialType>
operator_prod(HermitianProxy<OperatorComponent> const& M,
              HermitianProxy<StateComponent> const& A,
              std::vector<MomentumPolynomialType> const& E,
              StateComponent const& B)
{
   //   DEBUG_PRECONDITION_EQUAL(M.base().LocalBasis2(), A.base().Basis1());
   DEBUG_PRECONDITION_EQUAL(M.base().LocalBasis1(), B.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.base().Basis1().size(), E.size());

   std::vector<MomentumPolynomialType> Result(M.base().Basis2().size());

   // Iterate over the components in M, first index
   for (LinearAlgebra::const_iterator<OperatorComponent>::type I = iterate(M.base()); I; ++I)
   {
      // second index in M
      for (LinearAlgebra::const_inner_iterator<OperatorComponent>::type J = iterate(I); J; ++J)
      {
         // Iterate over the irreducible components of M(I,J)
         for (SimpleRedOperator::const_iterator k = J->begin(); k != J->end(); ++k)
         {
            // *k is an irreducible operator.  Iterate over the components of this operator
            for (LinearAlgebra::const_iterator<SimpleOperator>::type R = iterate(*k); R; ++R)
            {
               for (LinearAlgebra::const_inner_iterator<SimpleOperator>::type
                       S = iterate(R); S; ++S)
               {
                  add_triple_prod(Result[J.index2()], herm(*S),
                                  herm(A.base()[S.index1()]),
                                  E[J.index1()],
                                  B[S.index2()],
                                  k->TransformsAs(),
                                  M.base().Basis2()[J.index2()]);
               }
            }
         }
      }
   }
   return Result;
}

std::vector<MomentumPolynomialType>
operator_prod(HermitianProxy<OperatorComponent> const& M,
              HermitianProxy<StateComponent> const& A,
              std::vector<MomentumPolynomialType> const& E,
              StateComponent const& B,
              std::vector<int> const& OutMask,
              std::vector<int> const& InMask)
{
   //   DEBUG_PRECONDITION_EQUAL(M.base().LocalBasis2(), A.base().Basis2());
   DEBUG_PRECONDITION_EQUAL(M.base().LocalBasis1(), B.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.base().Basis1().size(), E.size());

   std::vector<MomentumPolynomialType> Result(M.base().Basis2().size());

   // Iterate over the components in M, first index
   for (LinearAlgebra::const_iterator<OperatorComponent>::type I = iterate(M.base()); I; ++I)
   {
      // skip over masked components
      if (!InMask[I.index()])
         continue;

      // second index in M
      for (LinearAlgebra::const_inner_iterator<OperatorComponent>::type J = iterate(I); J; ++J)
      {
         // skip over masked components
         if (!OutMask[J.index2()])
            continue;

         // Iterate over the irreducible components of M(I,J)
         for (SimpleRedOperator::const_iterator k = J->begin(); k != J->end(); ++k)
         {
            // *k is an irreducible operator.  Iterate over the components of this operator
            for (LinearAlgebra::const_iterator<SimpleOperator>::type R = iterate(*k); R; ++R)
            {
               for (LinearAlgebra::const_inner_iterator<SimpleOperator>::type
                       S = iterate(R); S; ++S)
               {
                  add_triple_prod(Result[J.index2()], herm(*S),
                                  herm(A.base()[S.index1()]),
                                  E[J.index1()],
                                  B[S.index2()],
                                  k->TransformsAs(),
                                  M.base().Basis2()[J.index2()]);
               }
            }
         }
      }
   }
   return Result;
}

std::vector<MomentumPolynomialType>
inject_left(std::vector<MomentumPolynomialType> const& In,
            LinearWavefunction const& Psi1,
            MPOperator const& Op,
            LinearWavefunction const& Psi2)
{
   PRECONDITION_EQUAL(Psi1.size(), Op.size());
   PRECONDITION_EQUAL(Psi1.size(), Psi2.size());
   PRECONDITION_EQUAL(Op.Basis1().size(), In.size());

   LinearWavefunction::const_iterator I1 = Psi1.begin();
   LinearWavefunction::const_iterator I2 = Psi2.begin();
   MPOperator::const_iterator OpIter = Op.begin();

   std::vector<MomentumPolynomialType> E;
   std::vector<MomentumPolynomialType> Result(In);

   while (OpIter != Op.end())
   {
      std::swap(E, Result);
      //      std::vector<MomentumPolynomialType>(OpIter->Basis2().size()).swap(Result);

      Result = operator_prod(herm(*OpIter), herm(*I1), E, *I2);

      ++I1; ++I2; ++OpIter;
   }
   return Result;
}

std::vector<MomentumPolynomialType>
inject_left_mask(std::vector<MomentumPolynomialType> const& In,
                 LinearWavefunction const& Psi1,
                 QuantumNumber const& QShift,
                 MPOperator const& Op,
                 LinearWavefunction const& Psi2,
                 std::vector<std::vector<int> > const& Mask)
{
   PRECONDITION_EQUAL(Psi1.size(), Op.size());
   PRECONDITION_EQUAL(Psi1.size(), Psi2.size());
   PRECONDITION_EQUAL(Op.Basis1().size(), In.size());

   LinearWavefunction::const_iterator I1 = Psi1.begin();
   LinearWavefunction::const_iterator I2 = Psi2.begin();
   MPOperator::const_iterator OpIter = Op.begin();
   std::vector<std::vector<int> >::const_iterator MaskIter = Mask.begin();

   std::vector<MomentumPolynomialType> E;
   std::vector<MomentumPolynomialType> Result(In);

   while (OpIter != Op.end())
   {
      std::swap(E, Result);
      //      std::vector<MomentumPolynomialType>(OpIter->Basis2().size()).swap(Result);

      Result = operator_prod(herm(*OpIter), herm(*I1), E, *I2, *(MaskIter+1), *MaskIter);

      ++I1; ++I2; ++OpIter; ++MaskIter;
   }
   return delta_shift(Result, QShift);
}

// Solve an MPO in the left-handed sense, as x_L * Op = lambda * x_L
// We currently assume there is only one eigenvalue 1 of the transfer operator
MomentumPolynomialType
SolveMPO_Left(LinearWavefunction const& Psi, QuantumNumber const& QShift,
              TriangularOperator const& Op, MatrixOperator const& Rho,
              MatrixOperator const& Identity, bool Verbose = false)
{
   typedef Polynomial<MatrixOperator> PolyType;
   typedef Polynomial<std::complex<double> > ComplexPolyType;
   int Dim = Op.Basis1().size();  // dimension of the MPO

   typedef Polynomial<std::complex<double> > ComplexPolyType;
   typedef std::map<std::complex<double>, ComplexPolyType, CompareComplex> KComplexPolyType;

   typedef std::map<std::complex<double>, PolynomialType, CompareComplex> KPolyType;
   typedef std::vector<KPolyType> KPolyVecType;
   KPolyVecType EMatK(Dim);

   // Initialize the first E matrix
   EMatK[Dim-1][1.0] = PolyType(Identity);

   // solve recursively
   int Col = Dim-2;
   while (Col >= 0)
   {
      if (Verbose)
      {
         std::cerr << "Solving column " << Col << '\n';
      }

      // Generate the next C matrices, C(n) = sum_{j>Col} Op(j,Col) E_j(n)
      KPolyType C;

#if 0
      MPOperator M = extract_lower_column(Op, Col);
      C = inject_left(EMatK, Psi, M, Psi)[0];
#else
      std::vector<std::vector<int> > Mask;
      mask_lower_column(Op, Col, Mask);
      C = inject_left_mask(EMatK, Psi, QShift, Op.data(), Psi, Mask)[Col];
#endif

      // Now do the classification, based on the properties of the diagonal operator
      MPOperator Diag = Op(Col, Col);
      MPOperatorClassification Classification = classify(Diag);
      if (Classification.is_null())
      {
         DEBUG_TRACE("Zero diagonal element")(Col)(Diag);
         // the operator is zero.  In this case the solution is equation (20)
         // E(n) = C(n-1)
         // which expands to, for the degree d component of the polynomial with momentum k,
         // E^k_d = e^{-ik} C^k_d - \sum_{j=d+1}^p (j d) E^k_j

         KPolyType E;
         for (KPolyType::const_iterator I = C.begin(); I != C.end(); ++I) // sum over momenta
         {
            std::complex<double> K = I->first;  // the momentum (complex phase)
            int MaxDegree = I->second.degree();
            for (int i = MaxDegree; i >= 0; --i)
            {
               if (I->second.has_term(i))
                  E[K][i] = conj(K) * I->second[i];
               for (int j = i+1; j <= MaxDegree; ++j)
               {
                  E[K][i] -= double(Binomial(j,i)) * E[K][j];
               }
            }
         }
         EMatK[Col] = E;
      }
      else if (Classification.is_prop_identity())
      {
         // Operator is proportional to the identity
         DEBUG_TRACE("diagonal element proportional to identity")(Col)(Diag)(Classification.factor());

         KComplexPolyType EParallel;  // components parallel to the identity at momentum factor(), may be zero

         if (Classification.is_string())
         {
            // diagonal element is the identity, up to a unitary factor
            DEBUG_TRACE("Unit diagonal element")(Col)(Classification.factor());

            // decompose C into components parallel and perpendicular to the identity
            // The only part we have to care about is a component with the same momentum as our unit operator
            for (KPolyType::iterator Ck = C.begin(); Ck != C.end(); ++Ck) // sum over momenta
            {
               ComplexPolyType CParallel;
               std::complex<double> K = Ck->first;  // the momentum (complex phase)
               for (PolynomialType::iterator I = Ck->second.begin(); I != Ck->second.end(); ++I)
               {
                  if (!is_scalar(I->second.TransformsAs()))
                     continue;

                  std::complex<double> Overlap = inner_prod(I->second, Rho);
                  I->second -= Overlap*Identity;
                  if (norm_frob(Overlap) > 1E-16)
                     CParallel[I->first] = Overlap;
               }

               // Is this the same momentum as our unit operator?
               if (norm_frob(K - Classification.factor()) < 1E-12)
               {
                  DEBUG_TRACE("Component at equal momenta")(K);
                  // same momenta, these components diverge
                  for (int m = CParallel.degree(); m >= 0; --m)
                  {
                     if (CParallel.has_term(m))
                     {
                        EParallel[Classification.factor()][m+1] = conj(Classification.factor()) * CParallel[m]; // include momentum
                        for (int k = m+2; k <= CParallel.degree()+1; ++k)
                        {
                           if (EParallel[Classification.factor()].has_term(k))
                           {
                              EParallel[Classification.factor()][m+1] -= double(Binomial(k,m))
                                 * EParallel[Classification.factor()][k];
                           }
                        }
                        EParallel[Classification.factor()][m+1] *= 1.0 / (1.0 + m);
                     }
                  }
               }
               else
               {
                  for (int m = CParallel.degree(); m >= 0; --m)
                  {
                     // different momenta, we get a contribution both at K and factor()
                     if (CParallel.has_term(m))
                     {
                        std::complex<double> Term = CParallel[m] / (K - Classification.factor());
                        EParallel[K][m] += Term;
                        EParallel[Classification.factor()][m] -= Term;
                     }
                  }
               }
            }
         }

         // Now the remaining components

         // Components perpendicular to the identity satisfy equation (24)
         KPolyType E;
         for (KPolyType::const_iterator I = C.begin(); I != C.end(); ++I) // sum over momenta
         {
            std::complex<double> K = I->first;  // the momentum (complex phase)
            DEBUG_TRACE("Momentum")(K);
            for (int m = I->second.degree(); m >= 0; --m)
            {
               DEBUG_TRACE("degree")(m);
               MatrixOperator Rhs = conj(K) * I->second[m];
               for (int k = m+1; k <= I->second.degree(); ++k)
               {
                  Rhs -= double(Binomial(k,m)) * E[K][k];
               }

               // orthogonalize Rhs against the identity again, which is a kind of
               // DGKS correction
               if (is_scalar(Rhs.TransformsAs()))
                  Rhs -= inner_prod(Rho, Rhs) * Identity;

               double RhsNorm2 = norm_frob_sq(Rhs);
               RhsNorm2 = RhsNorm2 / (Rhs.Basis1().total_degree()*Rhs.Basis2().total_degree());
               //TRACE(RhsNorm2);
               if (RhsNorm2 > 1E-22)
               {
                  E[K][m] = LinearSolve(OneMinusTransfer(K*conj(Classification.factor()), Psi, QShift, Rho, Identity), Rhs);
                  // do another orthogonalization
                  if (is_scalar(E[K][m].TransformsAs()))
                     E[K][m] -= inner_prod(Rho, E[K][m]) * Identity;
               }
            }
         }

         // Reinsert the components parallel to the identity
         for (KComplexPolyType::const_iterator I = EParallel.begin(); I != EParallel.end(); ++I)
         {
            for (ComplexPolyType::const_iterator J = I->second.begin(); J != I->second.end(); ++J)
            {
               E[I->first][J->first] += J->second * Identity;
            }
         }

         EMatK[Col] = E;
      }

      --Col;
   }

   return EMatK[0];
}

void remove_redundant(OperatorComponent& Op);

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

#if 1
   SiteBlock Site = CreateSpinSite(0.5);
   double J = 1.0;
   double Lambda = 1.0;
   TriangularOperator Op;
   Op =  TriangularOneSite(Site["Sz"]);
   Op = Op * Op;


   //   Op = J * 4.0 * TriangularTwoSite(Site["Sz"], Site["Sz"])
   //      + Lambda * 2.0 * TriangularOneSite(Site["Sx"]);
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

#if 0
   SiteBlock Site = CreateSpinSite(1);
   TriangularOperator Op;
   double const k = math_const::pi*0.75;
   Op = TriangularOneSite(Site["Sp"], k) * TriangularOneSite(Site["Sm"], -k);

   double J =  1.0;
   double Jz = 1;
   TriangularOperator H = J * 0.5 * (TriangularTwoSite(Site["Sp"], Site["Sm"])
                                     + TriangularTwoSite(Site["Sm"], Site["Sp"]))
      + Jz * TriangularTwoSite(Site["Sz"], Site["Sz"]);

   TriangularOperator HB = TriangularOneSite(Site["Sp"], k) * H * TriangularOneSite(Site["Sm"], -k);

   //Op = HB;
#endif

#if 0
   // These parameters are defined in this form to match the Lieb-Liniger model, U = m g dz / hbar^2, energy scaled by hbar^2 / 2 m dz^2
   //   std::cout << "Hamiltonian is spinless Bose-Hubbard, U1 symmetry, T=1, U=" << U << ", Nmax=" << NMax << "\n";
   int NMax = 3;
   SiteBlock Site = CreateBoseHubbardSpinlessU1Site(NMax);
   TriangularOperator Op;
#endif

   //   MpOpTriangular Ham;
   //   Ham = -1.0 * TriangularTwoSite(Site["BH"], Site["B"]) - 1.0 * TriangularTwoSite(Site["B"], Site["BH"])
   //      + U * TriangularOneSite(Site["N2"]) + 2.0 * TriangularOneSite(Site["N"]);
   //   HamMPO.push_back(Ham.data());

   MatrixOperator Identity = MatrixOperator::make_identity(Psi.Psi.Basis1());
   MatrixOperator Rho = scalar_prod(Psi.C_right, herm(Psi.C_right));
   LinearWavefunction Phi = Psi.Psi; // no need to bugger around with C_old,C_right

   MatrixOperator LambdaSqrt = SqrtDiagonal(Psi.C_old);
   MatrixOperator LambdaInvSqrt = InvertDiagonal(LambdaSqrt, InverseTol);

   bool Verbose = true;

   Phi.set_front(prod(LambdaInvSqrt, Phi.get_front()));
   Phi.set_back(prod(Phi.get_back(), delta_shift(LambdaSqrt, adjoint(Psi.QShift))));
   Rho = Psi.C_old;
   Identity = Rho;

   //   TRACE(norm_frob(operator_prod(herm(Phi), Identity, Phi) - Identity));
   //   TRACE(norm_frob(operator_prod(Phi, Rho, herm(Phi)) - Rho));

   std::cout.precision(14);

#if 0
   double Delta = 0.001;
   for (double k = 0; k < 1+Delta/2; k += Delta)
   {
      Op = TriangularOneSite(Site["BH"], k*math_const::pi) * TriangularOneSite(Site["B"], -k*math_const::pi);

      TriangularOperator H = J * 0.5 * (TriangularTwoSite(Site["Sp"], Site["Sm"])
                                        + TriangularTwoSite(Site["Sm"], Site["Sp"]))
         + Jz * TriangularTwoSite(Site["Sz"], Site["Sz"]);

      TriangularOperator HSmCom = H * TriangularOneSite(Site["Sm"], -k*math_const::pi) - TriangularOneSite(Site["Sm"], -k*math_const::pi) * H;

      remove_redundant(HSmCom.data().front());

      HSmCom = H*HSmCom - HSmCom*H;
      remove_redundant(HSmCom.data().front());
      std::cout << HSmCom.Basis1().size() << '\n';
#endif

#if 0
      if (HSmCom.Basis1().size() != 12 || k == Delta)
      {
         std::cout << HSmCom << '\n';
      }
#endif

      //      HB =  TriangularOneSite(Site["Sp"], k*math_const::pi) * HSmCom;

      //HB = TriangularOneSite(Site["Sp"], k*math_const::pi)
      // * H * H * TriangularOneSite(Site["Sm"], -k*math_const::pi);

      Op = extend(Op, 2);
      //      HB = extend(HB, 2);

      std::map<std::complex<double>, Polynomial<MatrixOperator>, CompareComplex>
         E = SolveMPO_Left(Phi, Psi.QShift, Op, Rho, Identity, Verbose);
      Polynomial<std::complex<double> > aNorm = ExtractOverlap(E[1.0], Rho);

      //      E = SolveMPO_Left(Phi, Psi.QShift, HB, Rho, Identity, Verbose);
      //      Polynomial<std::complex<double> > aEnergy = ExtractOverlap(E[1.0], Rho);

      //      std::cout << k;
      for (int i = 0; i <= aNorm.degree(); ++i)
      {
         std::cout << "Coefficient of n^" << i << " term is " << aNorm[i].real() << '\n';
      }
      std::cout << std::endl;

      //      std::cout << "at k=" << k << " pi, E = "
      //                << (aEnergy[2].real() / aNorm[1].real()) << "*N + " << (aEnergy[1].real() / aNorm[1].real()) << '\n';
      //   }

   pheap::Shutdown();
}
