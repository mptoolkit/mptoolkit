// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mpo/e-optim.cpp
//
// Copyright (C) 2015-2022 Ian McCulloch <ian@qusim.net>
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
//
// Optimized version of the E-matrix operator_prod
//

#include "mpo/operator_component.h"
#include "common/openmp.h"

// note: we could parallelize the construction of the G and H indices over jP

namespace
{

typedef std::complex<double> NumberType;
typedef LinearAlgebra::Matrix<NumberType> MatrixType;

// This version uses a Matrix* directly
struct JMatrixRef
{
   MatrixType const* AH;
   MatrixType const* E;

   JMatrixRef(MatrixType const& AH_, MatrixType const& E_) : AH(&AH_), E(&E_) {}
};

bool
operator<(JMatrixRef const& x, JMatrixRef const& y)
{
   // FIXME: should use std::less here, but in practice it surely doesn't matter
   return x.AH < y.AH || (x.AH == y.AH && x.E < y.E);
}

class JMatrixRefList
{
   public:
      JMatrixRefList() {}

      // we can reserve a maximum size
      void reserve(int MaxSize);

      // returns the index of a G-matrix (which is inserted if it didn't exist previously)
      int Lookup(MatrixType const& AH, MatrixType const& E);

      // do the actual evaluations (in parallel)
      void Evaluate();

      // after the evaluation, return the n'th matrix
      MatrixType const& operator[](int n) const
      {
         return Data[n].second;
      }

   private:
      typedef std::map<JMatrixRef, int> IndexType;
      IndexType Indices;
      std::vector<std::pair<JMatrixRef, MatrixType> > Data;
};

void
JMatrixRefList::reserve(int MaxSize)
{
   Data.reserve(MaxSize);
}

int
JMatrixRefList::Lookup(MatrixType const& AH, MatrixType const& E)
{
   JMatrixRef Index(AH, E);
   IndexType::const_iterator I = Indices.find(Index);
   if (I == Indices.end())
   {
      int r = Data.size();
      Indices[Index] = r;
      Data.push_back(std::make_pair(Index, MatrixType()));
      return r;
   }
   return I->second;
}

void
JMatrixRefList::Evaluate()
{
   for (unsigned n = 0; n < Data.size(); ++n)
   {
      Data[n].second = herm(*Data[n].first.AH) * (*Data[n].first.E);
   }
}

// array of coefficient * index into the GMatrices array
typedef std::vector<std::pair<NumberType, int> > KMatrixDescriptor;

MatrixType EvaluateK(KMatrixDescriptor const& K, JMatrixRefList const& J)
{
   CHECK(!K.empty())("K Matrix descriptor should not be empty!");

   MatrixType Result = K[0].first * J[K[0].second];
   for (unsigned n = 1; n < K.size(); ++n)
   {
      Result += K[n].first * J[K[n].second];
   }
   return Result;
}

struct ElementRec
{
   ElementRec(int s_, int jP_, KMatrixDescriptor const& K_)
      : s(s_), jP(jP_), K(K_) {}

   int s, jP;
   KMatrixDescriptor K;
};

struct OuterIndex
{
   OuterIndex(int i_, int a_, int j_)
      : i(i_), a(a_), j(j_) {}

   bool empty() const { return Components.empty(); }

   // do the final evaluation.  Precondition: !empty()
   MatrixType Evaluate(StateComponent const& A, JMatrixRefList const& J) const;

   int i, a, j;
   std::vector<ElementRec> Components;
};

MatrixType
OuterIndex::Evaluate(StateComponent const& B, JMatrixRefList const& J) const
{
   DEBUG_CHECK(!Components.empty());
#if defined(USE_OPENMP_OPTIM)
   int Sz = Components.size();
   std::vector<MatrixType> Result(Sz);
   //   std::cout << "Size:" << Sz << '\n';
#pragma omp parallel for schedule(dynamic) num_threads(omp::threads_to_use(Sz))
   for (int n = 0; n < Sz; ++n)
   {
      //      std::cout << omp_get_num_threads() << '\n';
      MatrixOperator::const_inner_iterator x = iterate_at(B[Components[n].s].data(), Components[n].jP, j);
      DEBUG_CHECK(x);
      Result[n] = EvaluateK(Components[n].K, J) * (*x);
   }
   return omp::parallel_sum(std::move(Result));
#else
   MatrixOperator::const_inner_iterator x = iterate_at(B[Components[0].s].data(), Components[0].jP, j);
   MatrixType Result = EvaluateK(Components[0].K, J) * (*x);
   for (unsigned n = 1; n < Components.size(); ++n)
   {
      x = iterate_at(B[Components[n].s].data(), Components[n].jP, j);
      DEBUG_CHECK(x);
      Result += EvaluateK(Components[n].K, J) * (*x);
   }
   return Result;
#endif
}

} // namespace

StateComponent
operator_prod(HermitianProxy<OperatorComponent> const& M,
              HermitianProxy<StateComponent> const& A,
              StateComponent const& E,
              StateComponent const& B)
{
   PRECONDITION_EQUAL(M.base().LocalBasis2(), A.base().LocalBasis());
   PRECONDITION_EQUAL(M.base().LocalBasis1(), B.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.base().Basis1(), E.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(A.base().Basis1(), E.Basis1());
   DEBUG_PRECONDITION_EQUAL(E.Basis2(), B.Basis1());

   // firstly, assemble the required H matrix descriptors
   JMatrixRefList JList;

   // this describes how to construct the final matrices
   std::vector<OuterIndex> C;

   // Firstly, iterate over the outer indices - we want to be careful to only
   // add indices that are definitely needed, the MPS might have more sparseness than
   // the minimum imposed by the symmetries.
   StateComponent Result(M.base().Basis2(), A.base().Basis2(), B.Basis2());

   // iterate over all possible output matrix elements E'[a](i, j)
   for (unsigned a = 0; a < M.base().Basis2().size(); ++a)
   {
      for (unsigned i = 0; i < A.base().Basis2().size(); ++i)
      {
         double degree_i = degree(A.base().Basis2()[i]);

         for (unsigned j = 0; j < B.Basis2().size(); ++j)
         {
            if (!is_transform_target(B.Basis2()[j], M.base().Basis2()[a], A.base().Basis2()[i]))
               continue;

            // Now we have enough to construct the output descriptor
            C.push_back(OuterIndex(i, a, j));

            // Iterate over B[s](jP,j)
            for (unsigned s = 0; s < A.base().LocalBasis().size(); ++s)
            {
               // We already know j, unlike the F-matrix case we don't know the leading index
               // so we need to iterate
               for (unsigned jP = 0; jP < B.Basis1().size(); ++jP)
               {
                  MatrixOperator::const_inner_iterator J = iterate_at(B[s].data(), jP, j);
                  if (!J)
                     continue;

                  KMatrixDescriptor KMat;

                  // Iterate over the components in M[a'a]
                  // Unlike the F-matrix case, we don't know the leading index so we need to iterate
                  // over both indices
                  for (unsigned aP = 0; aP < M.base().Basis1().size(); ++aP)
                  {
                     OperatorComponent::const_inner_iterator AA = iterate_at(M.base().data(), aP, a);
                     if (!AA)
                        continue;

                     // Iterate over the irreducible components of M(aP,a)
                     for (SimpleRedOperator::const_iterator k = AA->begin(); k != AA->end(); ++k)
                     {
                        // *k is an irreducible operator.  Iterate over the components of this operator.
                        // Unlike the F-matrix case, we don't know the leading index so we need to iterate
                        for (unsigned sP = 0; sP < M.base().LocalBasis1().size(); ++sP)
                        {
                           SimpleOperator::const_inner_iterator S = iterate_at(k->data(), sP, s);
                           if (!S)
                              continue;

                           // The final index is i' - we only need this if the
                           // element exists in both A.base()[s'][i', i] and
                           // E[a'](i',j')

                           for (unsigned iP = 0; iP < E.Basis1().size(); ++iP)
                           {
                              MatrixOperator::const_inner_iterator Ai = iterate_at(A.base()[sP].data(), iP, i);
                              if (!Ai)
                                 continue;
                              MatrixOperator::const_inner_iterator Ei = iterate_at(E[aP].data(), iP, jP);
                              if (!Ei)
                                 continue;

                              double degree_iP = degree(A.base().Basis1()[iP]);

                              // now assemble the component
                              double Coeff = tensor_coefficient(B.Basis2()[j],
                                                                B.LocalBasis()[s],
                                                                B.Basis1()[jP],

                                                                M.base().Basis2()[a],
                                                                k->TransformsAs(),
                                                                M.base().Basis1()[aP],

                                                                A.base().Basis2()[i],
                                                                A.base().LocalBasis()[sP],
                                                                A.base().Basis1()[iP]);

                              if (LinearAlgebra::norm_frob(Coeff) > 1E-14)
                              {
                                 int JIndex = JList.Lookup(*Ai, *Ei);
                                 KMat.push_back(std::make_pair(Coeff * (degree_iP / degree_i) * herm(*S), JIndex));
                              }
                           }
                        }
                     }
                  }
                  if (!KMat.empty())
                     C.back().Components.push_back(ElementRec(s, jP, KMat));
               }
            }
         }
      }
   }

   // do the evaluations
   JList.Evaluate();

   for (unsigned n = 0; n < C.size(); ++n)
   {
      if (!C[n].empty())
      {
         set_element(Result[C[n].a].data(), C[n].i, C[n].j, C[n].Evaluate(B,JList));
      }
   }

   Result.debug_check_structure();
   return Result;
}

#if 0
// we probably want a _shift version to handle an implicit qshift
StateComponent
contract_from_left(OperatorComponent const& M,
                   HermitianProxy<StateComponent> const& A,
                   StateComponent const& E,
                   StateComponent const& B)
{
   QuantumNumbers::QuantumNumber I(M.GetSymmetryList());
   return contract_from_left_shift(M, A, I, E, B, I);
}

StateComponent
contract_from_left_shift(OperatorComponent const& M,
                         HermitianProxy<StateComponent> const& A, QuantumNumber const& QShiftA,
                         StateComponent const& E,
                         StateComponent const& B, QuantumNumber const& QShiftB)
#endif

StateComponent
contract_from_left(OperatorComponent const& M,
                   HermitianProxy<StateComponent> const& A,
                   StateComponent const& E,
                   StateComponent const& B)
{
   PRECONDITION_EQUAL(M.LocalBasis1(), A.base().LocalBasis());
   PRECONDITION_EQUAL(M.LocalBasis2(), B.LocalBasis());
   PRECONDITION_EQUAL(M.Basis1(), E.LocalBasis());
   //   DEBUG_PRECONDITION_EQUAL(qshift(A.base().Basis1(), QShiftA), E.Basis1());
   //   DEBUG_PRECONDITION_EQUAL(E.Basis2(), delta_shift(B.Basis1(), QShiftB));
   PRECONDITION_EQUAL(A.base().Basis1(), E.Basis1());
   PRECONDITION_EQUAL(E.Basis2(), B.Basis1());

   // firstly, assemble the required H matrix descriptors
   JMatrixRefList JList;

   // this describes how to construct the final matrices
   std::vector<OuterIndex> C;

   // Firstly, iterate over the outer indices - we want to be careful to only
   // add indices that are definitely needed, the MPS might have more sparseness than
   // the minimum imposed by the symmetries.
   StateComponent Result(M.Basis2(), A.base().Basis2(), B.Basis2());

   // iterate over all possible output matrix elements E'[a](i, j)
   for (unsigned a = 0; a < M.Basis2().size(); ++a)
   {
      for (unsigned i = 0; i < A.base().Basis2().size(); ++i)
      {
         double degree_i = degree(A.base().Basis2()[i]);

         for (unsigned j = 0; j < B.Basis2().size(); ++j)
         {
            if (!is_transform_target(B.Basis2()[j], M.Basis2()[a], A.base().Basis2()[i]))
               continue;

            // Now we have enough to construct the output descriptor
            C.push_back(OuterIndex(i, a, j));

            // Iterate over B[s](jP,j)
            for (unsigned s = 0; s < B.LocalBasis().size(); ++s)
            {
               // We already know j, unlike the F-matrix case we don't know the leading index
               // so we need to iterate
               for (unsigned jP = 0; jP < B.Basis1().size(); ++jP)
               {
                  MatrixOperator::const_inner_iterator J = iterate_at(B[s].data(), jP, j);
                  if (!J)
                     continue;

                  KMatrixDescriptor KMat;

                  // Iterate over the components in M[a'a]
                  // Unlike the F-matrix case, we don't know the leading index so we need to iterate
                  // over both indices
                  for (unsigned aP = 0; aP < M.Basis1().size(); ++aP)
                  {
                     OperatorComponent::const_inner_iterator AA = iterate_at(M.data(), aP, a);
                     if (!AA)
                        continue;

                     // Iterate over the irreducible components of M(aP,a)
                     for (SimpleRedOperator::const_iterator k = AA->begin(); k != AA->end(); ++k)
                     {
                        // *k is an irreducible operator.  Iterate over the components of this operator.
                        // Unlike the F-matrix case, we don't know the leading index so we need to iterate
                        for (unsigned sP = 0; sP < M.LocalBasis1().size(); ++sP)
                        {
                           SimpleOperator::const_inner_iterator S = iterate_at(k->data(), sP, s);
                           if (!S)
                              continue;

                           // The final index is i' - we only need this if the
                           // element exists in both A.base()[s'][i', i] and
                           // E[a'](i',j')

                           for (unsigned iP = 0; iP < E.Basis1().size(); ++iP)
                           {
                              MatrixOperator::const_inner_iterator Ai = iterate_at(A.base()[sP].data(), iP, i);
                              if (!Ai)
                                 continue;
                              MatrixOperator::const_inner_iterator Ei = iterate_at(E[aP].data(), iP, jP);
                              if (!Ei)
                                 continue;

                              double degree_iP = degree(A.base().Basis1()[iP]);

                              // now assemble the component
                              double Coeff = tensor_coefficient(B.Basis2()[j],
                                                                B.LocalBasis()[s],
                                                                B.Basis1()[jP],

                                                                M.Basis2()[a],
                                                                k->TransformsAs(),
                                                                M.Basis1()[aP],

                                                                A.base().Basis2()[i],
                                                                A.base().LocalBasis()[sP],
                                                                A.base().Basis1()[iP]);

                              if (LinearAlgebra::norm_frob(Coeff) > 1E-14)
                              {
                                 int JIndex = JList.Lookup(*Ai, *Ei);
                                 KMat.push_back(std::make_pair(Coeff * (degree_iP / degree_i) * (*S), JIndex));
                              }
                           }
                        }
                     }
                  }
                  if (!KMat.empty())
                     C.back().Components.push_back(ElementRec(s, jP, KMat));
               }
            }
         }
      }
   }

   // do the evaluations
   JList.Evaluate();

   int Sz = C.size();
#pragma omp parallel for schedule(dynamic) num_threads(omp::threads_to_use(Sz))
   for (int n = 0; n < Sz; ++n)
   {
      if (!C[n].empty())
      {
         set_element_lock(Result[C[n].a].data(), C[n].i, C[n].j, C[n].Evaluate(B,JList));
      }
   }

   Result.debug_check_structure();
   return Result;
}
