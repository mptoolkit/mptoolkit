// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mpo/e-optim.cpp
//
// Copyright (C) 2015-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
#include "tensor/tensor_types.h"

// note: we could parallelize the construction of the G and H indices over jP

namespace
{

// This version uses a Matrix* directly
struct JMatrixRef
{
   Matrix const* AH;
   Matrix const* E;

   JMatrixRef(Matrix const& AH_, Matrix const& E_) : AH(&AH_), E(&E_) {}
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
      int Lookup(Matrix const& AH, Matrix const& E);

      // do the actual evaluations (in parallel)
      void Evaluate();

      // after the evaluation, return the n'th matrix
      Matrix const& operator[](int n) const
      {
         return Data[n].second;
      }

   private:
      typedef std::map<JMatrixRef, int> IndexType;
      IndexType Indices;
      std::vector<std::pair<JMatrixRef, Matrix> > Data;
};

void
JMatrixRefList::reserve(int MaxSize)
{
   Data.reserve(MaxSize);
}

int
JMatrixRefList::Lookup(Matrix const& AH, Matrix const& E)
{
   JMatrixRef Index(AH, E);
   IndexType::const_iterator I = Indices.find(Index);
   if (I == Indices.end())
   {
      int r = Data.size();
      Indices[Index] = r;
      TRACE(AH.cols());
      TRACE(E.cols());
      Data.push_back(std::make_pair(Index, Matrix(AH.cols(), E.cols())));
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
typedef std::vector<std::pair<complex, int> > KMatrixDescriptor;

Matrix EvaluateK(KMatrixDescriptor const& K, JMatrixRefList const& J)
{
   CHECK(!K.empty())("K Matrix descriptor should not be empty!");

   Matrix Result = K[0].first * J[K[0].second];
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
   Matrix Evaluate(StateComponent const& A, JMatrixRefList const& J) const;

   int i, a, j;
   std::vector<ElementRec> Components;
};

Matrix
OuterIndex::Evaluate(StateComponent const& B, JMatrixRefList const& J) const
{
   DEBUG_CHECK(!Components.empty());
   auto x = B[Components[0].s].row(Components[0].jP).find(j);
   DEBUG_CHECK(x != B[Components[0].s].row(Components[0].jP).end());

   Matrix Result = EvaluateK(Components[0].K, J) * (*x).value;

   for (unsigned n = 1; n < Components.size(); ++n)
   {
      x = B[Components[n].s].row(Components[n].jP).find(j);
      DEBUG_CHECK(x != B[Components[n].s].row(Components[n].jP).end());
      Result += EvaluateK(Components[n].K, J) * (*x).value;
   }
   return Result;
}

} // namespace

StateComponent
contract_from_left(OperatorComponent const& M,
                   HermitianProxy<StateComponent> const& A,
                   StateComponent const& E,
                   StateComponent const& B)
{
   PRECONDITION_EQUAL(M.LocalBasis1(), A.base().LocalBasis());
   PRECONDITION_EQUAL(M.LocalBasis2(), B.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.Basis1(), E.LocalBasis());
   //   DEBUG_PRECONDITION_EQUAL(qshift(A.base().Basis1(), QShiftA), E.Basis1());
   //   DEBUG_PRECONDITION_EQUAL(E.Basis2(), delta_shift(B.Basis1(), QShiftB));
   DEBUG_PRECONDITION_EQUAL(A.base().Basis1(), E.Basis1());
   DEBUG_PRECONDITION_EQUAL(E.Basis2(), B.Basis1());

   using ::norm_frob;

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
                  auto J = B[s].row(jP).find(j);
                  if (J == B[s].row(jP).end())
                     continue;

                  KMatrixDescriptor KMat;

                  // Iterate over the components in M[a'a]
                  // Unlike the F-matrix case, we don't know the leading index so we need to iterate
                  // over both indices
                  for (unsigned aP = 0; aP < M.Basis1().size(); ++aP)
                  {
                     auto AA = M.row(aP).find(a);
                     if (AA == M.row(aP).end())
                        continue;

                     // Iterate over the irreducible components of M(aP,a)
                     for (auto const& k : (*AA).value)
                     {
                        // k is an irreducible operator.  Iterate over the components of this operator.
                        // Unlike the F-matrix case, we don't know the leading index so we need to iterate
                        for (unsigned sP = 0; sP < M.LocalBasis1().size(); ++sP)
                        {
                           auto S = k.row(sP).find(s);
                           if (S == k.row(sP).end())
                              continue;

                           // The final index is i' - we only need this if the
                           // element exists in both A.base()[s'][i', i] and
                           // E[a'](i',j')

                           for (unsigned iP = 0; iP < E.Basis1().size(); ++iP)
                           {
                              auto Ai = A.base()[sP].row(iP).find(i);
                              if (Ai == A.base()[sP].row(iP).end())
                                 continue;
                              auto Ei = E[aP].row(iP).find(jP);
                              if (Ei == E[aP].row(iP).end())
                                 continue;

                              TRACE(aP)(iP)(jP);
                              TRACE(&((*Ei).value));

                              double degree_iP = degree(A.base().Basis1()[iP]);

                              // now assemble the component
                              double Coeff = tensor_coefficient(B.Basis2()[j],
                                                                B.LocalBasis()[s],
                                                                B.Basis1()[jP],

                                                                M.Basis2()[a],
                                                                k.TransformsAs(),
                                                                M.Basis1()[aP],

                                                                A.base().Basis2()[i],
                                                                A.base().LocalBasis()[sP],
                                                                A.base().Basis1()[iP]);

                              if (norm_frob(Coeff) > 1E-14)
                              {
                                 int JIndex = JList.Lookup(*Ai, *Ei);
                                 KMat.push_back(std::make_pair(Coeff * (degree_iP / degree_i) * (*S).value, JIndex));
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
         Result[C[n].a].insert(C[n].i, C[n].j, C[n].Evaluate(B,JList));
      }
   }

   Result.debug_check_structure();
   return Result;
}
