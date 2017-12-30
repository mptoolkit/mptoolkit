// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mpo/f-optim.cpp
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
// Optimized version of the F-matrix operator_prod
//

#include "mpo/operator_component.h"
#include "tensor/tensor_types.h"

// note: we could parallelize the construction of the G and H indices over jP

namespace
{

// This version uses a Matrix* directly
struct GMatrixRef
{
   Matrix const* F;
   Matrix const* BH;

   GMatrixRef(Matrix const& F_, Matrix const& BH_) : F(&F_), BH(&BH_) {}
};

bool
operator<(GMatrixRef const& x, GMatrixRef const& y)
{
   // FIXME: should use std::less here, but in practice it surely doesn't matter
   return x.F < y.F || (x.F == y.F && x.BH < y.BH);
}

class GMatrixRefList
{
   public:
      GMatrixRefList() {}

      // we can reserve a maximum size
      void reserve(int MaxSize);

      // returns the index of a G-matrix (which is inserted if it didn't exist previously)
      int Lookup(Matrix const& F, Matrix const& BH);

      // do the actual evaluations (in parallel)
      void Evaluate();

      // after the evaluation, return the n'th matrix
      Matrix const& operator[](int n) const
      {
         return Data[n].second;
      }

   private:
      typedef std::map<GMatrixRef, int> IndexType;
      IndexType Indices;
      std::vector<std::pair<GMatrixRef, Matrix> > Data;
};

void
GMatrixRefList::reserve(int MaxSize)
{
   Data.reserve(MaxSize);
}

int
GMatrixRefList::Lookup(Matrix const& F, Matrix const& BH)
{
   GMatrixRef Index(F, BH);
   IndexType::const_iterator I = Indices.find(Index);
   if (I == Indices.end())
   {
      int r = Data.size();
      Indices[Index] = r;
      Data.push_back(std::make_pair(Index, Matrix(F.rows(), BH.cols())));
      return r;
   }
   return I->second;
}

void
GMatrixRefList::Evaluate()
{
   for (unsigned n = 0; n < Data.size(); ++n)
   {
      Data[n].second = (*Data[n].first.F) * herm(*Data[n].first.BH);
   }
}

// array of coefficient * index into the GMatrices array
typedef std::vector<std::pair<complex, int> > HMatrixDescriptor;

Matrix EvaluateH(HMatrixDescriptor const& H, GMatrixRefList const& G)
{
   CHECK(!H.empty())("H Matrix descriptor should not be empty!");

   Matrix Result = H[0].first * G[H[0].second];
   for (unsigned n = 1; n < H.size(); ++n)
   {
      Result += H[n].first * G[H[n].second];
   }
   return Result;
}

struct ElementRec
{
   ElementRec(int sP_, int i_, HMatrixDescriptor const& H_)
      : sP(sP_), i(i_), H(H_) {}

   int sP, i;
   HMatrixDescriptor H;
};

struct OuterIndex
{
   OuterIndex(int iP_, int aP_, int jP_)
      : iP(iP_), aP(aP_), jP(jP_) {}

   bool empty() const { return Components.empty(); }

   // do the final evaluation.  Precondition: !empty()
   Matrix Evaluate(StateComponent const& A, GMatrixRefList const& G) const;

   int iP, aP, jP;
   std::vector<ElementRec> Components;
};

Matrix
OuterIndex::Evaluate(StateComponent const& A, GMatrixRefList const& G) const
{
   DEBUG_CHECK(!Components.empty());
   auto x = A[Components[0].sP].row(iP).find(Components[0].i);
   DEBUG_CHECK(x != A[Components[0].sP].row(iP).end());

   Matrix Result = (*x).value * EvaluateH(Components[0].H, G);

   for (unsigned n = 1; n < Components.size(); ++n)
   {
      x = A[Components[n].sP].row(iP).find(Components[n].i);
      DEBUG_CHECK(x != A[Components[n].sP].row(iP).end());
      Result += (*x).value * EvaluateH(Components[n].H, G);
   }
   return Result;
}

} // namespace


StateComponent
contract_from_right_mask(HermitianProxy<OperatorComponent> const& M,
                         StateComponent const& A,
                         StateComponent const& F,
                         HermitianProxy<StateComponent> const& B,
                         std::set<int> const& Mask1,
                         std::set<int> const& Mask2)
{
   PRECONDITION_EQUAL(M.base().LocalBasis1(), A.LocalBasis());
   PRECONDITION_EQUAL(M.base().LocalBasis2(), B.base().LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.base().Basis2(), F.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), F.Basis1());
   DEBUG_PRECONDITION_EQUAL(F.Basis2(), B.base().Basis2());

   using ::norm_frob;

   // firstly, assemble the required H matrix descriptors
   GMatrixRefList G;

   // this describes how to construct the final matrices
   std::vector<OuterIndex> C;

   // Firstly, iterate over the outer indices - we want to be careful to only
   // add indices that are definitely needed, the MPS might have more sparseness than
   // the minimum imposed by the symmetries.
   StateComponent Result(M.base().Basis1(), A.Basis1(), B.base().Basis1());

   // iterate over all possible output matrix elements F'[aP](iP, jP)
   for (unsigned aP = 0; aP < M.base().Basis1().size(); ++aP)
   {
      // skip over masked components
      if (Mask1.find(aP) != Mask1.end())
         continue;

      for (unsigned iP = 0; iP < A.Basis1().size(); ++iP)
      {
         for (unsigned jP = 0; jP < B.base().Basis1().size(); ++jP)
         {
            if (!is_transform_target(B.base().Basis1()[jP], M.base().Basis1()[aP], A.Basis1()[iP]))
               continue;

            // Now we have enough to construct the output descriptor
            C.push_back(OuterIndex(iP, aP, jP));


            // Iterate over A[sP](iP,i) (we don't care about the actual values, just the used elements)
            for (unsigned sP = 0; sP < A.LocalBasis().size(); ++sP)
            {
               // We already know iP, so iterate over the row
               for (auto const& I : A[sP].row(iP))
               {
                  int i = I.col();

                  HMatrixDescriptor HMat;

                  // TODO: we already know AP
                  // Iterate over the components in M.  We already know aP, so we can iterate over
                  // the row of the matrix
                  for (auto const& AA : M.base().row(aP))
                  {
                     int a = AA.col();

                     // skip over masked components
                     if (Mask2.find(a) != Mask2.end())
                        continue;

                     // Iterate over the irreducible components of M(aP,a)
                     for (auto const& k : AA.value)
                     {
                        // k is an irreducible operator.  Iterate over the components of this operator.
                        // We already know the index sP, so iterate over that row
                        for (auto const& S : k.row(sP))
                        {
                           int s = S.col();

                           // The final index is j - we only need this if the
                           // element exists in both B.base()[s][jP, j] and
                           // F[a](i,j)
                           auto BjEnd = B.base()[s].row(jP).end();
                           auto FjEnd = F[a].row(i).end();

                           auto Bj = B.base()[s].row(jP).begin();
                           auto Fj = F[a].row(i).begin();

                           while (Bj != BjEnd && Fj != FjEnd)
                           {
                              if (Bj.col() < Fj.col())
                              {
                                 ++Bj;
                                 continue;
                              }
                              else if (Fj.col() < Bj.col())
                              {
                                 ++Fj;
                                 continue;
                              }

                              int j = Bj.col();
                              DEBUG_CHECK_EQUAL(j, int(Fj.col()));

                              // now assemble the component
                              double Coeff = tensor_coefficient(B.base().Basis2()[j],
                                                                B.base().LocalBasis()[s],
                                                                B.base().Basis1()[jP],

                                                                M.base().Basis2()[a],
                                                                k.TransformsAs(),
                                                                M.base().Basis1()[aP],

                                                                A.Basis2()[i],
                                                                A.LocalBasis()[sP],
                                                                A.Basis1()[iP]);

                              if (norm_frob(Coeff) > 1E-14)
                              {
                                 int GIndex = G.Lookup(*Fj, *Bj);
                                 HMat.push_back(std::make_pair(Coeff * herm(S.value), GIndex));
                              }

                              ++Bj;
                              ++Fj;
                           }
                        }
                     }
                  }
                  if (!HMat.empty())
                     C.back().Components.push_back(ElementRec(sP, i, HMat));
               }
            }
         }
      }
   }

   // do the evaluations
   G.Evaluate();

   for (unsigned n = 0; n < C.size(); ++n)
   {
      if (!C[n].empty())
      {
         Result[C[n].aP].insert(C[n].iP, C[n].jP, C[n].Evaluate(A,G));
      }
   }

   Result.debug_check_structure();
   return Result;
}

StateComponent
contract_from_right(HermitianProxy<OperatorComponent> const& M,
                    StateComponent const& A,
                    StateComponent const& F,
                    HermitianProxy<StateComponent> const& B)
{
   return contract_from_right_mask(M, A, F, B, std::set<int>(), std::set<int>());
}
