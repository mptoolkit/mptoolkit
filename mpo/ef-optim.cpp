// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mpo/ef-optim.cpp
//
// Copyright (C) 2015-2016 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER
//
// Optimized version of the E-matrix operator_prod
//

#include "mpo/operator_component.h"

#include "common/environment.h"
#include <fstream>

// note: we could parallelize the construction of the G and H indices over jP

typedef std::complex<double> NumberType;
typedef LinearAlgebra::Matrix<NumberType> MatrixType;

#if 0

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
   MatrixOperator::const_inner_iterator x = iterate_at(B[Components[0].s].data(), Components[0].jP, j);
   DEBUG_CHECK(x);

   MatrixType Result = EvaluateK(Components[0].K, J) * (*x);

   for (unsigned n = 1; n < Components.size(); ++n)
   {
      x = iterate_at(B[Components[n].s].data(), Components[n].jP, j);
      DEBUG_CHECK(x);
      Result += EvaluateK(Components[n].K, J) * (*x);
   }
   return Result;
}

#endif

bool DoLogElement = false;

std::map<std::tuple<int,int,int>, int> ResultMap, AMap, BMap, CMap;
std::vector<std::pair<int,int>> ResultSize, ASize, BSize, CSize;
std::map<std::tuple<int,int,int,int>, std::complex<double>> F;

void StartLogElement()
{
   ResultMap.clear();
   ResultSize.clear();
   AMap.clear();
   ASize.clear();
   BMap.clear();
   BSize.clear();
   CMap.clear();
   CSize.clear();
   F.clear();
}

int AddToMap(std::map<std::tuple<int,int,int>, int>& Map, std::tuple<int,int,int> x, std::vector<std::pair<int,int>>& SizeMap, int n, int m)
{
   auto i = Map.find(x);
   if (i == Map.end())
   {
      int r = Map.size();
      Map[x] = r;
      SizeMap.push_back({n,m});
      return r;
   }
   return i->second;
}

void LogElement(std::complex<double> Factor, std::tuple<int,int,int> R, std::tuple<int,int,int> A, std::tuple<int,int,int> B, std::tuple<int,int,int> C, int S1, int S2, int S3, int S4)
{
   int Ri = AddToMap(ResultMap, R, ResultSize, S1, S4);
   int Ai = AddToMap(AMap, A, ASize, S1, S2);
   int Bi = AddToMap(BMap, B, BSize, S2, S3);
   int Ci = AddToMap(CMap, C, CSize, S3, S4);
   F[{Ri,Ai,Bi,Ci}] = Factor;
}

std::ofstream ElementFile(getenv_or_default("MP_ELEMENTFLE", ""), std::ios::trunc | std::ios::out);

void ShowVec(std::ostream& out, std::vector<std::pair<int,int>> const& v)
{
   int i = 0;
   for (auto const& x : v)
   {
      out << (i++) << ' ' << x.first << ' ' << x.second << '\n';
   }
}

void FinishLogElement()
{
   ElementFile << "\nRSize:\n";
   ShowVec(ElementFile, ResultSize);

   ElementFile << "\nASize:\n";
   ShowVec(ElementFile, ASize);

   ElementFile << "\nBSize:\n";
   ShowVec(ElementFile, BSize);

   ElementFile << "\nCSize:\n";
   ShowVec(ElementFile, CSize);

   ElementFile << "\nF:\n";
   for (auto const& x : F)
   {
      ElementFile << std::get<0>(x.first) << ' ' << std::get<1>(x.first) << ' ' << std::get<2>(x.first) << ' ' << std::get<3>(x.first) << ' '
         << x.second.real() << '\n';
   }
   ElementFile << std::endl;
}

// Result[s'](i',i) = M(s',s)[a',a] E[a'](i',j') B[s](j',j) herm(F[a](i,j))
StateComponent
operator_prod_inner(OperatorComponent const& M,
                    StateComponent const& E,
                    StateComponent const& B,
                    LinearAlgebra::HermitianProxy<StateComponent> const& F)
{
   DEBUG_PRECONDITION_EQUAL(M.Basis1(), E.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.Basis2(), F.base().LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.LocalBasis2(), B.LocalBasis());

   DEBUG_PRECONDITION_EQUAL(B.Basis2(), F.base().Basis2());
   DEBUG_PRECONDITION_EQUAL(E.Basis2(), B.Basis1());

   StartLogElement();

   // firstly, assemble the required H matrix descriptors
   //   JMatrixRefList JList;

   // this describes how to construct the final matrices
   //   std::vector<OuterIndex> C;

   // Firstly, iterate over the outer indices - we want to be careful to only
   // add indices that are definitely needed, the MPS might have more sparseness than
   // the minimum imposed by the symmetries.
   StateComponent Result(M.LocalBasis1(), E.Basis1(), F.base().Basis1());

   // iterate over all possible output matrix elements A'[s'](i, j)
   for (unsigned sP = 0; sP < M.LocalBasis1().size(); ++sP)
   {
      for (unsigned iP = 0; iP < E.Basis1().size(); ++iP)
      {
         //double degree_i = degree(E.base().Basis1()[i]);

         for (unsigned i = 0; i < F.base().Basis1().size(); ++i)
         {
            if (!is_transform_target(F.base().Basis1()[i], M.LocalBasis1()[sP], E.Basis1()[iP]))
               continue;

            // Now we have enough to construct the output descriptor
            //      C.push_back(OuterIndex(iP, sP, i));

            // Iterate over F[a](i,j)
            for (unsigned a = 0; a < F.base().LocalBasis().size(); ++a)
            {
               // We already know i, we don't need to iterate
               for (MatrixOperator::MatrixType::data_type::value_type::const_iterator
                       J = iterate(F.base()[a].data().vec()[i]); J; ++J)
               {
                  int j = J.index();

                  //KMatrixDescriptor KMat;

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
                        // We already know sP
                        // We already know the index sP, so iterate over that row
                        for (SimpleOperator::MatrixType::data_type::value_type::const_iterator
                                S = iterate(k->data().vec()[sP]); S; ++S)
                        {
                           int s = S.index();

                           // The final index is j' - we only need this if the
                           // element exists in both E.base()[a'](i',j')
                           //and B.base()[s](j',j)

                           for (unsigned jP = 0; jP < E.Basis2().size(); ++jP)
                           {
                              MatrixOperator::const_inner_iterator EjP = iterate_at(E[aP].data(), iP, jP);
                              if (!EjP)
                                 continue;
                              MatrixOperator::const_inner_iterator BjP = iterate_at(B[s].data(), jP, j);
                              if (!BjP)
                                 continue;

                              //double degree_iP = degree(A.base().Basis1()[iP]);

                              // now assemble the component
                              double Coeff = tensor_coefficient(B.Basis2()[j],
                                                                B.LocalBasis()[s],
                                                                B.Basis1()[jP],

                                                                M.Basis2()[a],
                                                                k->TransformsAs(),
                                                                M.Basis1()[aP],

                                                                F.base().Basis1()[i],
                                                                M.LocalBasis1()[sP],
                                                                E.Basis1()[iP]);

                              if (LinearAlgebra::norm_frob(Coeff) > 1E-14)
                              {
                                 if (!iterate_at(Result[sP].data(), iP, i))
                                 {
                                    set_element(Result[sP].data(), iP, i, Coeff * (*S) * (*EjP) * (*BjP) * herm(*J));
                                 }
                                 else
                                 {
                                    Result[sP](iP, i) += Coeff * (*S) * (*EjP) * (*BjP) * herm(*J);
                                 }
                                 if (DoLogElement)
                                    LogElement(Coeff*(*S), {sP, iP, i}, {aP, iP, jP}, {s, jP, j}, {a, i, j}, (*EjP).size1(), (*EjP).size2(),
                                       (*BjP).size2(), (*J).size1());
                              }
                           }
                        }
                     }
                  }
                  //              if (!KMat.empty())
                  //                 C.back().Components.push_back(ElementRec(s, jP, KMat));
               }
            }
         }
      }
   }
   if (DoLogElement)
      FinishLogElement();

#if 0
   // do the evaluations
   JList.Evaluate();

   for (unsigned n = 0; n < C.size(); ++n)
   {
      if (!C[n].empty())
      {
         set_element(Result[C[n].a].data(), C[n].i, C[n].j, C[n].Evaluate(B,JList));
      }
   }
#endif

   Result.debug_check_structure();
   return Result;
}
