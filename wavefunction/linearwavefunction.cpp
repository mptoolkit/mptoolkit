// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// wavefunction/linearwavefunction.cpp
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

#include "linearwavefunction.h"

VectorBasis
LinearWavefunction::Basis1() const
{
   return Data.front().lock()->Basis1();
}

VectorBasis
LinearWavefunction::Basis2() const
{
   return Data.back().lock()->Basis2();
}

QuantumNumbers::QuantumNumber
LinearWavefunction::TransformsAs() const
{
   CHECK(this->is_irreducible());
   return this->Basis1()[0];
}

bool
LinearWavefunction::is_irreducible() const
{
   return this->Basis1().size() == 1;
}

#if 0
void LinearWavefunction::normalize()
{
   pvalue_ptr<StateComponent> p = Data.front().lock();
   Data.front() = new StateComponent((1.0 / norm_frob(*p)) * (*p));
}
#endif

LinearWavefunction::value_type
LinearWavefunction::get_front() const
{
   return *Data.front().lock();
}

LinearWavefunction::value_type
LinearWavefunction::get_back() const
{
   return *Data.back().lock();
}

void LinearWavefunction::set_front(value_type const& x)
{
   Data.front() = handle_type(new value_type(x));
}

void LinearWavefunction::set_back(value_type const& x)
{
   Data.back() = handle_type(new value_type(x));
}

PStream::opstream& operator<<(PStream::opstream& out, LinearWavefunction const& psi)
{
   return out << psi.SList << psi.Data;
}

PStream::ipstream& operator>>(PStream::ipstream& in, LinearWavefunction& psi)
{
   in >> psi.SList >> psi.Data;
   // at some point in the past, the LinearWavefunction didn't track properly the SList.
   // so fix it up if it isn't defined properly
   if (psi.SList.is_null() && !psi.Data.empty())
      psi.SList = psi.Data.front().lock()->GetSymmetryList();
   return in;
}

LinearWavefunction operator*(double a, LinearWavefunction const& x)
{
   LinearWavefunction Result(x);
   Result *= a;
   return Result;
}

LinearWavefunction operator*(LinearWavefunction const& x, double a)
{
   LinearWavefunction Result(x);
   Result *= a;
   return Result;
}

LinearWavefunction operator*(std::complex<double> a, LinearWavefunction const& x)
{
   LinearWavefunction Result(x);
   Result *= a;
   return Result;
}

LinearWavefunction operator*(LinearWavefunction const& x, std::complex<double> a)
{
   LinearWavefunction Result(x);
   Result *= a;
   return Result;
}

LinearWavefunction& operator*=(LinearWavefunction& psi, double a)
{
   LinearWavefunction::iterator I = psi.begin();
   MatrixOperator c = MatrixOperator::make_identity(I->Basis1());
   c *= a;
   *I = prod(c, *I);
   return psi;
}

LinearWavefunction& operator*=(LinearWavefunction& psi, std::complex<double> a)
{
   LinearWavefunction::iterator I = psi.begin();
   MatrixOperator c = MatrixOperator::make_identity(I->Basis1());
   c *= a;
   *I = prod(c, *I);
   return psi;
}

MatrixOperator CollapseBasis(VectorBasis const& b)
{
   std::set<QuantumNumber> QN(b.Basis().begin(), b.Basis().end());
   BasisList NewB(QN.begin(), QN.end());
   VectorBasis vb(NewB);
   MatrixOperator C(vb, b);
   for (unsigned j = 0; j < b.size(); ++j)
   {
      unsigned i = std::find(NewB.begin(), NewB.end(), b[j]) - NewB.begin();
      C(i,j) = LinearAlgebra::Matrix<double>(1, b.dim(j), 1.0);
   }
   return C;
}

LinearWavefunction operator+(LinearWavefunction const& x, LinearWavefunction const& y)
{
   // quick return for sum involving a null wavefunction
   if (x.empty())
      return y;
   if (y.empty())
      return x;

   CHECK_EQUAL(x.GetSymmetryList(), y.GetSymmetryList());
   LinearWavefunction Result(x.GetSymmetryList());
   LinearWavefunction::const_iterator xi = x.begin(), yi = y.begin();

   SumBasis<VectorBasis> B1(xi->Basis1(), yi->Basis1());
   while (xi != x.end())
   {
      SumBasis<VectorBasis> B2(xi->Basis2(), yi->Basis2());
      Result.push_back(tensor_sum(*xi, *yi, B1, B2));
      B1 = B2;
      ++xi;
      ++yi;
   }

   MatrixOperator M = CollapseBasis(Result.Basis1());
   M = left_orthogonalize(M, Result);
   M = M * herm(CollapseBasis(M.Basis2()));
   M = right_orthogonalize(Result, M);
   LinearWavefunction::iterator It = Result.begin();
   *It = prod(M, *It);
   return Result;
}

void Conjugate(StateComponent& x)
{
   for (std::size_t i = 0; i < x.size(); ++i)
   {
      x[i] = conj(x[i]);
   }
}

void Conjugate(LinearWavefunction& Psi)
{
   for (LinearWavefunction::iterator I = Psi.begin(); I != Psi.end(); ++I)
   {
      Conjugate(*I);
   }
}

LinearWavefunction conj(LinearWavefunction const& x)
{
   LinearWavefunction Result = x;
   Conjugate(Result);
   return Result;
}

LinearWavefunction operator-(LinearWavefunction const& x, LinearWavefunction const& y)
{
   LinearWavefunction z = y;
   z *= -1.0;
   return x+z;
}

MatrixOperator
inject_left(MatrixOperator const& m, LinearWavefunction const& Psi)
{
   MatrixOperator Result = m;
   LinearWavefunction::const_iterator I = Psi.begin();
   while (I != Psi.end())
   {
      Result = operator_prod(herm(*I), Result, *I);
      ++I;
   }
   return Result;
}

MatrixOperator
inject_left(MatrixOperator const& m,
            LinearWavefunction const& Psi1,
            LinearWavefunction const& Psi2)
{
   CHECK_EQUAL(Psi1.size(), Psi2.size());
   MatrixOperator Result = m;
   LinearWavefunction::const_iterator I1 = Psi1.begin();
   LinearWavefunction::const_iterator I2 = Psi2.begin();
   while (I1 != Psi1.end())
   {
      Result = operator_prod(herm(*I1), Result, *I2);
      ++I1; ++I2;
   }
   return Result;
}

MatrixOperator
inject_right(MatrixOperator const& m, LinearWavefunction const& Psi)
{
   MatrixOperator Result = m;
   LinearWavefunction::const_iterator I = Psi.end();
   while (I != Psi.begin())
   {
      --I;
      Result = operator_prod(*I, Result, herm(*I));
   }
   return Result;
}

MatrixOperator
inject_right(MatrixOperator const& m,
            LinearWavefunction const& Psi1,
            LinearWavefunction const& Psi2)
{
   CHECK_EQUAL(Psi1.size(), Psi2.size());
   MatrixOperator Result = m;
   LinearWavefunction::const_iterator I1 = Psi1.end();
   LinearWavefunction::const_iterator I2 = Psi2.end();
   while (I1 != Psi1.begin())
   {
      --I1; --I2;
      Result = operator_prod(*I1, Result, herm(*I2));
   }
   return Result;
}

MatrixOperator
operator_prod(HermitianProxy<LinearWavefunction> const& A,
              MatrixOperator const& m,
              LinearWavefunction const& B)
{
   PRECONDITION_EQUAL(A.base().size(), B.size());
   DEBUG_PRECONDITION_EQUAL(m.Basis1(), A.base().Basis1());
   DEBUG_PRECONDITION_EQUAL(m.Basis2(), B.Basis1());
   MatrixOperator Result = m;
   LinearWavefunction::const_iterator AI = A.base().begin(), BI = B.begin();
   while (BI != B.end())
   {
      Result = operator_prod(herm(*AI), Result, *BI);
      ++AI; ++BI;
   }
   return Result;
}

MatrixOperator
operator_prod(LinearWavefunction const& A,
              MatrixOperator const& m,
              HermitianProxy<LinearWavefunction> const& B)
{
   PRECONDITION_EQUAL(A.size(), B.base().size());
   DEBUG_PRECONDITION_EQUAL(m.Basis1(), A.Basis2());
   DEBUG_PRECONDITION_EQUAL(m.Basis2(), B.base().Basis2());
   MatrixOperator Result = m;
   LinearWavefunction::const_iterator AI = A.end(), BI = B.base().end();
   while (AI != A.begin())
   {
      --AI; --BI;
      Result = operator_prod(*AI, Result, herm(*BI));
   }
   return Result;
}

std::complex<double>
overlap(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2)
{
   CHECK_EQUAL(Psi1.GetSymmetryList(), Psi2.GetSymmetryList());
   CHECK_EQUAL(Psi1.size(), Psi2.size());

   LinearWavefunction::const_iterator I1 = Psi1.end(), I2 = Psi2.end();
   MatrixOperator E = make_vacuum_matrix(Psi1.GetSymmetryList());
   while (I1 != Psi1.begin())
   {
      --I1; --I2;
      E = operator_prod(*I2, E, herm(*I1));
   }
   return trace(E);
}

std::complex<double>
overlap_conj(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2)
{
   CHECK_EQUAL(Psi1.GetSymmetryList(), Psi2.GetSymmetryList());
   CHECK_EQUAL(Psi1.size(), Psi2.size());

   LinearWavefunction::const_iterator I1 = Psi1.end(), I2 = Psi2.end();
   MatrixOperator E = make_vacuum_matrix(Psi1.GetSymmetryList());
   while (I1 != Psi1.begin())
   {
      --I1; --I2;
      E = operator_prod(conj(*I2), E, herm(*I1));
   }
   return trace(E);
}

#if 0
StateComponent
reduced_matrix_element(LinearWavefunction const& Psi1,
                       LinearOperator const& M,
                       LinearWavefunction const& Psi2)
{
   CHECK_EQUAL(Psi1.GetSymmetryList(), Psi2.GetSymmetryList());
   CHECK_EQUAL(Psi1.GetSymmetryList(), M.GetSymmetryList());
   CHECK_EQUAL(Psi1.size(), Psi2.size());
   CHECK_EQUAL(int(Psi1.size()), int(M.size()));
   LinearWavefunction::const_iterator I1 = Psi1.end(), I2 = Psi2.end();
   MPOperator::const_iterator Mi = M.end();
   StateComponent E = make_vacuum_state(Psi1.GetSymmetryList());
   while (I1 != Psi1.begin())
   {
      --I1; --I2; --Mi;
      E = operator_prod(*Mi, *I1, E, herm(*I2));
   }
   return E;
}

StateComponent
reduced_matrix_element_conj(LinearWavefunction const& Psi1,
                            LinearOperator const& M,
                            LinearWavefunction const& Psi2)
{
   CHECK_EQUAL(Psi1.GetSymmetryList(), Psi2.GetSymmetryList());
   CHECK_EQUAL(Psi1.GetSymmetryList(), M.GetSymmetryList());
   CHECK_EQUAL(Psi1.size(), Psi2.size());
   CHECK_EQUAL(int(Psi1.size()), int(M.size()));
   LinearWavefunction::const_iterator I1 = Psi1.end(), I2 = Psi2.end();
   MPOperator::const_iterator Mi = M.end();
   StateComponent E = make_vacuum_state(Psi1.GetSymmetryList());
   while (I1 != Psi1.begin())
   {
      --I1; --I2; --Mi;
      E = operator_prod(*Mi, *I1, E, herm(conj(*I2)));
   }
   return E;
}

// These have been updated on trunk/.  Don't use these versions
// as they are here.

std::complex<double>
expectation(LinearWavefunction const& Psi1,
            LinearOperator const& M,
            LinearWavefunction const& Psi2)
{
   StateComponent E = reduced_matrix_element(Psi1, M, Psi2);
   //   CHECK_EQUAL(E.size(), 1)("expectation: not a scalar value");
   //   CHECK(is_scalar(E.SiteBasis()[0]));
   // return the scalar part
   std::complex<double> Result = 0;
   for (unsigned i = 0; i < E.SiteBasis().size(); ++i)
   {
      if (is_scalar(E.SiteBasis()[i]))
         Result += trace(E[i]);
   }
   return Result;
}

std::complex<double>
expectation_conj(LinearWavefunction const& Psi1,
                 LinearOperator const& M,
                 LinearWavefunction const& Psi2)
{
   StateComponent E = reduced_matrix_element_conj(Psi1, M, Psi2);
   //   CHECK_EQUAL(E.size(), 1)("expectation: not a scalar value");
   //   CHECK(is_scalar(E.SiteBasis()[0]));
   // return the scalar part
   std::complex<double> Result = 0;
   for (unsigned i = 0; i < E.SiteBasis().size(); ++i)
   {
      if (is_scalar(E.SiteBasis()[i]))
         Result += trace(E[i]);
   }
   return Result;
}
#endif

MatrixOperator
left_orthogonalize(MatrixOperator const& Mat, LinearWavefunction& Psi, int Verbose)
{
   MatrixOperator M = Mat;
   LinearWavefunction::iterator Pi = Psi.begin();
   int n = 0;
   while (Pi != Psi.end())
   {
      if (Verbose)
         std::cout << "orthogonalizing site " << n << std::endl;
      StateComponent x = prod(M, *Pi);
      M = TruncateBasis2(x);
      *Pi = x;
      ++Pi; ++n;
   }
   return M;
}

MatrixOperator
right_orthogonalize(LinearWavefunction& Psi, MatrixOperator const& Mat, int Verbose)
{
   MatrixOperator M = Mat;
   LinearWavefunction Result(Psi.GetSymmetryList());
   LinearWavefunction::iterator Pi = Psi.end();
   int n = Psi.size();
   while (Pi != Psi.begin())
   {
      --Pi; --n;
      if (Verbose)
         std::cout << "orthogonalizing site " << n << std::endl;
      StateComponent x = prod(*Pi, M);
      M = TruncateBasis1(x);
      *Pi = x;
   }
   return M;
}

#if 0
LinearWavefunction prod_helper(MPOperator::const_iterator Mi,
                               LinearWavefunction::const_iterator Pi,
                               MPOperator::const_iterator Mend,
                               LinearWavefunction::const_iterator Pend,
                               ProductBasis<BasisList, VectorBasis> B1)
{
   LinearWavefunction Result(Pi->GetSymmetryList());
   while (Pi != Pend)
   {
      ProductBasis<BasisList, VectorBasis> B2 =
         make_product_basis(Mi->Basis2(), Pi->Basis2());
      Result.push_back(mp_prod(*Mi, *Pi, B1, B2));
      B1 = B2;
      ++Mi;
      ++Pi;
   }
   CHECK(Mi == Mend); // make sure our two iterators are in lock-step
   return Result;
}

LinearWavefunction prod(LinearOperator const& Op, LinearWavefunction const& Psi)
{
   MPOperator::const_iterator Mi = Op.begin();
   LinearWavefunction::const_iterator Pi = Psi.begin();
   ProductBasis<BasisList, VectorBasis> B1 =
      make_product_basis(Mi->Basis1(), Pi->Basis1());
   LinearWavefunction Result = prod_helper(Mi, Pi, Op.end(), Psi.end(), B1);

   // normalization
   MatrixOperator M = MatrixOperator::make_identity(Result.Basis1());
   M = left_orthogonalize(M, Result);
   M = M * herm(CollapseBasis(M.Basis2()));
   M = right_orthogonalize(Result, M);
   LinearWavefunction::iterator It = Result.begin();
   *It = prod(M, *It);
   return Result;
}

LinearWavefunction prod(LinearOperator const& Op, LinearWavefunction const& Psi,
                        QuantumNumbers::QuantumNumber const& q)
{
   MPOperator::const_iterator Mi = Op.begin();
   LinearWavefunction::const_iterator Pi = Psi.begin();
   ProductBasis<BasisList, VectorBasis> B1 =
      make_product_basis(Mi->Basis1(), Pi->Basis1(), q);
   LinearWavefunction Result = prod_helper(Mi, Pi, Op.end(), Psi.end(), B1);

   // normalization
   MatrixOperator M = MatrixOperator::make_identity(Result.Basis1());
   M = left_orthogonalize(M, Result);
   M = M * herm(CollapseBasis(M.Basis2()));
   M = right_orthogonalize(Result, M);
   LinearWavefunction::iterator It = Result.begin();
   *It = prod(M, *It);
   return Result;
}
#endif

double norm_2_sq(LinearWavefunction const& Psi)
{
   return norm_frob_sq(*Psi.begin());
}

void truncate(LinearWavefunction& Psi, StatesInfo const& SInfo, bool ShowStates)
{
   LinearWavefunction::iterator I = Psi.begin();
   MatrixOperator M = MatrixOperator::make_identity(I->Basis1());
   int BondNr = 1;
   while (I != Psi.end())
   {
      *I = prod(M, *I);
      M = ExpandBasis2(*I);
      DensityMatrix<MatrixOperator> DM(scalar_prod(M, herm(M)));
      TruncationInfo Info;
      MatrixOperator U =
         DM.ConstructTruncator(DM.begin(),
                               TruncateFixTruncationErrorAbsolute(DM.begin(),
                                                                  DM.end(),
                                                                  SInfo,
                                                                  Info));
      if (ShowStates)
         std::cerr << "bond=" << BondNr
                   << ", states=" << Info.KeptStates()
                   << ", trunc=" << Info.TruncationError()
                   << ", largest_discarded_evalue=" << Info.LargestDiscardedEigenvalue()
                   << '\n';
      *I = prod(*I, herm(U));
      M = U*M;
      ++I;
      ++BondNr;
   }

   if (ShowStates)
      std::cerr << "Orthogonalizing...\n";

   M = right_orthogonalize(Psi, M);
   I = Psi.begin();
   *I = prod(M, *I);

   if (ShowStates)
      std::cerr << "Done.\n";
}

void project(LinearWavefunction& x, QuantumNumbers::QuantumNumber const& q)
{
   MatrixOperator U = CollapseBasis(x.Basis1());
   unsigned i = std::distance(U.Basis1().Basis().begin(),
                              std::find(U.Basis1().Basis().begin(),
                                        U.Basis1().Basis().end(), q));

   if (i == U.Basis1().size())  // if q does not exist in the basis the result is zero
   {
      x = LinearWavefunction(x.GetSymmetryList());
      return;
   }

   VectorBasis FinalB(x.GetSymmetryList());
   FinalB.push_back(q, 1);
   MatrixOperator Up(FinalB, U.Basis1());
   Up(0, i) = LinearAlgebra::Matrix<std::complex<double> >(1,1,1.0);

   MatrixOperator C = Up*U;
   C = left_orthogonalize(C, x);
   C = C * herm(CollapseBasis(C.Basis2()));
   C = right_orthogonalize(x, C);
   //   C = RemoveEmptyRows(C);
   LinearWavefunction::iterator It = x.begin();
   *It = prod(C, *It);
   //   x = multiply_left(C, x);
}

std::vector<SimpleOperator>
make_identity_string_operator(std::vector<BasisList> const& Basis)
{
   std::vector<SimpleOperator> Result;
   Result.reserve(Basis.size());
   for (unsigned i = 0; i < Basis.size(); ++i)
   {
      Result.push_back(SimpleOperator::make_identity(Basis[i]));
   }
   return Result;
}

std::vector<BasisList>
ExtractLocalBasis(LinearWavefunction const& Psi)
{
   std::vector<BasisList> Basis(Psi.size());
   LinearWavefunction::const_iterator Li = Psi.begin();
   for (unsigned i = 0; i < Psi.size(); ++i, ++Li)
   {
      Basis[i] = Li->LocalBasis();
   }
   return Basis;
}

LinearWavefunction coarse_grain(LinearWavefunction const& x, int N)
{
   CHECK(x.size() % N == 0)(x.size())(N);
   LinearWavefunction Result(x.GetSymmetryList());
   LinearWavefunction::const_iterator I = x.begin();
   while (I != x.end())
   {
      StateComponent A = *I;
      for (int i = 1; i < N; ++i)
      {
	 ++I;
	 A = local_tensor_prod(A, *I);
      }
      ++I;
      Result.push_back(A);
   }
   return Result;
}
