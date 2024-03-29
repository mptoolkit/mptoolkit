// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/matrixproduct-obsolete/linearwavefunction.cpp
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
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

void LinearWavefunction::normalize()
{
   pvalue_ptr<MPStateComponent> p = Data.front().lock();
   Data.front() = new MPStateComponent((1.0 / norm_frob(*p)) * (*p));
}

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
   return out << psi.SList << psi.Data << psi.Attr;
}

PStream::ipstream& operator>>(PStream::ipstream& in, LinearWavefunction& psi)
{
   return in >> psi.SList >> psi.Data >> psi.Attr;
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
   Result = inject_left_old_interface(M, Result);
   M = M * herm(CollapseBasis(M.Basis2()));
   Result = inject_right_old_interface(Result, M);
   LinearWavefunction::iterator It = Result.begin();
   *It = prod(M, *It);
   return Result;
}

void Conjugate(MPStateComponent& x)
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

MatrixOperator transfer_from_left(MatrixOperator const& m, LinearWavefunction const& Psi)
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

MatrixOperator transfer_from_right(MatrixOperator const& m, LinearWavefunction const& Psi)
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

MPMatrix
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
   MPMatrix E = make_vacuum_state(Psi1.GetSymmetryList());
   while (I1 != Psi1.begin())
   {
      --I1; --I2; --Mi;
      E = operator_prod(*Mi, *I1, E, herm(*I2));
   }
   return conj(E);
}

MPMatrix
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
   MPMatrix E = make_vacuum_state(Psi1.GetSymmetryList());
   while (I1 != Psi1.begin())
   {
      --I1; --I2; --Mi;
      E = operator_prod(*Mi, *I1, E, herm(conj(*I2)));
   }
   return conj(E);
}

std::complex<double>
expectation(LinearWavefunction const& Psi1,
            LinearOperator const& M,
            LinearWavefunction const& Psi2)
{
   MPMatrix E = reduced_matrix_element(Psi1, M, Psi2);
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
   MPMatrix E = reduced_matrix_element_conj(Psi1, M, Psi2);
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

LinearWavefunction inject_left_old_interface(MatrixOperator& M,
                                             LinearWavefunction const& Psi)
{
   LinearWavefunction Result(Psi.GetSymmetryList());
   LinearWavefunction::const_iterator Pi = Psi.begin();
   while (Pi != Psi.end())
   {
      MPStateComponent x = prod(M, *Pi);
      M = TruncateBasis2(x);
      Result.push_back(x);
      ++Pi;
   }
   return Result;
}

LinearWavefunction inject_right_old_interface( LinearWavefunction const& Psi,
                                               MatrixOperator& M)
{
   LinearWavefunction Result(Psi.GetSymmetryList());
   LinearWavefunction::const_iterator Pi = Psi.end();
   while (Pi != Psi.begin())
   {
      --Pi;
      MPStateComponent x = prod(*Pi, M);
      M = TruncateBasis1(x);
      Result.push_front(x);
   }
   return Result;
}

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
   Result = inject_left_old_interface(M, Result);
   M = M * herm(CollapseBasis(M.Basis2()));
   Result = inject_right_old_interface(Result, M);
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
   Result = inject_left_old_interface(M, Result);
   M = M * herm(CollapseBasis(M.Basis2()));
   Result = inject_right_old_interface(Result, M);
   LinearWavefunction::iterator It = Result.begin();
   *It = prod(M, *It);
   return Result;
}

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

   Psi = inject_right_old_interface(Psi, M);
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
   x = inject_left_old_interface(C, x);
   C = C * herm(CollapseBasis(C.Basis2()));
   x = inject_right_old_interface(x, C);
   //   C = RemoveEmptyRows(C);
   LinearWavefunction::iterator It = x.begin();
   *It = prod(C, *It);
   //   x = multiply_left(C, x);
}
