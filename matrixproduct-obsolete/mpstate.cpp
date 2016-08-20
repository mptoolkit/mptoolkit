// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/mpstate.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "mpstate.h"
#include "quantumnumbers/quantumnumber.h"
#include "quantumnumbers/u1.h"
#include "tensor/tensorproduct.h"
#include <tuple>
#include "linearalgebra/matrix_utility.h"

MPStateComponent make_vacuum_state(QuantumNumbers::SymmetryList const& S)
{
   BasisList Vacuum = make_vacuum_basis(S);
   MPStateComponent Result(Vacuum, VectorBasis(Vacuum), VectorBasis(Vacuum));
   Result[0](0,0) = LinearAlgebra::Matrix<double>(1,1,1);
   return Result;
}

MPStateComponent make_vacuum_state(QuantumNumbers::QuantumNumber const& Q)
{
   BasisList Vacuum = make_vacuum_basis(Q.GetSymmetryList());
   BasisList vq(Q.GetSymmetryList());
   vq.push_back(Q);
   VectorBasis V(vq);
   MPStateComponent Result(Vacuum, V, V);
   Result[0](0,0) = LinearAlgebra::Matrix<double>(1,1,1);
   return Result;
}

MPStateComponent&
operator+=(MPStateComponent& x, MPStateComponent const& y)
{
   CHECK_EQUAL(x.SiteBasis(), y.SiteBasis());
   CHECK_EQUAL(x.Basis1(), y.Basis1());
   CHECK_EQUAL(x.Basis2(), y.Basis2());
   for (unsigned i = 0; i < x.size(); ++i)
   {
      x[i] += y[i];
   }
   return x;
}

MPStateComponent
operator+(MPStateComponent const& x, MPStateComponent const& y)
{
   CHECK_EQUAL(x.SiteBasis(), y.SiteBasis());
   CHECK_EQUAL(x.Basis1(), y.Basis1());
   CHECK_EQUAL(x.Basis2(), y.Basis2());
   MPStateComponent Result(x);
   for (unsigned i = 0; i < Result.size(); ++i)
   {
      Result[i] += y[i];
   }
   return Result;
}

MPStateComponent&
operator-=(MPStateComponent& x, MPStateComponent const& y)
{
   CHECK_EQUAL(x.SiteBasis(), y.SiteBasis());
   CHECK_EQUAL(x.Basis1(), y.Basis1());
   CHECK_EQUAL(x.Basis2(), y.Basis2());
   for (unsigned i = 0; i < x.size(); ++i)
   {
      x[i] -= y[i];
   }
   return x;
}

MPStateComponent
operator-(MPStateComponent const& x, MPStateComponent const& y)
{
   CHECK_EQUAL(x.SiteBasis(), y.SiteBasis());
   CHECK_EQUAL(x.Basis1(), y.Basis1());
   CHECK_EQUAL(x.Basis2(), y.Basis2());
   MPStateComponent Result(x);
   for (unsigned i = 0; i < Result.size(); ++i)
   {
      Result[i] -= y[i];
   }
   return Result;
}

namespace LinearAlgebra
{

MatrixOperator
ScalarProd<MPStateComponent, HermitianProxy<MPStateComponent> >::
operator()(MPStateComponent const& A, HermitianProxy<MPStateComponent> const& B) const
{
   DEBUG_PRECONDITION_EQUAL(A.SiteBasis(), B.base().SiteBasis());
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), B.base().Basis2());

   QuantumNumbers::QuantumNumber Ident(A.GetSymmetryList());
   MatrixOperator Result(A.Basis1(), B.base().Basis1(), Ident);

   for (std::size_t s = 0; s < A.SiteBasis().size(); ++s)
   {
      Result += scalar_prod(A[s], herm(B.base()[s]));
   }
   return Result;
}

MatrixOperator
ScalarProd<HermitianProxy<MPStateComponent>, MPStateComponent>::
operator()(HermitianProxy<MPStateComponent> const& A, MPStateComponent const& B) const
{
   DEBUG_PRECONDITION_EQUAL(A.base().SiteBasis(), B.SiteBasis());
   DEBUG_PRECONDITION_EQUAL(A.base().Basis1(), B.Basis1());

   QuantumNumbers::QuantumNumber Ident(B.GetSymmetryList());
   MatrixOperator Result(A.base().Basis2(), B.Basis2(), Ident);

   for (std::size_t s = 0; s < B.SiteBasis().size(); ++s)
   {
      Result += scalar_prod(herm(A.base()[s]), B[s]);
   }
   return Result;
}

} // namespace LinearAlgebra

SimpleOperator trace_prod(LinearAlgebra::HermitianProxy<MPStateComponent> const& A,
                          MPStateComponent const& B)
{
   QuantumNumbers::QuantumNumber Ident(B.GetSymmetryList());
   SimpleOperator Result(A.base().SiteBasis(), B.SiteBasis(), Ident);

   BasisList const& Sa = Result.Basis1();
   BasisList const& Sb = Result.Basis2();

   for (std::size_t sa = 0; sa < Sa.size(); ++sa)
   {
      for (std::size_t sb = 0; sb < Sb.size(); ++sb)
      {
         //if (is_transform_target(Sb[sb], adjoint(Sa[sa]), Ident))
         if (Sb[sb] == Sa[sa])
            Result(sa, sb) = inner_prod(A.base()[sa], B[sb]);
      }
   }
   return Result;
}

SimpleOperator trace_prod(MPStateComponent const& A,
                          LinearAlgebra::HermitianProxy<MPStateComponent> const& B)
{
   QuantumNumbers::QuantumNumber Ident(A.GetSymmetryList());

   SimpleOperator Result(A.SiteBasis(), B.base().SiteBasis(), Ident);

   BasisList const& Sa = Result.Basis1();
   BasisList const& Sb = Result.Basis2();

   for (std::size_t sa = 0; sa < Sa.size(); ++sa)
   {
      for (std::size_t sb = 0; sb < Sb.size(); ++sb)
      {
         //if (is_transform_target(adjoint(Sb[sb]), Sa[sa], Ident))
         if (Sb[sb] == Sa[sa])
         {
            Result(sa, sb) = inner_prod(A[sa], B.base()[sb]);
         }
      }
   }
   return Result;
}

// does Result' = sum_{s',s} A[s'] * herm(B[s]) * M(s',s)
// M must transform as the identity representation.
// TODO: generalize this beyond identity op
// FIXME: This probably depends on the matrix elements of M having the correct norm convention
MatrixOperator operator_prod(SimpleOperator const& M,
                             MPStateComponent const& A,
                             LinearAlgebra::HermitianProxy<MPStateComponent> const& B)
{
   DEBUG_PRECONDITION_EQUAL(M.Basis1(), A.SiteBasis());
   DEBUG_PRECONDITION_EQUAL(M.Basis2(), B.base().SiteBasis());
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), B.base().Basis2());
   //DEBUG_PRECONDITION_EQUAL(M.TransformsAs(), QuantumNumbers::QuantumNumber(M.GetSymmetryList()));

   QuantumNumbers::QuantumNumber Ident(A.GetSymmetryList());
   MatrixOperator Result(A.Basis1(), B.base().Basis1(), M.TransformsAs());

   for (const_iterator<SimpleOperator>::type I = iterate(M); I; ++I)
   {
      for (const_inner_iterator<SimpleOperator>::type J = iterate(I); J; ++J)
      {
         Result += (*J) * triple_prod(A[J.index1()],
                                      MatrixOperator::make_identity(A.Basis2()),
                                      herm(B.base()[J.index2()]),
                                      Result.TransformsAs(),
                                      Result.TransformsAs());
         //Result += (*J) * scalar_prod(A[J.index1()], herm(B.base()[J.index2()]));
      }
   }
   return Result;
}

MatrixOperator operator_prod(LinearAlgebra::HermitianProxy<SimpleOperator> const& M,
                             LinearAlgebra::HermitianProxy<MPStateComponent> const& A,
                             MPStateComponent const& B)
{
   DEBUG_PRECONDITION_EQUAL(M.base().Basis1(), A.base().SiteBasis());
   DEBUG_PRECONDITION_EQUAL(M.base().Basis2(), B.SiteBasis());
   DEBUG_PRECONDITION_EQUAL(A.base().Basis1(), B.Basis1());
   //DEBUG_PRECONDITION_EQUAL(M.TransformsAs(), QuantumNumbers::QuantumNumber(M.GetSymmetryList()));

   QuantumNumbers::QuantumNumber Ident(B.GetSymmetryList());
   MatrixOperator Result(A.base().Basis2(), B.Basis2(), adjoint(M.base().TransformsAs()));

   for (const_iterator<SimpleOperator>::type I = iterate(M.base()); I; ++I)
   {
      for (const_inner_iterator<SimpleOperator>::type J = iterate(I); J; ++J)
      {
         Result += herm(*J) * triple_prod(herm(A.base()[J.index1()]),
                                          MatrixOperator::make_identity(B.Basis1()),
                                          B[J.index2()],
                                          Result.TransformsAs(),
                                          Result.TransformsAs());
         //Result += (*J) * scalar_prod(herm(A.base()[J.index1()]), B[J.index2()]);
      }
   }
   return Result;
}

std::ostream& operator<<(std::ostream& out, MPStateComponent const& Psi)
{
  out << "Site basis:\n" << Psi.SiteBasis() << "Basis1:\n" << Psi.Basis1()
      << "Basis2:\n" << Psi.Basis2() << '\n' << Psi.Data;
  return out;
}

MatrixOperator operator_prod(SimpleOperator const& M,
                             MPStateComponent const& A,
                             MatrixOperator const& E,
                             LinearAlgebra::HermitianProxy<MPStateComponent> const& B,
                             QuantumNumbers::QuantumNumber const& q)
{
   DEBUG_PRECONDITION_EQUAL(M.Basis1(), A.SiteBasis());
   DEBUG_PRECONDITION_EQUAL(M.Basis2(), B.base().SiteBasis());
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), E.Basis1());
   DEBUG_PRECONDITION_EQUAL(E.Basis2(), B.base().Basis2());

   MatrixOperator Result(A.Basis1(), B.base().Basis1(), q);

   for (const_iterator<SimpleOperator>::type I = iterate(M); I; ++I)
   {
      for (const_inner_iterator<SimpleOperator>::type J = iterate(I); J; ++J)
      {
         Result += (*J) * triple_prod(A[J.index1()], E, herm(B.base()[J.index2()]),
                                      M.TransformsAs(), q);
      }
   }
   return Result;
}

MatrixOperator operator_prod(LinearAlgebra::HermitianProxy<SimpleOperator> const& M,
                             LinearAlgebra::HermitianProxy<MPStateComponent> const& A,
                             MatrixOperator const& E,
                             MPStateComponent const& B,
                             QuantumNumbers::QuantumNumber const& q)
{
   DEBUG_PRECONDITION_EQUAL(M.base().Basis1(), A.base().SiteBasis());
   DEBUG_PRECONDITION_EQUAL(M.base().Basis2(), B.SiteBasis());
   DEBUG_PRECONDITION_EQUAL(A.base().Basis1(), E.Basis1());
   DEBUG_PRECONDITION_EQUAL(E.Basis2(), B.Basis1());

   MatrixOperator Result(A.base().Basis2(), B.Basis2(), q);

   for (const_iterator<SimpleOperator>::type I = iterate(M.base()); I; ++I)
   {
      for (const_inner_iterator<SimpleOperator>::type J = iterate(I); J; ++J)
      {
         Result += (*J) * triple_prod(herm(A.base()[J.index1()]), E, B[J.index2()],
                                      M.base().TransformsAs(), q);
      }
   }
   return Result;
}

MatrixOperator operator_prod(LinearAlgebra::HermitianProxy<MPStateComponent> const& A,
                             MatrixOperator const& E,
                             MPStateComponent const& B)
{
   DEBUG_PRECONDITION_EQUAL(A.base().Basis1(), E.Basis1());
   DEBUG_PRECONDITION_EQUAL(E.Basis2(), B.Basis1());
   DEBUG_PRECONDITION_EQUAL(A.base().SiteBasis(), B.SiteBasis());

   MatrixOperator Result(A.base().Basis2(), B.Basis2(), E.TransformsAs());
   QuantumNumbers::QuantumNumber Ident(B.GetSymmetryList());

   for (std::size_t i = 0; i < B.size(); ++i)
   {
      Result += triple_prod(herm(A.base()[i]), E, B[i], Ident, E.TransformsAs());
   }
   return Result;
}

MatrixOperator operator_prod(MPStateComponent const& A,
                             MatrixOperator const& E,
                             LinearAlgebra::HermitianProxy<MPStateComponent> const& B)
{
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), E.Basis1());
   DEBUG_PRECONDITION_EQUAL(E.Basis2(), B.base().Basis2());
   DEBUG_PRECONDITION_EQUAL(A.SiteBasis(), B.base().SiteBasis());

   MatrixOperator Result(A.Basis1(), B.base().Basis1(), E.TransformsAs());
   QuantumNumbers::QuantumNumber Ident(A.GetSymmetryList());

   for (std::size_t i = 0; i < A.size(); ++i)
   {
      Result += triple_prod(A[i], E, herm(B.base()[i]), Ident, E.TransformsAs());
   }
   return Result;
}

MatrixOperator operator_prod(LinearAlgebra::HermitianProxy<SimpleOperator> const& M,
                             LinearAlgebra::HermitianProxy<MPStateComponent> const& A,
                             MatrixOperator const& E,
                             MPStateComponent const& B)
{
   QuantumNumbers::QuantumNumberList ql = transform_targets(adjoint(M.base().TransformsAs()),
                                                            E.TransformsAs());
   PRECONDITION_EQUAL(ql.size(), 1)("Quantum number is ambiguous, must be specified.");

   return operator_prod(M, A, E, B, ql[0]);
}

MatrixOperator operator_prod(SimpleOperator const& M,
                             MPStateComponent const& A,
                             MatrixOperator const& E,
                             LinearAlgebra::HermitianProxy<MPStateComponent> const& B)
{
   QuantumNumbers::QuantumNumberList ql = transform_targets(M.TransformsAs(), E.TransformsAs());
   PRECONDITION_EQUAL(ql.size(), 1)("Quantum number is ambiguous, must be specified.")(ql);

   return operator_prod(M, A, E, B, ql[0]);
}

// tensor_sum

MPStateComponent tensor_sum(MPStateComponent const& A, MPStateComponent const& B,
                            SumBasis<VectorBasis> const& B1, SumBasis<VectorBasis> const& B2)
{
   PRECONDITION_EQUAL(A.SiteBasis(), B.SiteBasis());

   MPStateComponent Result(A.SiteBasis(), B1, B2);
   for (std::size_t s = 0; s < A.SiteBasis().size(); ++s)
   {
      Result[s] = tensor_sum(A[s], B[s], B1, B2);
   }
   return Result;
}

MPStateComponent tensor_row_sum(MPStateComponent const& A,
                                MPStateComponent const& B,
                                SumBasis<VectorBasis> const& B2)
{
   DEBUG_PRECONDITION_EQUAL(A.SiteBasis(), B.SiteBasis());
   DEBUG_PRECONDITION_EQUAL(A.Basis1(), B.Basis1());

   MPStateComponent Result(A.SiteBasis(), A.Basis1(), B2);
   for (std::size_t s = 0; s < A.SiteBasis().size(); ++s)
   {
      Result[s] = tensor_row_sum(A[s], B[s], B2);
   }
   return Result;
}

MPStateComponent tensor_col_sum(MPStateComponent const& A,
                                MPStateComponent const& B,
                                SumBasis<VectorBasis> const& B1)
{
   DEBUG_PRECONDITION_EQUAL(A.SiteBasis(), B.SiteBasis());
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), B.Basis2());

   MPStateComponent Result(A.SiteBasis(), B1, A.Basis2());
   for (std::size_t s = 0; s < A.SiteBasis().size(); ++s)
   {
      Result[s] = tensor_col_sum(A[s], B[s], B1);
   }
   return Result;
}

// prod

MPStateComponent prod(MPStateComponent const& A, MatrixOperator const& Op)
{
   PRECONDITION(is_scalar(Op.TransformsAs()))(Op.TransformsAs());
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), Op.Basis1());

   MPStateComponent Result(A.SiteBasis(), A.Basis1(), Op.Basis2());

   for (std::size_t s = 0; s < A.size(); ++s)
   {
      Result[s] = prod(A[s], Op, A[s].TransformsAs());
   }
   return Result;
}

MPStateComponent prod(MPStateComponent const& A, HermitianProxy<MatrixOperator> const& Op)
{
   PRECONDITION(is_scalar(Op.base().TransformsAs()))(Op.base().TransformsAs());
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), Op.base().Basis2());

   MPStateComponent Result(A.SiteBasis(), A.Basis1(), Op.base().Basis1());

   for (std::size_t s = 0; s < A.size(); ++s)
   {
      Result[s] = A[s] * Op;
   }
   return Result;
}

MPStateComponent prod(MatrixOperator const& Op, MPStateComponent const& A)
{
   PRECONDITION(is_scalar(Op.TransformsAs()))(Op.TransformsAs());
   DEBUG_PRECONDITION_EQUAL(Op.Basis2(), A.Basis1());

   MPStateComponent Result(A.SiteBasis(), Op.Basis1(), A.Basis2());

   for (std::size_t s = 0; s < A.size(); ++s)
   {
      Result[s] = prod(Op, A[s], A[s].TransformsAs());
   }
   return Result;
}

MPStateComponent prod(HermitianProxy<MatrixOperator> const& Op, MPStateComponent const& A)
{
   PRECONDITION(is_scalar(Op.base().TransformsAs()))(Op.base().TransformsAs());
   DEBUG_PRECONDITION_EQUAL(Op.base().Basis1(), A.Basis1());

   MPStateComponent Result(A.SiteBasis(), Op.base().Basis2(), A.Basis2());

   for (std::size_t s = 0; s < A.size(); ++s)
   {
      Result[s] = Op * A[s];
   }
   return Result;
}

// local_prod

MPStateComponent local_prod(MPStateComponent const& A, SimpleOperator const& x)
{
   DEBUG_CHECK_EQUAL(A.SiteBasis(), x.Basis1());
   MPStateComponent Result(x.Basis2(), A.Basis1(), A.Basis2());
   for (SimpleOperator::const_iterator I = iterate(x); I; ++I)
   {
      for (SimpleOperator::const_inner_iterator J = iterate(I); J; ++J)
      {
         Result[J.index2()] += (*J) * A[J.index1()];
      }
   }
   return Result;
}

MPStateComponent local_prod(SimpleOperator const& x, MPStateComponent const& A)
{
   // The conj_phase enters here for the same reason as prod(mp-operator, mp-state)
   DEBUG_CHECK_EQUAL(x.Basis2(), A.SiteBasis());
   MPStateComponent Result(x.Basis1(), A.Basis1(), A.Basis2());
   for (SimpleOperator::const_iterator I = iterate(x); I; ++I)
   {
      for (SimpleOperator::const_inner_iterator J = iterate(I); J; ++J)
      {
         Result[J.index1()] += conj_phase(x.Basis2()[J.index2()], x.TransformsAs(), x.Basis1()[J.index1()])
            * (*J) * A[J.index2()];
      }
   }
   return Result;
}

// triple_prod

MPStateComponent triple_prod(MatrixOperator const& Op1,
                             MPStateComponent const& A,
                             LinearAlgebra::HermitianProxy<MatrixOperator> const&Op2)
{
   MPStateComponent Result(A.SiteBasis(), Op1.Basis1(), Op2.base().Basis1());
   for (std::size_t s = 0; s < A.size(); ++s)
   {
      Result[s] = triple_prod(Op1, A[s], Op2);
   }
   return Result;
}

MPStateComponent triple_prod(LinearAlgebra::HermitianProxy<MatrixOperator> const& Op1,
                             MPStateComponent const& A,
                             MatrixOperator const&Op2)
{
   MPStateComponent Result(A.SiteBasis(), Op1.base().Basis2(), Op2.Basis2());
   for (std::size_t s = 0; s < A.size(); ++s)
   {
      Result[s] = triple_prod(Op1, A[s], Op2);
   }
   return Result;
}

MatrixOperator extract_diagonal(MPStateComponent const& A,
                                LinearAlgebra::HermitianProxy<MPStateComponent> const& B)
{
   MatrixOperator Result(A.Basis2(), B.base().Basis2());
   for (unsigned i = 0; i < A.Basis2().size(); ++i)
   {
      QuantumNumbers::QuantumNumber qi = A.Basis2()[i];
      for (unsigned j = 0; j < B.base().Basis2().size(); ++j)
      {
         QuantumNumbers::QuantumNumber qj = B.base().Basis2()[j];
         if (qi != qj)
            continue;

         for (unsigned a = 0; a < A.size(); ++a)
         {
            LinearAlgebra::const_inner_iterator<MatrixOperator>::type Aii = iterate_at(A[a].data(), i,i);
            LinearAlgebra::const_inner_iterator<MatrixOperator>::type Bjj = iterate_at(B.base()[a].data(), j,j);
            if (Aii && Bjj)
            {
               LinearAlgebra::Matrix<std::complex<double> > AA = *Aii, BB = *Bjj;
               LinearAlgebra::Matrix<std::complex<double> > R(size2(AA), size2(BB));
               for (unsigned x = 0; x < size1(R); ++x)
               {
                  for (unsigned y = 0; y < size2(R); ++y)
                  {
                     R(x,y) = AA(x,x) * herm(BB(y,y));
                  }
               }
               Result(i,j) = R;
            }
         }
      }
   }
   return Result;
}

MatrixOperator ExpandBasis1(MPStateComponent& A, Normalization n)
{
   ProductBasis<BasisList, VectorBasis> FullBasis1(A.SiteBasis(), A.Basis2());
   QuantumNumber Ident(A.GetSymmetryList());
   MPStateComponent Result(A.SiteBasis(), FullBasis1.Basis(), A.Basis2());
   for (std::size_t t = 0; t < FullBasis1.size(); ++t)
   {
      int s, b2;
      std::tie(s,b2) = FullBasis1.rmap(t);

      int Dim = FullBasis1.dim(t);
      DEBUG_CHECK_EQUAL(Dim, A.Basis2().dim(b2));

      // Make an identity matrix of the correct size
      LinearAlgebra::set_element(Result[s], t, b2, LinearAlgebra::identity_matrix<double>(Dim));
   }

   // check the normalization
   DEBUG_CHECK_CLOSE(norm_frob_sq(scalar_prod(Result, herm(Result))),
                     FullBasis1.total_degree());

   MatrixOperator Res = scalar_prod(A, herm(Result));
   A = Result;
   return Res;
}

MatrixOperator ExpandBasis2(MPStateComponent& A, Normalization n)
{
   ProductBasis<VectorBasis, BasisList> FullBasis2(A.Basis1(), adjoint(A.SiteBasis()));
   MPStateComponent Result(A.SiteBasis(), A.Basis1(), FullBasis2.Basis());
   for (std::size_t t = 0; t < FullBasis2.size(); ++t)
   {
      int s, b1;
      std::tie(b1,s) = FullBasis2.rmap(t);

      int Dim = FullBasis2.dim(t);
      DEBUG_CHECK_EQUAL(Dim, A.Basis1().dim(b1));

      // Make an identity matrix of the correct size
      set_element(Result[s], b1, t,
                                 std::sqrt(double(degree(FullBasis2[t])) / degree(A.Basis1()[b1]))
                                 * LinearAlgebra::identity_matrix<double>(Dim));
   }

   // check the normalization
   DEBUG_CHECK_CLOSE(norm_frob_sq(scalar_prod(herm(Result), Result)),
                     FullBasis2.total_degree());

   MatrixOperator Res = scalar_prod(herm(Result), A);
   A = Result;
   return Res;
}

MPStateComponent ConstructFromRightBasis(BasisList const& LocalBasis,
                                         VectorBasis const& RightBasis)
{
   ProductBasis<BasisList, VectorBasis> FullBasis1(LocalBasis, RightBasis);
   MPStateComponent Result(LocalBasis, FullBasis1.Basis(), RightBasis);
   for (std::size_t t = 0; t < FullBasis1.size(); ++t)
   {
      int s, b2;
      std::tie(s,b2) = FullBasis1.rmap(t);

      int Dim = FullBasis1.dim(t);
      DEBUG_CHECK_EQUAL(Dim, RightBasis.dim(b2));

      // Make an identity matrix of the correct size
      LinearAlgebra::set_element(Result[s], t, b2, LinearAlgebra::identity_matrix<double>(Dim));
   }
   return Result;
}

MPStateComponent ConstructFromLeftBasis(BasisList const& LocalBasis,
                                        VectorBasis const& LeftBasis)
{
   ProductBasis<VectorBasis, BasisList> FullBasis2(LeftBasis, adjoint(LocalBasis));
   MPStateComponent Result(LocalBasis, LeftBasis, FullBasis2.Basis());
   for (std::size_t t = 0; t < FullBasis2.size(); ++t)
   {
      int s, b1;
      std::tie(b1,s) = FullBasis2.rmap(t);

      int Dim = FullBasis2.dim(t);
      DEBUG_CHECK_EQUAL(Dim, LeftBasis.dim(b1));

      // Make an identity matrix of the correct size
      set_element(Result[s], b1, t,
                                 std::sqrt(double(degree(FullBasis2[t])) / degree(LeftBasis[b1]))
                                 * LinearAlgebra::identity_matrix<double>(Dim));
   }
   return Result;
}

MPStateComponent ShiftLocalBasis(MPStateComponent const& Op,
                                 QuantumNumber QL,
                                 QuantumNumber QM)
{
   CHECK_EQUAL(degree(QL), 1);
   CHECK_EQUAL(degree(QM), 1);

   // update the local basis
   BasisList NewLocal(Op.GetSymmetryList());
   for (unsigned i = 0; i < Op.LocalBasis().size(); ++i)
      NewLocal.push_back(transform_targets(Op.LocalBasis()[i], QL)[0]);

   // Update basis 2
   VectorBasis NewBasis2(Op.GetSymmetryList());
   for (unsigned i = 0; i < Op.Basis2().size(); ++i)
      NewBasis2.push_back(transform_targets(Op.Basis2()[i], QM)[0], Op.Basis2().dim(i));

   // Update basis 1 - in this case the quantum number shift is QM+QL
   QuantumNumber Q1 = transform_targets(QL,QM)[0];
   VectorBasis NewBasis1(Op.GetSymmetryList());
   for (unsigned i = 0; i < Op.Basis1().size(); ++i)
      NewBasis1.push_back(transform_targets(Op.Basis1()[i], Q1)[0], Op.Basis1().dim(i));

   // new component, and set the matrix elements
   MPStateComponent Result(NewLocal, NewBasis1, NewBasis2);
   for (unsigned i = 0; i < NewLocal.size(); ++i)
      Result[i].data() = Op[i].data();

   return Result;
}

MPStateComponent delta_shift(MPStateComponent const& Op, QuantumNumber const& q)
{
   CHECK_EQUAL(degree(q), 1);

   // Update basis 2
   VectorBasis NewBasis2(Op.GetSymmetryList());
   for (unsigned i = 0; i < Op.Basis2().size(); ++i)
      NewBasis2.push_back(transform_targets(Op.Basis2()[i], q)[0], Op.Basis2().dim(i));

   // Update basis 1
   VectorBasis NewBasis1(Op.GetSymmetryList());
   for (unsigned i = 0; i < Op.Basis1().size(); ++i)
      NewBasis1.push_back(transform_targets(Op.Basis1()[i], q)[0], Op.Basis1().dim(i));

   // new component, and set the matrix elements
   MPStateComponent Result(Op.LocalBasis(), NewBasis1, NewBasis2);
   for (unsigned i = 0; i < Op.size(); ++i)
      Result[i].data() = Op[i].data();

   return Result;
}

MPStateComponent ScaleBasisU1(MPStateComponent const& Op,
                              std::string const& Name,
                              double Factor)
{
   using QuantumNumbers::U1;

   int const s = Op.GetSymmetryList().WhichSymmetry(Name);
   if (s < 0)      // early return if Name does not name a quantum number
      return Op;

   CHECK_EQUAL(Op.GetSymmetryList().SymmetryType(s), "U(1)");

   // update the local basis
   BasisList NewLocal(Op.GetSymmetryList());
   for (unsigned i = 0; i < Op.LocalBasis().size(); ++i)
   {
      QuantumNumber q = Op.LocalBasis()[i];
      q.set<U1>(Name, U1(q.get<U1>(Name).x.to_double() * Factor));
      NewLocal.push_back(q);
   }

   // Update basis 2
   VectorBasis NewBasis2(Op.GetSymmetryList());
   for (unsigned i = 0; i < Op.Basis2().size(); ++i)
   {
      QuantumNumber q = Op.Basis2()[i];
      q.set<U1>(Name, U1(q.get<U1>(Name).x.to_double() * Factor));
      NewBasis2.push_back(q, Op.Basis2().dim(i));
   }

   // Update basis 1 - in this case the quantum number shift is QM+QL
   VectorBasis NewBasis1(Op.GetSymmetryList());
   for (unsigned i = 0; i < Op.Basis1().size(); ++i)
   {
      QuantumNumber q = Op.Basis1()[i];
      q.set<U1>(Name, U1(q.get<U1>(Name).x.to_double() * Factor));
      NewBasis1.push_back(q, Op.Basis1().dim(i));
   }

   // new component, and set the matrix elements
   MPStateComponent Result(NewLocal, NewBasis1, NewBasis2);
   for (unsigned i = 0; i < NewLocal.size(); ++i)
      Result[i].data() = Op[i].data();

   return Result;
}

MPStateComponent RenameSymmetry(MPStateComponent const& Op, SymmetryList const& NewSL)
{
   MPStateComponent Result(RenameSymmetry(Op.LocalBasis(), NewSL),
                           RenameSymmetry(Op.Basis1(), NewSL),
                           RenameSymmetry(Op.Basis2(), NewSL));

   // we assume here that the reduced matrix elements are identical
   for (unsigned i = 0; i < Result.LocalBasis().size(); ++i)
      Result[i].data() = Op[i].data();
   return Result;
}

MPStateComponent ReorderLocalBasis(MPStateComponent const& Op, std::list<int> const& NewOrder)
{
   // We could use the assert here, but perhaps it is useful to
   // allow things that are not strict permutations, eg miss out or repeat some basis states.
   // Note that the state won't be orthogonalized properly in that case.
   //   CHECK_EQUAL(Op.LocalBasis().size(), NewOrder.size());

   // update the local basis
   BasisList NewLocal(Op.GetSymmetryList());
   for (std::list<int>::const_iterator I = NewOrder.begin(); I != NewOrder.end(); ++I)
   {
      CHECK(*I < int(Op.LocalBasis().size()))(*I)(Op.LocalBasis().size());
      NewLocal.push_back(Op.LocalBasis()[*I]);
   }

   // construct the Result component
   MPStateComponent Result(NewLocal, Op.Basis1(), Op.Basis2());
   int i = 0;
   for (std::list<int>::const_iterator I = NewOrder.begin(); I != NewOrder.end(); ++I)
   {
      Result[i++] = Op[*I];
   }

   return Result;
}

MPStateComponent CoerceSymmetryList(MPStateComponent const& Op, SymmetryList const& NewSL)
{
   BasisList NewLocal = Op.LocalBasis();
   VectorBasis NewBasis1 = Op.Basis1();
   VectorBasis NewBasis2 = Op.Basis2();

   CoerceSymmetryList(NewLocal, NewSL);
   CoerceSymmetryList(NewBasis1, NewSL);
   CoerceSymmetryList(NewBasis2, NewSL);

   MPStateComponent Result(NewLocal, NewBasis1, NewBasis2);
   for (unsigned i = 0; i < Result.size(); ++i)
      Result[i].data() = Op[i].data();

   return Result;
}

MPStateComponent local_tensor_prod(MPStateComponent const& A, MPStateComponent const& B)
{
   Tensor::ProductBasis<BasisList, BasisList> PB(A.LocalBasis(), B.LocalBasis());
   MPStateComponent Result(PB.Basis(), A.Basis1(), B.Basis2());
   for (unsigned i = 0; i < A.LocalBasis().size(); ++i)
   {
      for (unsigned j = 0; j < B.LocalBasis().size(); ++j)
      {
         ProductBasis<BasisList, BasisList>::const_iterator IEnd = PB.end(i,j);
         ProductBasis<BasisList, BasisList>::const_iterator I = PB.begin(i,j);
         for ( ; I != IEnd; ++I)
         {
            Result[*I] = prod(A[i], B[j], PB.Basis()[*I]);
         }
      }
   }

   return Result;
}

MatrixOperator MakeRandomMatrixOperator(VectorBasis const& B1,
                                        VectorBasis const& B2,
                                        QuantumNumber q)
{
   MatrixOperator Result(B1, B2, q);
   for (std::size_t i = 0; i < B1.size(); ++i)
   {
      for (std::size_t j = 0; j < B2.size(); ++j)
      {
         if (is_transform_target(B2[j], q, B1[i]))
         {
            Result(i,j) = LinearAlgebra::random_matrix<double>(B1.dim(i), B2.dim(j));
         }
      }
   }
   return Result;
}
