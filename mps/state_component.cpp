// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mps/state_component.cpp
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

#include "state_component.h"
#include "quantumnumbers/quantumnumber.h"
#include "quantumnumbers/u1.h"
#include "tensor/tensorproduct.h"
#include <tuple>
#include "tensor/tensor_eigen.h"
#include "blas/matrix_utility.h"

StateComponent make_vacuum_state(QuantumNumbers::SymmetryList const& S)
{
   BasisList Vacuum = make_vacuum_basis(S);
   StateComponent Result(Vacuum, VectorBasis(Vacuum), VectorBasis(Vacuum));
   Result[0].emplace(0,0, 1,1,1.0); // MatrixOperator::value_type(1,1,1.0);
   return Result;
}

StateComponent make_vacuum_state(QuantumNumbers::QuantumNumber const& Q)
{
   BasisList Vacuum = make_vacuum_basis(Q.GetSymmetryList());
   BasisList vq(Q.GetSymmetryList());
   vq.push_back(Q);
   VectorBasis V(vq);
   StateComponent Result(Vacuum, V, V);
   Result[0].emplace(0,0, 1,1,1.0); // MatrixOperator::value_type(1,1,1.0);
   return Result;
}

StateComponent&
operator+=(StateComponent& x, StateComponent const& y)
{
   CHECK_EQUAL(x.LocalBasis(), y.LocalBasis());
   CHECK_EQUAL(x.Basis1(), y.Basis1());
   CHECK_EQUAL(x.Basis2(), y.Basis2());
   for (unsigned i = 0; i < x.size(); ++i)
   {
      x[i] += y[i];
   }
   return x;
}

StateComponent
operator+(StateComponent&& x, StateComponent const& y)
{
   CHECK_EQUAL(x.LocalBasis(), y.LocalBasis());
   CHECK_EQUAL(x.Basis1(), y.Basis1());
   CHECK_EQUAL(x.Basis2(), y.Basis2());
   StateComponent Result(std::move(x));
   for (unsigned i = 0; i < Result.size(); ++i)
   {
      Result[i] += y[i];
   }
   return Result;
}

StateComponent&
operator-=(StateComponent& x, StateComponent const& y)
{
   CHECK_EQUAL(x.LocalBasis(), y.LocalBasis());
   CHECK_EQUAL(x.Basis1(), y.Basis1());
   CHECK_EQUAL(x.Basis2(), y.Basis2());

   x.check_structure();
   y.check_structure();
   for (unsigned i = 0; i < x.size(); ++i)
   {
      x[i] -= y[i];
   }
   return x;
}

StateComponent
operator-(StateComponent&& x, StateComponent const& y)
{
   CHECK_EQUAL(x.LocalBasis(), y.LocalBasis());
   CHECK_EQUAL(x.Basis1(), y.Basis1());
   CHECK_EQUAL(x.Basis2(), y.Basis2());
   StateComponent Result(std::move(x));
   for (unsigned i = 0; i < Result.size(); ++i)
   {
      Result[i] -= y[i];
   }
   return Result;
}

StateComponent reflect(StateComponent const& S)
{
   StateComponent Result(S.LocalBasis(), adjoint(S.Basis2()), adjoint(S.Basis1()));
   int n=0;
   for (StateComponent::const_iterator I = S.begin(); I != S.end(); ++I)
   {
      Result[n++] = flip_conj(adjoint(*I));
   }
   return Result;
}


MatrixOperator
scalar_prod(StateComponent const& A, HermitianProxy<StateComponent> const& B)
{
   DEBUG_PRECONDITION_EQUAL(A.LocalBasis(), B.base().LocalBasis());
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), B.base().Basis2());

   QuantumNumbers::QuantumNumber Ident(A.GetSymmetryList());
   MatrixOperator Result(A.Basis1(), B.base().Basis1(), Ident);

   for (std::size_t s = 0; s < A.LocalBasis().size(); ++s)
   {
      Result += scalar_prod(A[s], herm(B.base()[s]));
   }
   return Result;
}

MatrixOperator
scalar_prod(HermitianProxy<StateComponent> const& A, StateComponent const& B)
{
   DEBUG_PRECONDITION_EQUAL(A.base().LocalBasis(), B.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(A.base().Basis1(), B.Basis1());

   QuantumNumbers::QuantumNumber Ident(B.GetSymmetryList());
   MatrixOperator Result(A.base().Basis2(), B.Basis2(), Ident);

   for (std::size_t s = 0; s < B.LocalBasis().size(); ++s)
   {
      Result += scalar_prod(herm(A.base()[s]), B[s]);
   }
   return Result;
}

SimpleOperator trace_prod(Tensor::HermitianProxy<StateComponent> const& A,
                          StateComponent const& B)
{
   QuantumNumbers::QuantumNumber Ident(B.GetSymmetryList());
   SimpleOperator Result(A.base().LocalBasis(), B.LocalBasis(), Ident);

   BasisList const& Sa = Result.Basis1();
   BasisList const& Sb = Result.Basis2();

   // we need a temporary Matrix to exist in device memory
   Matrix Temp(Sa.size(), Sb.size());
   for (std::size_t sa = 0; sa < Sa.size(); ++sa)
   {
      for (std::size_t sb = 0; sb < Sb.size(); ++sb)
      {
         //if (is_transform_target(Sb[sb], adjoint(Sa[sa]), Ident))
         if (Sb[sb] == Sa[sa])
            inner_prod(A.base()[sa], B[sb], Temp(sa, sb));
      }
   }

   // copy from device memory to the SimpleOperator
   for (std::size_t sa = 0; sa < Sa.size(); ++sa)
   {
      for (std::size_t sb = 0; sb < Sb.size(); ++sb)
      {
         if (Sb[sb] == Sa[sa])
	    Result(sa,sb) = get_wait(Temp(sa,sb));
      }
   }
   return Result;
}

SimpleOperator trace_prod(StateComponent const& A,
                          Tensor::HermitianProxy<StateComponent> const& B)
{
   QuantumNumbers::QuantumNumber Ident(A.GetSymmetryList());

   SimpleOperator Result(A.LocalBasis(), B.base().LocalBasis(), Ident);

   BasisList const& Sa = Result.Basis1();
   BasisList const& Sb = Result.Basis2();

   Matrix Temp(Sa.size(), Sb.size());
   for (std::size_t sa = 0; sa < Sa.size(); ++sa)
   {
      for (std::size_t sb = 0; sb < Sb.size(); ++sb)
      {
         //if (is_transform_target(adjoint(Sb[sb]), Sa[sa], Ident))
         if (Sb[sb] == Sa[sa])
         {
            inner_prod(B.base()[sb], A[sa], Temp(sa,sb));
         }
      }
   }
   for (std::size_t sa = 0; sa < Sa.size(); ++sa)
   {
      for (std::size_t sb = 0; sb < Sb.size(); ++sb)
      {
         if (Sb[sb] == Sa[sa])
	    Result(sa,sb) = get_wait(Temp(sa,sb));
      }
   }
   return Result;
}

#if 0
// does Result' = sum_{s',s} A[s'] * herm(B[s]) * M(s',s)
// M must transform as the identity representation.
// TODO: generalize this beyond identity op
// FIXME: This probably depends on the matrix elements of M having the correct norm convention
MatrixOperator operator_prod(SimpleOperator const& M,
                             StateComponent const& A,
                             Tensor::HermitianProxy<StateComponent> const& B)
{
   DEBUG_PRECONDITION_EQUAL(M.Basis1(), A.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.Basis2(), B.base().LocalBasis());
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
      }
   }
   return Result;
}
#endif

#if 0
MatrixOperator operator_prod(SimpleRedOperator const& M,
                             StateComponent const& A,
                             Tensor::HermitianProxy<StateComponent> const& B)
{
   MatrixOperator Result;
   for (SimpleRedOperator::const_iterator I = M.begin(); I != M.end(); ++I)
   {
      Result += operator_prod(*I, A, B);
   }
   return Result;
}

MatrixOperator operator_prod(Tensor::HermitianProxy<SimpleOperator> const& M,
                             Tensor::HermitianProxy<StateComponent> const& A,
                             StateComponent const& B)
{
   DEBUG_PRECONDITION_EQUAL(M.base().Basis1(), A.base().LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.base().Basis2(), B.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(A.base().Basis1(), B.Basis1());
   //DEBUG_PRECONDITION_EQUAL(M.TransformsAs(), QuantumNumbers::QuantumNumber(M.GetSymmetryList()));

   QuantumNumbers::QuantumNumber Ident(B.GetSymmetryList());
   MatrixOperator Result(A.base().Basis2(), B.Basis2(), adjoint(M.base().TransformsAs()));

   for (auto const& I : M.base())
   {
      for (auto const& J : I)
      {
         Result += std::conj(J.value) * triple_prod(herm(A.base()[I.row()]),
						    MatrixOperator::make_identity(B.Basis1()),
						    B[J.col()],
						    adjoint(Result.TransformsAs()),
						    Result.TransformsAs());
      }
   }
   return Result;
}

MatrixOperator operator_prod(Tensor::HermitianProxy<SimpleRedOperator> const& M,
                             Tensor::HermitianProxy<StateComponent> const& A,
                             StateComponent const& B)
{
   MatrixOperator Result;
   for (SimpleRedOperator::const_iterator I = M.base().begin(); I != M.base().end(); ++I)
   {
      Result += operator_prod(herm(*I), A, B);
   }
   return Result;
}

#endif

std::ostream& operator<<(std::ostream& out, StateComponent const& Psi)
{
  out << "Site basis:\n" << Psi.LocalBasis() << "Basis1:\n" << Psi.Basis1()
      << "Basis2:\n" << Psi.Basis2() << "\nNumber of components: " << Psi.Data.size() << '\n';
  for (int i = 0; i < Psi.Data.size(); ++i)
  {
     out << "Component " << i << ":\n";
     out << Psi.Data[i] << '\n';
  }
  return out;
}

MatrixOperator operator_prod(SimpleOperator const& M,
                             StateComponent const& A,
                             MatrixOperator const& E,
                             Tensor::HermitianProxy<StateComponent> const& B,
                             QuantumNumbers::QuantumNumber const& q)
{
   DEBUG_PRECONDITION_EQUAL(M.Basis1(), A.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.Basis2(), B.base().LocalBasis());
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), E.Basis1());
   DEBUG_PRECONDITION_EQUAL(E.Basis2(), B.base().Basis2());

   MatrixOperator Result(A.Basis1(), B.base().Basis1(), q);

   for (auto const& I : M)
   {
      for (auto const& J : I)
      {
         Result += J.value * triple_prod(A[I.row()], E, herm(B.base()[J.col()]),
					 M.TransformsAs(), q);
      }
   }
   return Result;
}

MatrixOperator operator_prod(SimpleRedOperator const& M,
                             StateComponent const& A,
                             MatrixOperator const& E,
                             Tensor::HermitianProxy<StateComponent> const& B,
                             QuantumNumbers::QuantumNumber const& q)
{
   MatrixOperator Result;
   for (SimpleRedOperator::const_iterator I = M.begin(); I != M.end(); ++I)
   {
      Result += operator_prod(*I, A, E, B, q);
   }
   return Result;
}


MatrixOperator operator_prod(Tensor::HermitianProxy<SimpleOperator> const& M,
                             Tensor::HermitianProxy<StateComponent> const& A,
                             MatrixOperator const& E,
                             StateComponent const& B,
                             QuantumNumbers::QuantumNumber const& q)
{
   DEBUG_PRECONDITION_EQUAL(M.base().Basis1(), A.base().LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.base().Basis2(), B.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(A.base().Basis1(), E.Basis1());
   DEBUG_PRECONDITION_EQUAL(E.Basis2(), B.Basis1());

   MatrixOperator Result(A.base().Basis2(), B.Basis2(), q);

   for (auto const& I : M.base())
   {
      for (auto const& J : I)
      {
         Result += herm(J.value) * triple_prod(herm(A.base()[I.row()]), E, B[J.col()],
					       M.base().TransformsAs(), q);
      }
   }
   return Result;
}

MatrixOperator operator_prod(Tensor::HermitianProxy<SimpleRedOperator> const& M,
                             Tensor::HermitianProxy<StateComponent> const& A,
                             MatrixOperator const& E,
                             StateComponent const& B,
                             QuantumNumbers::QuantumNumber const& q)
{
   MatrixOperator Result;
   for (SimpleRedOperator::const_iterator I = M.base().begin(); I != M.base().end(); ++I)
   {
      Result += operator_prod(herm(*I), A, E, B, q);
   }
   return Result;
}

MatrixOperator operator_prod(Tensor::HermitianProxy<StateComponent> const& A,
                             MatrixOperator const& E,
                             StateComponent const& B)
{
   DEBUG_PRECONDITION_EQUAL(A.base().Basis1(), E.Basis1());
   DEBUG_PRECONDITION_EQUAL(E.Basis2(), B.Basis1());
   DEBUG_PRECONDITION_EQUAL(A.base().LocalBasis(), B.LocalBasis());

   MatrixOperator Result(A.base().Basis2(), B.Basis2(), E.TransformsAs());
   QuantumNumbers::QuantumNumber Ident(B.GetSymmetryList());

   for (std::size_t i = 0; i < B.size(); ++i)
   {
      Result += triple_prod(herm(A.base()[i]), E, B[i], Ident, E.TransformsAs());
   }
   return Result;
}

MatrixOperator operator_prod(StateComponent const& A,
                             MatrixOperator const& E,
                             Tensor::HermitianProxy<StateComponent> const& B)
{
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), E.Basis1());
   DEBUG_PRECONDITION_EQUAL(E.Basis2(), B.base().Basis2());
   DEBUG_PRECONDITION_EQUAL(A.LocalBasis(), B.base().LocalBasis());

   MatrixOperator Result(A.Basis1(), B.base().Basis1(), E.TransformsAs());
   QuantumNumbers::QuantumNumber Ident(A.GetSymmetryList());

   for (std::size_t i = 0; i < A.size(); ++i)
   {
      Result += triple_prod(A[i], E, herm(B.base()[i]), Ident, E.TransformsAs());
   }
   return Result;
}

StateComponent operator_prod(Tensor::HermitianProxy<StateComponent> const& A,
                             StateComponent const& E,
                             StateComponent const& B)
{
   StateComponent Result(E.LocalBasis(), A.base().Basis2(), B.Basis2());
   for (unsigned q = 0; q < E.size(); ++q)
   {
      Result[q] = operator_prod(A, E[q], B);
   }
   return Result;
}

StateComponent operator_prod(StateComponent const& A,
                             StateComponent const& E,
                             Tensor::HermitianProxy<StateComponent> const& B)
{
   StateComponent Result(E.LocalBasis(), A.Basis1(), B.base().Basis1());
   for (unsigned q = 0; q < E.size(); ++q)
   {
      Result[q] = operator_prod(A, E[q], B);
   }
   return Result;
}

MatrixOperator
operator_prod_regular(StateComponent const& A,
                      MatrixOperator const& E,
                      Tensor::HermitianProxy<StateComponent> const& B)
{
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), E.Basis1());
   DEBUG_PRECONDITION_EQUAL(E.Basis2(), B.base().Basis2());
   DEBUG_PRECONDITION_EQUAL(A.LocalBasis(), B.base().LocalBasis());
   DEBUG_CHECK(norm_frob(A.front() - MatrixOperator::make_identity(A.front().Basis1())) < 1E-10);
   DEBUG_CHECK(norm_frob(B.base().back() - MatrixOperator::make_identity(B.base().back().Basis1())) < 1E-10);

   MatrixOperator Result(A.Basis1(), B.base().Basis1(), E.TransformsAs());
   QuantumNumbers::QuantumNumber Ident(A.GetSymmetryList());

   Result += scalar_prod(E, herm(B.base().front()));
   for (std::size_t i = 1; i < A.size()-1; ++i)
   {
      Result += triple_prod(A[i], E, herm(B.base()[i]), Ident, E.TransformsAs());
   }
   Result += prod(A.back(), E, E.TransformsAs());
   return Result;
}

MatrixOperator operator_prod(Tensor::HermitianProxy<SimpleOperator> const& M,
                             Tensor::HermitianProxy<StateComponent> const& A,
                             MatrixOperator const& E,
                             StateComponent const& B)
{
   QuantumNumbers::QuantumNumberList ql = transform_targets(adjoint(M.base().TransformsAs()),
                                                            E.TransformsAs());
   PRECONDITION_EQUAL(ql.size(), 1)("Quantum number is ambiguous, must be specified.");

   return operator_prod(M, A, E, B, ql[0]);
}

MatrixOperator operator_prod(SimpleOperator const& M,
                             StateComponent const& A,
                             MatrixOperator const& E,
                             Tensor::HermitianProxy<StateComponent> const& B)
{
   QuantumNumbers::QuantumNumberList ql = transform_targets(M.TransformsAs(), E.TransformsAs());
   PRECONDITION_EQUAL(ql.size(), 1)("Quantum number is ambiguous, must be specified.")(ql);

   return operator_prod(M, A, E, B, ql[0]);
}

// tensor_sum

StateComponent tensor_sum(StateComponent const& A, StateComponent const& B,
                            SumBasis<VectorBasis> const& B1, SumBasis<VectorBasis> const& B2)
{
   PRECONDITION_EQUAL(A.LocalBasis(), B.LocalBasis());

   StateComponent Result(A.LocalBasis(), B1, B2);
   for (std::size_t s = 0; s < A.LocalBasis().size(); ++s)
   {
      Result[s] = tensor_sum(A[s], B[s], B1, B2);
   }
   return Result;
}

StateComponent tensor_row_sum(StateComponent const& A,
                                StateComponent const& B,
                                SumBasis<VectorBasis> const& B2)
{
   DEBUG_PRECONDITION_EQUAL(A.LocalBasis(), B.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(A.Basis1(), B.Basis1());

   StateComponent Result(A.LocalBasis(), A.Basis1(), B2);
   for (std::size_t s = 0; s < A.LocalBasis().size(); ++s)
   {
      Result[s] = tensor_row_sum(A[s], B[s], B2);
   }
   return Result;
}

StateComponent tensor_col_sum(StateComponent const& A,
                                StateComponent const& B,
                                SumBasis<VectorBasis> const& B1)
{
   DEBUG_PRECONDITION_EQUAL(A.LocalBasis(), B.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), B.Basis2());

   StateComponent Result(A.LocalBasis(), B1, A.Basis2());
   for (std::size_t s = 0; s < A.LocalBasis().size(); ++s)
   {
      Result[s] = tensor_col_sum(A[s], B[s], B1);
   }
   return Result;
}

// prod

StateComponent prod(StateComponent const& A, MatrixOperator const& Op)
{
   PRECONDITION(is_scalar(Op.TransformsAs()))(Op.TransformsAs());
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), Op.Basis1());

   StateComponent Result(A.LocalBasis(), A.Basis1(), Op.Basis2());

   for (std::size_t s = 0; s < A.size(); ++s)
   {
      Result[s] = prod(A[s], Op, A[s].TransformsAs());
   }
   return Result;
}

StateComponent prod(StateComponent const& A, HermitianProxy<MatrixOperator> const& Op)
{
   PRECONDITION(is_scalar(Op.base().TransformsAs()))(Op.base().TransformsAs());
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), Op.base().Basis2());

   StateComponent Result(A.LocalBasis(), A.Basis1(), Op.base().Basis1());

   for (std::size_t s = 0; s < A.size(); ++s)
   {
      Result[s] = A[s] * Op;
   }
   return Result;
}

StateComponent prod(MatrixOperator const& Op, StateComponent const& A)
{
   PRECONDITION(is_scalar(Op.TransformsAs()))(Op.TransformsAs());
   DEBUG_PRECONDITION_EQUAL(Op.Basis2(), A.Basis1());

   StateComponent Result(A.LocalBasis(), Op.Basis1(), A.Basis2());

   for (std::size_t s = 0; s < A.size(); ++s)
   {
      Result[s] = prod(Op, A[s], A[s].TransformsAs());
   }
   return Result;
}

StateComponent prod(HermitianProxy<MatrixOperator> const& Op, StateComponent const& A)
{
   PRECONDITION(is_scalar(Op.base().TransformsAs()))(Op.base().TransformsAs());
   DEBUG_PRECONDITION_EQUAL(Op.base().Basis1(), A.Basis1());

   StateComponent Result(A.LocalBasis(), Op.base().Basis2(), A.Basis2());

   for (std::size_t s = 0; s < A.size(); ++s)
   {
      Result[s] = Op * A[s];
   }
   return Result;
}

// local_prod

StateComponent local_prod(StateComponent const& A, SimpleOperator const& x)
{
   DEBUG_CHECK_EQUAL(A.LocalBasis(), x.Basis1());
   StateComponent Result(x.Basis2(), A.Basis1(), A.Basis2());
   for (auto const& I : x)
   {
      for (auto const& J : I)
      {
         Result[J.col()] += J.value * A[I.row()];
      }
   }
   return Result;
}

StateComponent local_prod(SimpleOperator const& x, StateComponent const& A)
{
   // The conj_phase enters here for the same reason as prod(mp-operator, mp-state)
   DEBUG_CHECK_EQUAL(x.Basis2(), A.LocalBasis());
   StateComponent Result(x.Basis1(), A.Basis1(), A.Basis2());
   for (auto const& I : x)
   {
      for (auto const& J : I)
      {
	 Result[I.row()] += conj_phase(x.Basis2()[J.col()], x.TransformsAs(), x.Basis1()[I.row()])
            * J.value * A[J.col()];
      }
   }
   return Result;
}

// triple_prod

StateComponent triple_prod(MatrixOperator const& Op1,
                             StateComponent const& A,
                             Tensor::HermitianProxy<MatrixOperator> const&Op2)
{
   StateComponent Result(A.LocalBasis(), Op1.Basis1(), Op2.base().Basis1());
   for (std::size_t s = 0; s < A.size(); ++s)
   {
      Result[s] = triple_prod(Op1, A[s], Op2);
   }
   return Result;
}

StateComponent triple_prod(Tensor::HermitianProxy<MatrixOperator> const& Op1,
                             StateComponent const& A,
                             MatrixOperator const&Op2)
{
   StateComponent Result(A.LocalBasis(), Op1.base().Basis2(), Op2.Basis2());
   for (std::size_t s = 0; s < A.size(); ++s)
   {
      Result[s] = triple_prod(Op1, A[s], Op2);
   }
   return Result;
}

#if 0
MatrixOperator extract_diagonal(StateComponent const& A,
                                Tensor::HermitianProxy<StateComponent> const& B)
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
#endif

MatrixOperator ExpandBasis1(StateComponent& A)
{
   ProductBasis<BasisList, VectorBasis> FullBasis1(A.LocalBasis(), A.Basis2());
   QuantumNumber Ident(A.GetSymmetryList());
   StateComponent Result(A.LocalBasis(), FullBasis1.Basis(), A.Basis2());
   for (std::size_t t = 0; t < FullBasis1.size(); ++t)
   {
      int s, b2;
      std::tie(s,b2) = FullBasis1.rmap(t);

      int Dim = FullBasis1.dim(t);
      DEBUG_CHECK_EQUAL(Dim, A.Basis2().dim(b2));

      // Make an identity matrix of the correct size
      Result[s].insert(t, b2, MatrixOperator::value_type::make_identity(Dim));
   }

   // check the normalization
   DEBUG_CHECK_CLOSE(norm_frob_sq(scalar_prod(Result, herm(Result))),
                     FullBasis1.total_degree())(Result)(scalar_prod(Result, herm(Result)));

   MatrixOperator Res = scalar_prod(A, herm(Result));
   A = std::move(Result);
   return Res;
}

MatrixOperator ExpandBasis2(StateComponent& A)
{
   ProductBasis<VectorBasis, BasisList> FullBasis2(A.Basis1(), adjoint(A.LocalBasis()));
   StateComponent Result(A.LocalBasis(), A.Basis1(), FullBasis2.Basis());
   for (std::size_t t = 0; t < FullBasis2.size(); ++t)
   {
      int s, b1;
      std::tie(b1,s) = FullBasis2.rmap(t);

      int Dim = FullBasis2.dim(t);
      DEBUG_CHECK_EQUAL(Dim, A.Basis1().dim(b1));

      // Make an identity matrix of the correct size
      Result[s].insert(b1, t, std::sqrt(double(degree(FullBasis2[t])) / degree(A.Basis1()[b1]))
		       * MatrixOperator::value_type::make_identity(Dim));
   }

   // check the normalization
   DEBUG_CHECK_CLOSE(norm_frob_sq(scalar_prod(herm(Result), Result)),
                     FullBasis2.total_degree())(scalar_prod(herm(Result), Result));

   MatrixOperator Res = scalar_prod(herm(Result), A);
   A = std::move(Result);
   return Res;
}

MatrixOperator ExpandBasis1Used(StateComponent& A, std::vector<int> const& Used)
{
   // get a set of the local basis states that are actually used
   std::vector<int> UsedBasisElements;
   BasisList UsedBasis(A.GetSymmetryList());
   for (unsigned i = 0; i < A.size(); ++i)
   {
      if (Used[i])
      {
         UsedBasisElements.push_back(i);
         UsedBasis.push_back(A.LocalBasis()[i]);
      }
   }
   ProductBasis<BasisList, VectorBasis> FullBasis1(UsedBasis, A.Basis2());
   QuantumNumber Ident(A.GetSymmetryList());
   StateComponent Result(A.LocalBasis(), FullBasis1.Basis(), A.Basis2());
   for (std::size_t t = 0; t < FullBasis1.size(); ++t)
   {
      int s, b2;
      std::tie(s,b2) = FullBasis1.rmap(t);
      s = UsedBasisElements[s]; // reverse map to the complete local basis

      int Dim = FullBasis1.dim(t);
      DEBUG_CHECK_EQUAL(Dim, A.Basis2().dim(b2));

      // Make an identity matrix of the correct size
      Result[s].insert(t,b2, MatrixOperator::value_type::make_identity(Dim));
   }

   // check the normalization
   DEBUG_CHECK_CLOSE(norm_frob_sq(scalar_prod(Result, herm(Result))),
                     FullBasis1.total_degree());

   MatrixOperator Res = scalar_prod(A, herm(Result));
   A = std::move(Result);
   return Res;
}

MatrixOperator ExpandBasis2Used(StateComponent& A, std::vector<int> const& Used)
{
   // get a set of the local basis states that are actually used
   std::vector<int> UsedBasisElements;
   BasisList UsedBasis(A.GetSymmetryList());
   for (unsigned i = 0; i < A.size(); ++i)
   {
      if (Used[i])
      {
         UsedBasisElements.push_back(i);
         UsedBasis.push_back(A.LocalBasis()[i]);
      }
   }
   CHECK(UsedBasis.size() != 0)(UsedBasis);

   ProductBasis<VectorBasis, BasisList> FullBasis2(A.Basis1(), adjoint(UsedBasis));
   StateComponent Result(A.LocalBasis(), A.Basis1(), FullBasis2.Basis());
   for (std::size_t t = 0; t < FullBasis2.size(); ++t)
   {
      int s, b1;
      std::tie(b1,s) = FullBasis2.rmap(t);
      s = UsedBasisElements[s];

      int Dim = FullBasis2.dim(t);
      DEBUG_CHECK_EQUAL(Dim, A.Basis1().dim(b1));

      // Make an identity matrix of the correct size
      Result[s].insert(b1,t, std::sqrt(double(degree(FullBasis2[t])) / degree(A.Basis1()[b1]))
		       * MatrixOperator::value_type::make_identity(Dim));
   }

   // check the normalization
   DEBUG_CHECK_CLOSE(norm_frob_sq(scalar_prod(herm(Result), Result)),
                     FullBasis2.total_degree());

   MatrixOperator Res = scalar_prod(herm(Result), A);
   A = std::move(Result);
   return Res;
}

#if 0
std::pair<MatrixOperator, SimpleStateComponent>
ExpandBasis1_(StateComponent const& A)
{
   ProductBasis<BasisList, VectorBasis> FullBasis1(A.LocalBasis(), A.Basis2());
   QuantumNumber Ident(A.GetSymmetryList());
   SimpleStateComponent Result(A.LocalBasis(), FullBasis1.Basis(), A.Basis2());
   for (std::size_t t = 0; t < FullBasis1.size(); ++t)
   {
      int s, b2;
      std::tie(s,b2) = FullBasis1.rmap(t);

      int Dim = FullBasis1.dim(t);
      DEBUG_CHECK_EQUAL(Dim, A.Basis2().dim(b2));

      // Make an identity matrix of the correct size
      Result[s].insert(t, b2, LinearAlgebra::ScalarMatrix<std::complex<double>>(Dim, Dim, 1.0));
   }

   // check the normalization
   DEBUG_CHECK_CLOSE(norm_frob_sq(scalar_prod(Result, herm(Result))),
                     FullBasis1.total_degree());

   MatrixOperator Res = scalar_prod(A, herm(Result));
   return std::make_pair(Res, Result);
}
#endif

std::pair<RealDiagonalOperator, MatrixOperator>
OrthogonalizeBasis2(StateComponent& A)
{
   // TODO: optimize this implementation

   MatrixOperator U, Vh;
   RealDiagonalOperator D;
   MatrixOperator M = ExpandBasis2(A);
   SVD(M, U, D, Vh);
   A = prod(A, U);
   return std::make_pair(std::move(D), std::move(Vh));
}

std::pair<MatrixOperator, RealDiagonalOperator>
OrthogonalizeBasis1(StateComponent& A)
{
   // TODO: optimize this implementation
   MatrixOperator U, Vh;
   RealDiagonalOperator D;
   MatrixOperator M = ExpandBasis1(A);
   SVD(M, U, D, Vh);
   A = prod(Vh, A);
   return std::make_pair(std::move(U), std::move(D));
}

StateComponent ConstructFromRightBasis(BasisList const& LocalBasis,
                                         VectorBasis const& RightBasis)
{
   ProductBasis<BasisList, VectorBasis> FullBasis1(LocalBasis, RightBasis);
   StateComponent Result(LocalBasis, FullBasis1.Basis(), RightBasis);
   for (std::size_t t = 0; t < FullBasis1.size(); ++t)
   {
      int s, b2;
      std::tie(s,b2) = FullBasis1.rmap(t);

      int Dim = FullBasis1.dim(t);
      DEBUG_CHECK_EQUAL(Dim, RightBasis.dim(b2));

      // Make an identity matrix of the correct size
      Result[s].insert(t,b2, MatrixOperator::value_type::make_identity(Dim));
   }
   return Result;
}

StateComponent ConstructFromLeftBasis(BasisList const& LocalBasis,
                                        VectorBasis const& LeftBasis)
{
   ProductBasis<VectorBasis, BasisList> FullBasis2(LeftBasis, adjoint(LocalBasis));
   StateComponent Result(LocalBasis, LeftBasis, FullBasis2.Basis());
   for (std::size_t t = 0; t < FullBasis2.size(); ++t)
   {
      int s, b1;
      std::tie(b1,s) = FullBasis2.rmap(t);

      int Dim = FullBasis2.dim(t);
      DEBUG_CHECK_EQUAL(Dim, LeftBasis.dim(b1));

      // Make an identity matrix of the correct size
      Result[s].insert(b1,t, std::sqrt(double(degree(FullBasis2[t])) / degree(LeftBasis[b1]))
		       * MatrixOperator::value_type::make_identity(Dim));
   }
   return Result;
}

StateComponent ShiftLocalBasis(StateComponent&& Op,
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
   StateComponent Result(NewLocal, NewBasis1, NewBasis2);
   for (unsigned i = 0; i < NewLocal.size(); ++i)
      Result[i].data() = std::move(Op[i].data());

   return Result;
}

StateComponent delta_shift(StateComponent const& Op, QuantumNumber const& q)
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
   StateComponent Result(Op.LocalBasis(), NewBasis1, NewBasis2);
   for (unsigned i = 0; i < Op.size(); ++i)
      Result[i].data() = copy(Op[i].data());

   return Result;
}

StateComponent delta_shift(StateComponent&& Op, QuantumNumber const& q)
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
   StateComponent Result(Op.LocalBasis(), NewBasis1, NewBasis2);
   for (unsigned i = 0; i < Op.size(); ++i)
      Result[i].data() = std::move(Op[i].data());

   return Result;
}

#if 0
StateComponent ScaleBasisU1(StateComponent const& Op,
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
   StateComponent Result(NewLocal, NewBasis1, NewBasis2);
   for (unsigned i = 0; i < NewLocal.size(); ++i)
      Result[i].data() = Op[i].data();

   return Result;
}
#endif

#if 0
MatrixOperator RenameSymmetry(MatrixOperator const& Op, SymmetryList const& NewSL)
{
   MatrixOperator Result(RenameSymmetry(Op.Basis1(), NewSL),
                         RenameSymmetry(Op.Basis2(), NewSL),
                         RenameSymmetry(Op.TransformsAs(), NewSL));
   Result.data() = Op.data();
   return Result;
}

StateComponent RenameSymmetry(StateComponent const& Op, SymmetryList const& NewSL)
{
   StateComponent Result(RenameSymmetry(Op.LocalBasis(), NewSL),
                           RenameSymmetry(Op.Basis1(), NewSL),
                           RenameSymmetry(Op.Basis2(), NewSL));

   // we assume here that the reduced matrix elements are identical
   for (unsigned i = 0; i < Result.LocalBasis().size(); ++i)
      Result[i].data() = Op[i].data();
   return Result;
}
#endif

StateComponent ReorderLocalBasis(StateComponent const& Op, std::list<int> const& NewOrder)
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
   StateComponent Result(NewLocal, Op.Basis1(), Op.Basis2());
   int i = 0;
   for (std::list<int>::const_iterator I = NewOrder.begin(); I != NewOrder.end(); ++I)
   {
      Result[i++] = copy(Op[*I]);
   }

   return Result;
}

StateComponent CoerceSymmetryList(StateComponent&& Op, SymmetryList const& NewSL)
{
   BasisList NewLocal = CoerceSymmetryList(Op.LocalBasis(), NewSL);
   VectorBasis NewBasis1 = CoerceSymmetryList(Op.Basis1(), NewSL);
   VectorBasis NewBasis2 = CoerceSymmetryList(Op.Basis2(), NewSL);

   StateComponent Result(NewLocal, NewBasis1, NewBasis2);
   for (unsigned i = 0; i < Result.size(); ++i)
      Result[i].data() = std::move(Op[i].data());

   return Result;
}

StateComponent CoerceSymmetryList(StateComponent const& Op, SymmetryList const& NewSL)
{
   BasisList NewLocal = CoerceSymmetryList(Op.LocalBasis(), NewSL);
   VectorBasis NewBasis1 = CoerceSymmetryList(Op.Basis1(), NewSL);
   VectorBasis NewBasis2 = CoerceSymmetryList(Op.Basis2(), NewSL);

   StateComponent Result(NewLocal, NewBasis1, NewBasis2);
   for (unsigned i = 0; i < Result.size(); ++i)
      Result[i].data() = copy(Op[i].data());

   return Result;
}

StateComponent local_tensor_prod(StateComponent const& A, StateComponent const& B)
{
   Tensor::ProductBasis<BasisList, BasisList> PB(A.LocalBasis(), B.LocalBasis());
   StateComponent Result(PB.Basis(), A.Basis1(), B.Basis2());
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
	    Matrix Temp(B1.dim(i), B2.dim(j));
	    set_wait(Temp, blas::random_matrix<blas::numeric_type_of_t<Matrix>>(B1.dim(i), B2.dim(j)));
	    Result.insert(i,j, Temp);
         }
      }
   }
   return Result;
}

StateComponent
MakeRandomStateComponent(BasisList const& Local, VectorBasis const& B1, VectorBasis const& B2)
{
   StateComponent Result(Local, B1, B2);
   for (unsigned i = 0; i < Result.size(); ++i)
   {
      Result[i] = MakeRandomMatrixOperator(B1, B2, Local[i]);
   }
   return Result;
}
