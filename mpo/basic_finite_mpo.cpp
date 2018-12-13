// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mpo/finite_mpo.cpp
//
// Copyright (C) 2013-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "basic_finite_mpo.h"
#include "pstream/variant.h"
#include "tensor/tensor_exponential.h"
#include "lattice/siteoperator-parser.h"
#include <algorithm>
#include <stack>

using QuantumNumbers::QuantumNumber;

bool BasicFiniteMPO::is_irreducible() const
{
   return this->is_null() || (this->Basis1().size() == 1 && QuantumNumbers::is_scalar(this->qn2()));
}

bool BasicFiniteMPO::is_scalar() const
{
   return this->is_null() || (this->Basis1().size() == 1 && QuantumNumbers::is_scalar(this->qn1())
                              && QuantumNumbers::is_scalar(this->qn2()));
}

#if !defined(DISABLE_FINITE_MPO_TRANSFORMS_AS)
QuantumNumbers::QuantumNumber
BasicFiniteMPO::TransformsAs() const
{
   if (this->is_null())
      return QuantumNumbers::QuantumNumber();
   CHECK(this->is_irreducible());
   return this->Basis1()[0];
}
#endif

QuantumNumbers::QuantumNumber
BasicFiniteMPO::qn1() const
{
   if (this->is_null())
      return QuantumNumbers::QuantumNumber();
   CHECK_EQUAL(this->Basis1().size(), 1);
   return this->Basis1()[0];
}

QuantumNumbers::QuantumNumber
BasicFiniteMPO::qn2() const
{
   if (this->is_null())
      return QuantumNumbers::QuantumNumber();
   DEBUG_CHECK_EQUAL(this->Basis2().size(), 1);  // DEBUG since this should always be true
   return this->Basis2()[0];
}

PStream::opstream& operator<<(PStream::opstream& out, BasicFiniteMPO const& op)
{
   return out << op.data();
}

PStream::ipstream& operator>>(PStream::ipstream& in, BasicFiniteMPO& op)
{
   return in >> op.data();
}

void print_structure(BasicFiniteMPO const& Op, std::ostream& out, double UnityEpsilon)
{
   out << "BasicFiniteMPO has " << Op.size() << " sites\n";
   for (int i = 0; i < Op.size(); ++i)
   {
      out << "Site " << i << " dimension " << Op[i].size1() << " x " << Op[i].size2() << '\n';
      print_structure(Op[i], out, UnityEpsilon);
   }
}

BasicFiniteMPO
join(BasicFiniteMPO const& Op1, BasicFiniteMPO const& Op2)
{
   if (Op1.is_null())
      return Op2;
   if (Op2.is_null())
      return Op1;
   CHECK_EQUAL(Op1.Basis2(), Op2.Basis1());
   BasicFiniteMPO Result(Op1.size() + Op2.size());
   for (int i = 0; i < Op1.size(); ++i)
   {
      Result[i] = copy(Op1[i]);
   }
   for (int i = 0; i < Op2.size(); ++i)
   {
      Result[i+Op1.size()] = copy(Op2[i]);
   }
   return Result;
}

BasicFiniteMPO
repeat(BasicFiniteMPO const& Op, int Count)
{
   CHECK_EQUAL(Op.Basis1(), Op.Basis2());
   BasicFiniteMPO Result(Op.size()*Count);
   for (int i = 0; i < Result.size(); ++i)
   {
      Result[i] = copy(Op[i%Op.size()]);
   }
   return Result;
}

void optimize(BasicFiniteMPO& Op)
{
   if (Op.size() < 2)
      return;

#if 0
#if !defined(NDEBUG)
   SimpleRedOperator X = coarse_grain(Op);
#endif
#endif

   bool Reduced = true; // flag to indicate that we reduced a dimension
   // loop until we do a complete sweep with no reduction in dimensions
   while (Reduced)
   {
      Reduced = false;

      // Working left to right, optimize the Basis2
      SimpleOperator T = TruncateBasis2(Op.front());
      if (T.size1() != T.size2())
         Reduced = true;
      for (int i = 1; i < Op.size()-1; ++i)
      {
         Op[i] = T * Op[i];
         T = TruncateBasis2(Op[i]);
         if (T.size1() != T.size2())
            Reduced = true;
      }
      Op.back() = T * Op.back();

      // Working right to left, optimize Basis1
      T = TruncateBasis1(Op.back());
      if (T.size1() != T.size2())
         Reduced = true;
      for (int i = Op.size()-2; i >= 1; --i)
      {
         Op[i] = Op[i] * T;
         T = TruncateBasis1(Op[i]);
         if (T.size1() != T.size2())
            Reduced = true;
      }
      Op.front() = Op.front() * T;
   }

#if 0
#if !defined(NDEBUG)
   SimpleRedOperator Y = coarse_grain(Op);
   CHECK(norm_frob(X-Y) < 1E-10*Y.Basis1().total_degree())("failure in MPO optimize()!")(X)(Y);
#endif
#endif
}

#if 0
void qr_optimize(BasicFiniteMPO& Op)
{
   if (Op.size() < 2)
      return;

#if 0
#if !defined(NDEBUG)
   SimpleRedOperator X = coarse_grain(Op);
#endif
#endif

   double const Eps = 1E-8;

   //   TRACE(Op);

   bool Reduced = true; // flag to indicate that we reduced a dimension
   // loop until we do a complete sweep with no reduction in dimensions
   bool First = true;
   bool Second = true;

   First = false;
   while (Reduced || Second)
   {
      Reduced = false;

      OperatorComponent Op2 = copy(Op.front());
      if (!First && Second)
      {
         TRACE("XXXXX");
      }
      SimpleOperator T2 = TruncateBasis2MkII(Op2, First ? 0.0 : Eps);
      //TRACE(norm_frob(Op.front() - Op2*T2));

      // Working left to right, optimize the Basis2
      SimpleOperator T = TruncateBasis2MkII(Op.front(), First ? 0.0 : Eps);
      if (T.size1() != T.size2())
         Reduced = true;
      for (int i = 1; i < Op.size()-1; ++i)
      {
         Op[i] = T * Op[i];
         T = TruncateBasis2MkII(Op[i]);
         if (T.size1() != T.size2())
            Reduced = true;
      }
      Op.back() = T * Op.back();

      // Working right to left, optimize Basis1
      T = TruncateBasis1MkII(Op.back(), Eps);
      if (T.size1() != T.size2())
         Reduced = true;
      for (int i = Op.size()-2; i >= 1; --i)
      {
         Op[i] = Op[i] * T;
         T = TruncateBasis1MkII(Op[i], Eps);
         if (T.size1() != T.size2())
            Reduced = true;
      }
      Op.front() = Op.front() * T;

      if (!First) Second = false;
      First = false;
   }

   SimpleOperator Overlaps = local_inner_prod(herm(Op[Op.size()/2]), Op[Op.size()/2]);
   TRACE(Overlaps);

   //TRACE(Op);

#if 0
#if !defined(NDEBUG)
   SimpleRedOperator Y = coarse_grain(Op);
   CHECK(norm_frob(X-Y) < 1E-10*Y.Basis1().total_degree())("failure in MPO optimize()!")(X)(Y);
#endif
#endif
}
#endif

BasicFiniteMPO&
operator+=(BasicFiniteMPO& x, BasicFiniteMPO const& y)
{
   if (x.empty())
   {
      x = y;
      return x;
   }
   if (y.empty())
      return x;

   SumBasis<BasisList> sb1(x.Basis1(), y.Basis1());
   SimpleOperator C = CollapseBasis(sb1.Basis());
   SumBasis<BasisList> sb2(x.front().Basis2(), y.front().Basis2());
   x.front() = C * tensor_sum(x.front(), y.front(), sb1, sb2);
   for (int i = 1; i < x.size(); ++i)
   {
      sb1 = sb2;
      sb2 = SumBasis<BasisList>(x[i].Basis2(), y[i].Basis2());
      x[i] = tensor_sum(x[i], y[i], sb1, sb2);
   }
   x.back() = x.back() * herm(CollapseBasis(x.back().Basis2()));
   optimize(x);
   return x;
}

BasicFiniteMPO
operator+(BasicFiniteMPO const& x, BasicFiniteMPO const& y)
{
   if (x.empty())
      return y;
   if (y.empty())
      return x;

   BasicFiniteMPO Result(x.size());

   SumBasis<BasisList> sb1(x.Basis1(), y.Basis1());
   SimpleOperator C = CollapseBasis(sb1.Basis());
   SumBasis<BasisList> sb2(x.front().Basis2(), y.front().Basis2());
   Result.front() = C * tensor_sum(x.front(), y.front(), sb1, sb2);
   for (int i = 1; i < x.size(); ++i)
   {
      sb1 = sb2;
      sb2 = SumBasis<BasisList>(x[i].Basis2(), y[i].Basis2());
      Result[i] = tensor_sum(x[i], y[i], sb1, sb2);
   }
   Result.back() = Result.back() * herm(CollapseBasis(Result.back().Basis2()));
   optimize(Result);
   return Result;
}

BasicFiniteMPO
operator-(BasicFiniteMPO const& x, BasicFiniteMPO const& y)
{
   if (x.empty())
      return y;
   if (y.empty())
      return x;

   BasicFiniteMPO Result(x.size());

   SumBasis<BasisList> sb1(x.Basis1(), y.Basis1());
   SimpleOperator C = CollapseBasis(sb1.Basis());
   SumBasis<BasisList> sb2(x.front().Basis2(), y.front().Basis2());
   Result.front() = C * tensor_sum(x.front(), -1.0 * y.front(), sb1, sb2);
   for (int i = 1; i < x.size(); ++i)
   {
      sb1 = sb2;
      sb2 = SumBasis<BasisList>(x[i].Basis2(), y[i].Basis2());
      Result[i] = tensor_sum(x[i], y[i], sb1, sb2);
   }
   Result.back() = Result.back() * herm(CollapseBasis(Result.back().Basis2()));
   optimize(Result);
   return Result;
}

BasicFiniteMPO
operator-(BasicFiniteMPO const& x)
{
   BasicFiniteMPO Result(x);
   if (!Result.empty())
      Result.front() *= -1.0;
   return Result;
}

BasicFiniteMPO&
operator-=(BasicFiniteMPO& x, BasicFiniteMPO const& y)
{
   if (x.empty())
   {
      x = y;
      return x;
   }
   if (y.empty())
      return x;

   SumBasis<BasisList> sb1(x.Basis1(), y.Basis1());
   SimpleOperator C = CollapseBasis(sb1.Basis());
   SumBasis<BasisList> sb2(x.front().Basis2(), y.front().Basis2());
   x.front() = C * tensor_sum(x.front(), -1.0 * y.front(), sb1, sb2);
   for (int i = 1; i < x.size(); ++i)
   {
      sb1 = sb2;
      sb2 = SumBasis<BasisList>(x[i].Basis2(), y[i].Basis2());
      x[i] = tensor_sum(x[i], y[i], sb1, sb2);
   }
   x.back() = x.back() * herm(CollapseBasis(x.back().Basis2()));
   optimize(x);
   return x;
}

BasicFiniteMPO&
operator*=(BasicFiniteMPO& m, double x)
{
   if (!m.empty())
      m.front() *= x;
   return m;
}

BasicFiniteMPO&
operator*=(BasicFiniteMPO& m, std::complex<double> x)
{
   if (!m.empty())
      m.front() *= x;
   return m;
}

BasicFiniteMPO
operator*(BasicFiniteMPO const& m, double x)
{
   BasicFiniteMPO Result(m);
   if (!Result.empty())
      Result.front() *= x;
   return Result;
}

BasicFiniteMPO
operator*(BasicFiniteMPO const& m, std::complex<double> x)
{
   BasicFiniteMPO Result(m);
   if (!Result.empty())
      Result.front() *= x;
   return Result;
}

BasicFiniteMPO
operator*(double x, BasicFiniteMPO const& m)
{
   BasicFiniteMPO Result(m);
   if (!Result.empty())
      Result.front() *= x;
   return Result;
}

BasicFiniteMPO
operator*(std::complex<double> x, BasicFiniteMPO const& m)
{
   BasicFiniteMPO Result(m);
   if (!Result.empty())
      Result.front() *= x;
   return Result;
}

BasicFiniteMPO
operator*(BasicFiniteMPO const& x, BasicFiniteMPO const& y)
{
   if (x.empty())
      return x;
   if (y.empty())
      return y;
   CHECK_EQUAL(x.size(), y.size());
   BasicFiniteMPO Result(x.size());
   for (int i = 0; i < x.size(); ++i)
   {
      Result[i] = aux_tensor_prod(x[i], y[i]);
   }
   Result.front() = CollapseBasis(Result.front().Basis1()) * Result.front();
   Result.back() = Result.back() * herm(CollapseBasis(Result.back().Basis2()));
   optimize(Result);
   return Result;
}

BasicFiniteMPO
prod(BasicFiniteMPO const& x, BasicFiniteMPO const& y)
{
   return x*y;
}

BasicFiniteMPO
prod(BasicFiniteMPO const& x, BasicFiniteMPO const& y, QuantumNumbers::QuantumNumber const& q)
{
   return project(x*y, q);
}

BasicFiniteMPO
project(BasicFiniteMPO const& x, QuantumNumbers::QuantumNumber const& q)
{
   BasicFiniteMPO Result(x);
   Result.front() = ProjectBasis(Result.front().Basis1(), q) * Result.front();
   optimize(Result);
   return Result;
}


BasicFiniteMPO dot(BasicFiniteMPO const& x, BasicFiniteMPO const& y)
{
   if (x.empty())
      return x;
   if (y.empty())
      return y;

   // If the MPO's are reducible then sum over every combination that
   // leads to a scalar.  I'm not sure if this is physically useful?
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());
   BasicFiniteMPO Result;
   for (BasisList::const_iterator xI = x.Basis1().begin(); xI != x.Basis1().end(); ++xI)
   {
      for (BasisList::const_iterator yI = y.Basis1().begin(); yI != y.Basis1().end(); ++yI)
      {
         if (is_transform_target(*xI, *yI, Ident))
         {
            BasicFiniteMPO Temp = prod(project(x, *xI), project(y, *yI), QuantumNumber(x.GetSymmetryList()));
            Temp *= std::sqrt(double(degree(*xI)));
            Result += Temp;
         }
      }
   }
   optimize(Result);
   return Result;
}

BasicFiniteMPO
cross(BasicFiniteMPO const& x, BasicFiniteMPO const& y)
{
   CHECK(cross_product_exists(x.TransformsAs(), y.TransformsAs()))
      ("Cross product does not exist for these operators")
      (x.TransformsAs())(y.TransformsAs());

   return cross_product_factor(x.TransformsAs(), y.TransformsAs())
      * prod(x, y, cross_product_transforms_as(x.TransformsAs(), y.TransformsAs()));
}

double outer_coefficient(int dx, int dy, int dq)
{
   // The coefficient here is set so that the outer product of spin operators
   // matches the form required when calculating
   // ABCD = inner(A,B)*inner(C,D) + inner(cross(A,B), cross(C,D)) + inner(outer(A,B), outer(C,D))
   // for spin-1 operators A,B,C,D.
   // For spin-1 the coefficient is outer(A,B) = sqrt(6/5) * prod(A,B,2)

   // For spin-1/2 operators, we want outer(CH,C) to be the spin operator itself.
   // sqrt(1/2)

   double Factor = std::sqrt(double(dx+dy) / double(dq));
   if (dx == 2 && dy == 2)
   {
      Factor = std::sqrt(0.5);
   }
   else if (dx == 3 && dy == 3)
   {
      Factor = std::sqrt(6.0 / 5.0);
   }
   return Factor;
}

BasicFiniteMPO outer(BasicFiniteMPO const& x, BasicFiniteMPO const& y)
{
   QuantumNumbers::QuantumNumberList L = transform_targets(x.TransformsAs(), y.TransformsAs());
   QuantumNumbers::QuantumNumber q = L[0];
   bool Unique = true;
   for (unsigned i = 1; i < L.size(); ++i)
   {
      if (degree(L[i]) > degree(q))
      {
         q = L[i];
         Unique = true;
      }
      else if (degree(L[i]) == degree(q))
      {
         Unique = false;
      }
   }
   CHECK(Unique)("outer product is not defined for these operators")
      (x.TransformsAs())(y.TransformsAs());

   int dx = degree(x.TransformsAs());
   int dy = degree(y.TransformsAs());
   int dq = degree(q);
   return outer_coefficient(dx, dy, dq) * prod(x,y,q);
}

BasicFiniteMPO conj(BasicFiniteMPO const& x)
{
   BasicFiniteMPO Result(x);
   for (int i = 0; i < x.size(); ++i)
   {
      Result[i] = conj(Result[i]);
   }
   return Result;
}

SimpleRedOperator coarse_grain(BasicFiniteMPO const& x)
{
   if (x.is_null())
      return SimpleRedOperator();
   // else

   BasicFiniteMPO::const_iterator I = x.begin();
   OperatorComponent Op = copy(*I);
   ++I;
   while (I != x.end())
   {
      Op = local_tensor_prod(Op, *I);
      ++I;
   }

   // Sum over the components to get the final result
   SimpleRedOperator Result(Op.LocalBasis1(), Op.LocalBasis2());
   for (unsigned n = 0; n < Op.Basis1().size(); ++n)
   {
      Result += Op(n,0);
   }
   return Result;
}

BasicFiniteMPO fine_grain(SimpleOperator const& x,
                     std::vector<BasisList> const& LocalBasis1,
                     std::vector<BasisList> const& LocalBasis2)
{
   CHECK_EQUAL(LocalBasis1.size(), LocalBasis2.size());
   // quick return if we're already fine enough
   if (LocalBasis1.size() == 0)
   {
      return BasicFiniteMPO();
   }
   if (LocalBasis1.size() == 1)
   {
      // Turn the operator into a 1x1 MPO
      BasicFiniteMPO Result(1);
      BasisList Vacuum = make_vacuum_basis(x.GetSymmetryList());
      BasisList B1(x.GetSymmetryList());
      B1.push_back(x.TransformsAs());
      Result[0] = OperatorComponent(x.Basis1(), x.Basis2(), B1, Vacuum);
      Result[0].insert(0,0, copy(x));
      return Result;
   }

   // To invert the decomposition, we need to obtain the sequence of ProductBasis and
   // iteratively apply decompose_tensor_prod.
   // The order of the tensor products is (((1*2)*3)*4)*5
   // so we need to reverse this order
   std::stack<ProductBasis<BasisList, BasisList> > TensorProdBasis1;
   std::stack<ProductBasis<BasisList, BasisList> > TensorProdBasis2;
   TensorProdBasis1.push(make_product_basis(LocalBasis1[0], LocalBasis1[1]));
   TensorProdBasis2.push(make_product_basis(LocalBasis2[0], LocalBasis2[1]));
   for (unsigned i = 2; i < LocalBasis1.size(); ++i)
   {
      TensorProdBasis1.push(make_product_basis(TensorProdBasis1.top().Basis(), LocalBasis1[i]));
      TensorProdBasis2.push(make_product_basis(TensorProdBasis2.top().Basis(), LocalBasis2[i]));
   }

   // And now reverse the decomposition
   BasicFiniteMPO Result(LocalBasis1.size());
   OperatorComponent R1, R2;
   std::tie(R1, R2) = decompose_local_tensor_prod(x, TensorProdBasis1.top(), TensorProdBasis2.top());
   int i = LocalBasis1.size()-1;
   Result[i] = std::move(R2);
   --i;
   TensorProdBasis1.pop();
   TensorProdBasis2.pop();
   while (!TensorProdBasis1.empty())
   {
      std::tie(R1, R2) = decompose_local_tensor_prod(R1, TensorProdBasis1.top(), TensorProdBasis2.top());
      Result[i] = std::move(R2);
      --i;
      TensorProdBasis1.pop();
      TensorProdBasis2.pop();
   }
   CHECK(i == 0);
   CHECK(TensorProdBasis2.empty());
   Result[0] = std::move(R1);

   optimize(Result);
   return Result;
}

BasicFiniteMPO exp(BasicFiniteMPO const& x)
{
   CHECK(x.is_scalar())("Must be a scalar operator to calculate the operator exponential!");

   SimpleOperator Op = copy(coarse_grain(x).scalar());
   Op = Tensor::exp(Op);
   BasicFiniteMPO Result = fine_grain(Op, x.LocalBasis1List(), x.LocalBasis2List());
   return Result;
}

BasicFiniteMPO
pow(BasicFiniteMPO const& x, int n)
{
   if (n == 0)
   {
      return MakeIdentityFrom(x);
   }
   else if (n%2 == 0)
   {
      return pow(x*x, n/2);
   }
   else if (n == 1)
   {
      return x;
   }
   else
   {
      return x*pow(x*x, (n-1)/2);
   }
}

BasicFiniteMPO
adjoint(BasicFiniteMPO const& x)
{
   BasicFiniteMPO Result(x);
   for (BasicFiniteMPO::iterator I = Result.begin(); I != Result.end(); ++I)
   {
      *I = adjoint(*I);
   }
   return Result;
}

BasicFiniteMPO
inv_adjoint(BasicFiniteMPO const& x)
{
   BasicFiniteMPO Result(x);
   for (BasicFiniteMPO::iterator I = Result.begin(); I != Result.end(); ++I)
   {
      *I = inv_adjoint(*I);
   }
   return Result;
}

BasicFiniteMPO
BasicFiniteMPO::make_identity(std::vector<BasisList> const& Basis)
{
   BasicFiniteMPO Result(Basis.size());
   if (Basis.empty())
      return Result;

   QuantumNumbers::QuantumNumber Ident(Basis[0].GetSymmetryList());
   BasisList b(Basis[0].GetSymmetryList());
   b.push_back(Ident);
   for (unsigned i = 0; i < Basis.size(); ++i)
   {
      Result[i] = OperatorComponent(Basis[i], Basis[i], b, b);
      Result[i].insert(0,0, SimpleOperator::make_identity(Basis[i]));
   }
   return Result;
}

BasicFiniteMPO
BasicFiniteMPO::make_identity(std::vector<BasisList> const& Basis, QuantumNumber const& q)
{
   BasicFiniteMPO Result(Basis.size());
   if (Basis.empty())
      return Result;
   BasisList b(Basis[0].GetSymmetryList());
   b.push_back(q);
   for (unsigned i = 0; i < Basis.size(); ++i)
   {
      Result[i] = OperatorComponent(Basis[i], Basis[i], b, b);
      Result[i].insert(0,0, SimpleOperator::make_identity(Basis[i]));
   }
   return Result;
}

BasicFiniteMPO
MakeIdentityFrom(BasicFiniteMPO const& x)
{
   BasicFiniteMPO Result(x.size());
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());
   BasisList b(x.GetSymmetryList());
   b.push_back(Ident);
   for (int i = 0; i < x.size(); ++i)
   {
      CHECK_EQUAL(x.LocalBasis1(i), x.LocalBasis2(i));
      Result[i] = OperatorComponent(x.LocalBasis1(i), x.LocalBasis1(i), b, b);
      Result[i].insert(0,0, SimpleOperator::make_identity(x.LocalBasis1(i)));
   }
   return Result;
}

BasicFiniteMPO
MakeIdentityFrom(BasicFiniteMPO const& x, QuantumNumber const& q)
{
   BasicFiniteMPO Result(x.size());
   BasisList b(x.GetSymmetryList());
   b.push_back(q);
   for (int i = 0; i < x.size(); ++i)
   {
      CHECK_EQUAL(x.LocalBasis1(i), x.LocalBasis2(i));
      Result[i] = OperatorComponent(x.LocalBasis1(i), x.LocalBasis1(i), b, b);
      Result[i].insert(0,0, SimpleOperator::make_identity(x.LocalBasis1(i)));
   }
   return Result;
}

BasicFiniteMPO identity_mpo(SiteListType const& SiteList, QuantumNumbers::QuantumNumber const& q)
{
   BasicFiniteMPO Result(SiteList.size());
   BasisList b = make_single_basis(q);
   for (unsigned i = 0; i < SiteList.size(); ++i)
   {
      Result[i] = OperatorComponent(SiteList[i].Basis1(), b, b);
      Result[i].insert(0,0, SiteList[i].identity());
   }
   return Result;
}

BasicFiniteMPO identity_mpo(SiteListType const& SiteList)
{
   if (SiteList.empty())
      return BasicFiniteMPO();
   return identity_mpo(SiteList, QuantumNumbers::QuantumNumber(SiteList[0].GetSymmetryList()));
}

BasicFiniteMPO string_mpo(SiteListType const& SiteList,
			  std::string const& OpName, QuantumNumbers::QuantumNumber const& Trans)
{
   BasicFiniteMPO Result(SiteList.size());

   BasisList Vacuum = make_single_basis(Trans);

   // Assemble the JW-string
   for (unsigned i = 0; i < SiteList.size(); ++i)
   {
      if (!SiteList[i].operator_exists(OpName))
      {
         WARNING("JW-string operator doesn't exist at a lattice site, using the identity")(i)(OpName);
      }
      SimpleOperator Op = SiteList[i].operator_exists(OpName) ? copy(SiteList[i][OpName])
	 : SiteList[i].identity();
      Result[i] = OperatorComponent(Op.Basis1(), Op.Basis2(), Vacuum, Vacuum);
      Result[i].insert(0,0, std::move(Op));
   }
   return Result;
}

BasicFiniteMPO string_mpo(SiteListType const& SiteList, std::string const& OpName)
{
   if (SiteList.empty())
      return BasicFiniteMPO();
   return string_mpo(SiteList, OpName, QuantumNumbers::QuantumNumber(SiteList[0].GetSymmetryList()));
}

BasicFiniteMPO
ParseStringOperator(SiteListType const& SiteList, std::string const& Expr, int Size)
{
   CHECK(Size % SiteList.size() == 0)
      ("The size of the string operator must be a multiple of the unit cell");

   BasicFiniteMPO Result(SiteList.size());

   BasisList Vacuum = make_vacuum_basis(SiteList[0].GetSymmetryList());

   for (unsigned i = 0; i < SiteList.size(); ++i)
   {
      SiteOperator Op = ParseSiteOperator(SiteList[i], Expr);
      CHECK_EQUAL(Op.Commute(), LatticeCommute::Bosonic)("String operator must have bosonic commutation");
      Result[i] = OperatorComponent(Op.Basis1(), Op.Basis2(), Vacuum, Vacuum);
      Result[i].insert(0,0, std::move(Op));
   }
   Result = repeat(Result, Size / SiteList.size());

   return Result;
}

double
log_norm_frob_sq(BasicFiniteMPO const& Op)
{
   return 0;
}
