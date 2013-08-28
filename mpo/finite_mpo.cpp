// -*- C++ -*- $Id$

#include "finite_mpo.h"
#include "pstream/variant.h"
#include <algorithm>

using QuantumNumbers::QuantumNumber;

bool FiniteMPO::is_irreducible() const
{
   //   return (!this->is_null() && this->Basis1().size() == 1);
   return (this->is_null() || this->Basis1().size() == 1);
}

QuantumNumbers::QuantumNumber 
FiniteMPO::TransformsAs() const
{
   if (this->is_null())
      return QuantumNumbers::QuantumNumber();
   CHECK(this->is_irreducible());
   return this->Basis1()[0];
}

PStream::opstream& operator<<(PStream::opstream& out, FiniteMPO const& op)
{
   return out << op.data();
}

PStream::ipstream& operator>>(PStream::ipstream& in, FiniteMPO& op)
{
   return in >> op.data();
}

void optimize(FiniteMPO& Op)
{
   if (Op.size() < 2)
      return;

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
}

FiniteMPO&
operator+=(FiniteMPO& x, FiniteMPO const& y)
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

FiniteMPO
operator+(FiniteMPO const& x, FiniteMPO const& y)
{
   if (x.empty())
      return y;
   if (y.empty())
      return x;

   FiniteMPO Result(x.size());

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

FiniteMPO
operator-(FiniteMPO const& x, FiniteMPO const& y)
{
   if (x.empty())
      return y;
   if (y.empty())
      return x;

   FiniteMPO Result(x.size());

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

FiniteMPO&
operator-=(FiniteMPO& x, FiniteMPO const& y)
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

FiniteMPO&
operator*=(FiniteMPO& m, double x)
{
   if (!m.empty())
      m.front() *= x;
   return m;
}

FiniteMPO&
operator*=(FiniteMPO& m, std::complex<double> x)
{
   if (!m.empty())
      m.front() *= x;
   return m;
}

FiniteMPO
operator*(FiniteMPO const& m, double x)
{
   FiniteMPO Result(m);
   if (!Result.empty())
      Result.front() *= x;
   return Result;
}

FiniteMPO
operator*(FiniteMPO const& m, std::complex<double> x)
{
   FiniteMPO Result(m);
   if (!Result.empty())
      Result.front() *= x;
   return Result;
}

FiniteMPO
operator*(double x, FiniteMPO const& m)
{
   FiniteMPO Result(m);
   if (!Result.empty())
      Result.front() *= x;
   return Result;
}

FiniteMPO
operator*(std::complex<double> x, FiniteMPO const& m)
{
   FiniteMPO Result(m);
   if (!Result.empty())
      Result.front() *= x;
   return Result;
}


FiniteMPO
operator*(FiniteMPO const& x, FiniteMPO const& y)
{
   if (x.empty())
      return x;
   if (y.empty())
      return y;

   FiniteMPO Result(x.size());
   for (int i = 0; i < x.size(); ++i)
   {
      Result[i] = aux_tensor_prod(x[i], y[i]);
   }
   Result.front() = CollapseBasis(Result.front().Basis1()) * Result.front();
   Result.back() = Result.back() * herm(CollapseBasis(Result.back().Basis2()));
   optimize(Result);
   return Result;
}

FiniteMPO
project(FiniteMPO const& x, QuantumNumbers::QuantumNumber const& q)
{
   FiniteMPO Result(x);
   Result.front() = ProjectBasis(Result.front().Basis1(), q) * Result.front();
   optimize(Result);
   return Result;
}


FiniteMPO dot(FiniteMPO const& x, FiniteMPO const& y)
{
   if (x.empty())
      return x;
   if (y.empty())
      return y;

   // If the MPO's are reducible then sum over every combination that
   // leads to a scalar.  I'm not sure if this is physically useful?
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());
   FiniteMPO Result;
   for (BasisList::const_iterator xI = x.Basis1().begin(); xI != x.Basis1().end(); ++xI)
   {
      for (BasisList::const_iterator yI = y.Basis1().begin(); yI != y.Basis1().end(); ++yI)
      {
         if (is_transform_target(*xI, *yI, Ident))
         {
            FiniteMPO Temp = project(x, *xI) * project(y, *yI);
            Temp *= std::sqrt(double(degree(*xI)));
            Result += Temp;
         }
      }
   }
   return Result;
}

FiniteMPO conj(FiniteMPO const& x)
{
   FiniteMPO Result(x);
   for (int i = 0; i < x.size(); ++i)
   {
      Result[i] = conj(Result[i]);
   }
   return Result;
}


