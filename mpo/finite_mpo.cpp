// -*- C++ -*- $Id$

#include "finite_mpo.h"
#include "pstream/variant.h"
#include "tensor/tensor_exponential.h"
#include "lattice/siteoperator-parser.h"
#include <algorithm>
#include <stack>

using QuantumNumbers::QuantumNumber;

bool FiniteMPO::is_irreducible() const
{
   return this->is_null() || (this->Basis1().size() == 1 && QuantumNumbers::is_scalar(this->qn2()));
}

bool FiniteMPO::is_scalar() const
{
   return this->is_null() || (this->Basis1().size() == 1 && QuantumNumbers::is_scalar(this->qn1()) 
			      && QuantumNumbers::is_scalar(this->qn2()));
}

#if !defined(DISABLE_FINITE_MPO_TRANSFORMS_AS)
QuantumNumbers::QuantumNumber 
FiniteMPO::TransformsAs() const
{
   if (this->is_null())
      return QuantumNumbers::QuantumNumber();
   CHECK(this->is_irreducible());
   return this->Basis1()[0];
}
#endif

QuantumNumbers::QuantumNumber 
FiniteMPO::qn1() const
{
   if (this->is_null())
      return QuantumNumbers::QuantumNumber();
   CHECK_EQUAL(this->Basis1().size(), 1);
   return this->Basis1()[0];
}

QuantumNumbers::QuantumNumber 
FiniteMPO::qn2() const
{
   if (this->is_null())
      return QuantumNumbers::QuantumNumber();
   DEBUG_CHECK_EQUAL(this->Basis2().size(), 1);  // DEBUG since this should always be true
   return this->Basis2()[0];
}

PStream::opstream& operator<<(PStream::opstream& out, FiniteMPO const& op)
{
   return out << op.data();
}

PStream::ipstream& operator>>(PStream::ipstream& in, FiniteMPO& op)
{
   return in >> op.data();
}

void print_structure(FiniteMPO const& Op, std::ostream& out)
{
   out << "FiniteMPO has " << Op.size() << " sites\n";
   for (int i = 0; i < Op.size(); ++i)
   {
      out << "Site " << i << " dimension " << Op[i].size1() << " x " << Op[i].size2() << '\n';
      print_structure(Op[i], out);
   }
}

FiniteMPO
join(FiniteMPO const& Op1, FiniteMPO const& Op2)
{
   if (Op1.is_null())
      return Op2;
   if (Op2.is_null())
      return Op1;
   CHECK_EQUAL(Op1.Basis2(), Op2.Basis1());
   FiniteMPO Result(Op1.size() + Op2.size());
   for (int i = 0; i < Op1.size(); ++i)
   {
      Result[i] = Op1[i];
   }
   for (int i = 0; i < Op2.size(); ++i)
   {
      Result[i+Op1.size()] = Op2[i];
   }
   return Result;
}

FiniteMPO
repeat(FiniteMPO const& Op, int Count)
{
   CHECK_EQUAL(Op.Basis1(), Op.Basis2());
   FiniteMPO Result(Op.size()*Count);
   for (int i = 0; i < Result.size(); ++i)
   {
      Result[i] = Op[i%Op.size()];
   }
   return Result;
}

void optimize(FiniteMPO& Op)
{
   if (Op.size() < 2)
      return;

   SimpleRedOperator X = coarse_grain(Op);

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

   SimpleRedOperator Y = coarse_grain(Op);
   CHECK(norm_frob(X-Y) < 1E-10*Y.Basis1().total_degree())("failure in MPO optimize()!")(X)(Y);
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

FiniteMPO
operator-(FiniteMPO const& x)
{
   FiniteMPO Result(x);
   if (!Result.empty())
      Result.front() *= -1.0;
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
   CHECK_EQUAL(x.size(), y.size());
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
prod(FiniteMPO const& x, FiniteMPO const& y)
{
   return x*y;
}

FiniteMPO
prod(FiniteMPO const& x, FiniteMPO const& y, QuantumNumbers::QuantumNumber const& q)
{
   return project(x*y, q);
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
            FiniteMPO Temp = prod(project(x, *xI), project(y, *yI), QuantumNumber(x.GetSymmetryList()));
            Temp *= std::sqrt(double(degree(*xI)));
            Result += Temp;
         }
      }
   }
   return Result;
}

FiniteMPO
cross(FiniteMPO const& x, FiniteMPO const& y)
{
   // Legitimate uses of FiniteMPO::TransformsAs()
#if !defined(DISABLE_FINITE_MPO_TRANSFORMS_AS)
   CHECK(cross_product_exists(x.TransformsAs(), y.TransformsAs()))
      ("Cross product does not exist for these operators")
      (x.TransformsAs())(y.TransformsAs());

   return cross_product_factor(x.TransformsAs(), y.TransformsAs())
      * prod(x, y, cross_product_transforms_as(x.TransformsAs(), y.TransformsAs()));
#else
   return x;
#endif
}

FiniteMPO outer(FiniteMPO const& x, FiniteMPO const& y)
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
   return std::sqrt(double(degree(x.TransformsAs()) + degree(y.TransformsAs())) / degree(q)) * prod(x,y,q);
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

SimpleRedOperator coarse_grain(FiniteMPO const& x)
{
   if (x.is_null())
      return SimpleRedOperator();
   // else

   FiniteMPO::const_iterator I = x.begin();
   OperatorComponent Op = *I;
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

FiniteMPO fine_grain(SimpleOperator const& x,
		     std::vector<BasisList> const& LocalBasis1,
		     std::vector<BasisList> const& LocalBasis2)
{
   CHECK_EQUAL(LocalBasis1.size(), LocalBasis2.size());
   // quick return if we're already fine enough
   if (LocalBasis1.size() == 0)
   {
      return FiniteMPO();
   }
   if (LocalBasis1.size() == 1)
   {
      // Turn the operator into a 1x1 MPO
      FiniteMPO Result(1);
      BasisList Vacuum = make_vacuum_basis(x.GetSymmetryList());
      BasisList B1(x.GetSymmetryList());
      B1.push_back(x.TransformsAs());
      Result[0] = OperatorComponent(x.Basis1(), x.Basis2(), B1, Vacuum);
      Result[0](0,0) = x;
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
   FiniteMPO Result(LocalBasis1.size());
   OperatorComponent R1, R2;
   boost::tie(R1, R2) = decompose_local_tensor_prod(x, TensorProdBasis1.top(), TensorProdBasis2.top());
   int i = LocalBasis1.size()-1;
   Result[i] = R2;
   --i;
   TensorProdBasis1.pop();
   TensorProdBasis2.pop();
   while (!TensorProdBasis1.empty())
   {
      boost::tie(R1, R2) = decompose_local_tensor_prod(R1, TensorProdBasis1.top(), TensorProdBasis2.top());
      Result[i] = R2;
      --i;
      TensorProdBasis1.pop();
      TensorProdBasis2.pop();
   }
   CHECK(i == 0);
   CHECK(TensorProdBasis2.empty());
   Result[0] = R1;
   return Result;
}

FiniteMPO exp(FiniteMPO const& x)
{
   CHECK(x.is_scalar())("Must be a scalar operator to calculate the operator exponential!");

   SimpleOperator Op = coarse_grain(x).scalar();
   Op = Tensor::Exponentiate(Op);
   FiniteMPO Result = fine_grain(Op, x.LocalBasis1List(), x.LocalBasis2List());
   return Result;
}

FiniteMPO
pow(FiniteMPO const& x, int n)
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

FiniteMPO
adjoint(FiniteMPO const& x)
{
   FiniteMPO Result(x);
   for (FiniteMPO::iterator I = Result.begin(); I != Result.end(); ++I)
   {
      *I = adjoint(*I);
   }
   return Result;
}

FiniteMPO
inv_adjoint(FiniteMPO const& x)
{
   FiniteMPO Result(x);
   for (FiniteMPO::iterator I = Result.begin(); I != Result.end(); ++I)
   {
      *I = inv_adjoint(*I);
   }
   return Result;
}

FiniteMPO
FiniteMPO::make_identity(std::vector<BasisList> const& Basis)
{
   FiniteMPO Result(Basis.size());
   if (Basis.empty())
      return Result;

   QuantumNumbers::QuantumNumber Ident(Basis[0].GetSymmetryList());
   BasisList b(Basis[0].GetSymmetryList());
   b.push_back(Ident);
   for (unsigned i = 0; i < Basis.size(); ++i)
   {
      Result[i] = OperatorComponent(Basis[i], Basis[i], b, b);
      Result[i](0,0) = SimpleOperator::make_identity(Basis[i]);
   }
   return Result;
}

FiniteMPO
FiniteMPO::make_identity(std::vector<BasisList> const& Basis, QuantumNumber const& q)
{
   FiniteMPO Result(Basis.size());
   if (Basis.empty())
      return Result;
   BasisList b(Basis[0].GetSymmetryList());
   b.push_back(q);
   for (unsigned i = 0; i < Basis.size(); ++i)
   {
      Result[i] = OperatorComponent(Basis[i], Basis[i], b, b);
      Result[i](0,0) = SimpleOperator::make_identity(Basis[i]);
   }
   return Result;
}

FiniteMPO
MakeIdentityFrom(FiniteMPO const& x)
{
   FiniteMPO Result(x.size());
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());
   BasisList b(x.GetSymmetryList());
   b.push_back(Ident);
   for (int i = 0; i < x.size(); ++i)
   {
      CHECK_EQUAL(x.LocalBasis1(i), x.LocalBasis2(i));
      Result[i] = OperatorComponent(x.LocalBasis1(i), x.LocalBasis1(i), b, b);
      Result[i](0,0) = SimpleOperator::make_identity(x.LocalBasis1(i));
   }
   return Result;
}

FiniteMPO
MakeIdentityFrom(FiniteMPO const& x, QuantumNumber const& q)
{
   FiniteMPO Result(x.size());
   BasisList b(x.GetSymmetryList());
   b.push_back(q);
   for (int i = 0; i < x.size(); ++i)
   {
      CHECK_EQUAL(x.LocalBasis1(i), x.LocalBasis2(i));
      Result[i] = OperatorComponent(x.LocalBasis1(i), x.LocalBasis1(i), b, b);
      Result[i](0,0) = SimpleOperator::make_identity(x.LocalBasis1(i));
   }
   return Result;
}

FiniteMPO identity_mpo(SiteListType const& SiteList, QuantumNumbers::QuantumNumber const& q)
{
   FiniteMPO Result(SiteList.size());
   BasisList b = make_single_basis(q);
   for (unsigned i = 0; i < SiteList.size(); ++i)
   {
      Result[i] = OperatorComponent(SiteList[i].Basis1(), b, b);
      Result[i](0,0) = SiteList[i].identity();
   }
   return Result;
}

FiniteMPO identity_mpo(SiteListType const& SiteList)
{
   if (SiteList.empty())
      return FiniteMPO();
   return identity_mpo(SiteList, QuantumNumbers::QuantumNumber(SiteList[0].GetSymmetryList()));
}

FiniteMPO string_mpo(SiteListType const& SiteList,
		     std::string const& OpName, QuantumNumbers::QuantumNumber const& Trans)
{
   FiniteMPO Result(SiteList.size());

   BasisList Vacuum = make_single_basis(Trans);

   // Assemble the JW-string
   for (unsigned i = 0; i < SiteList.size(); ++i)
   {
      if (!SiteList[i].operator_exists(OpName))
      {
	 WARNING("JW-string operator doesn't exist at a lattice site, using the identity")(i)(OpName);
      }
      SimpleOperator Op = SiteList[i].operator_exists(OpName) ? SiteList[i][OpName]
	 : SiteList[i].identity();
      Result[i] = OperatorComponent(Op.Basis1(), Op.Basis2(), Vacuum, Vacuum);
      Result[i](0,0) = Op;
   }
   return Result;
}   

FiniteMPO string_mpo(SiteListType const& SiteList, std::string const& OpName)
{
   if (SiteList.empty())
      return FiniteMPO();
   return string_mpo(SiteList, OpName, QuantumNumbers::QuantumNumber(SiteList[0].GetSymmetryList()));
}

FiniteMPO
ParseStringOperator(SiteListType const& SiteList, std::string const& Expr, int Size)
{
   CHECK(Size % SiteList.size() == 0)
      ("The size of the string operator must be a multiple of the unit cell");

   FiniteMPO Result(SiteList.size());
   
   BasisList Vacuum = make_vacuum_basis(SiteList[0].GetSymmetryList());

   for (unsigned i = 0; i < SiteList.size(); ++i)
   {
      SiteOperator Op = ParseSiteOperator(SiteList[i], Expr);
      CHECK_EQUAL(Op.Commute(), LatticeCommute::Bosonic)("String operator must have bosonic commutation");
      Result[i] = OperatorComponent(Op.Basis1(), Op.Basis2(), Vacuum, Vacuum);
      Result[i](0,0) = Op;
   }
   Result = repeat(Result, Size / SiteList.size());
   
   return Result;
}   

