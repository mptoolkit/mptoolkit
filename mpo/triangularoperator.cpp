// -*- C++ -*- $Id$

#include "triangularoperator.h"

MpOpTriangular::MpOpTriangular(OperatorComponent const& Data)
   : Data_(Data)
{
}

std::ostream&
operator<<(std::ostream& out, MpOpTriangular const& op)
{
   return out << op.data();
}

MpOpTriangular TriangularOneSite(SiteOperator const& x)
{
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());
   BasisList b(x.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(x.TransformsAs());

   SimpleOperator I = SimpleOperator::make_identity(x.Basis1().Basis());

   OperatorComponent data(x.Basis1().Basis(), b, b);
   data.set_operator(0,0,I);
   data.set_operator(1,1,I);
   data.set_operator(1,0,x);

   return MpOpTriangular(data);
}

MpOpTriangular TriangularTwoSite(SiteOperator const& x, SiteOperator const& y, QuantumNumber const& Trans)
{
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());
   BasisList b(x.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(y.TransformsAs());
   b.push_back(Trans);

   SimpleOperator I = SimpleOperator::make_identity(x.Basis1().Basis());

   OperatorComponent data(x.Basis1().Basis(), b, b);
   data.set_operator(0,0,I);
   data.set_operator(2,2,I);
   data.set_operator(1,0,y);
   data.set_operator(2,1,x);

   return MpOpTriangular(data);
}

MpOpTriangular TriangularTwoSite(SiteOperator const& x, SiteOperator const& y)
{
   QuantumNumber q;
   if (num_transform_targets(x.TransformsAs(), y.TransformsAs()) == 1)
      q = transform_targets(x.TransformsAs(), y.TransformsAs())[0];
   else
      q = QuantumNumber(x.GetSymmetryList());
   return TriangularTwoSite(x, y, q);
}

MpOpTriangular TriangularThreeSite(SiteOperator const& x, SiteOperator const& y, SiteOperator const& z, 
				   QuantumNumber const& yz_trans)
{
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());

   BasisList b(x.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(z.TransformsAs());
   b.push_back(yz_trans);
   b.push_back(Ident);

   SimpleOperator I = SimpleOperator::make_identity(x.Basis1().Basis());

   OperatorComponent data(x.Basis1().Basis(), b, b);
   data.set_operator(0,0,I);
   data.set_operator(3,3,I);
   data.set_operator(1,0,z);
   data.set_operator(2,1,y);
   data.set_operator(3,2,x);

   return MpOpTriangular(data);
}

MpOpTriangular TriangularStretchedTwoSite(SiteOperator const& x, int NumNeighbors,
					  SiteOperator const& y)
{
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());
   BasisList b(x.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(y.TransformsAs());
   for (int i = 0; i < NumNeighbors; ++i)
   {
      b.push_back(y.TransformsAs());
   }
   b.push_back(Ident);

   SimpleOperator I = SimpleOperator::make_identity(x.Basis1().Basis());

   OperatorComponent data(x.Basis1().Basis(), b, b);
   data.set_operator(0,0,I);
   data.set_operator(1,0,y);
   for (int i = 0; i < NumNeighbors; ++i)
   {
      data.set_operator(i+2, i+1, I);
   }
   data.set_operator(2+NumNeighbors,1+NumNeighbors,x);
   data.set_operator(2+NumNeighbors,2+NumNeighbors,I);

   return MpOpTriangular(data);
}


MpOpTriangular TriangularThreeSite(SiteOperator const& x, SiteOperator const& y, SiteOperator const& z)
{
   QuantumNumbers::QuantumNumberList ql = transform_targets(y.TransformsAs(), z.TransformsAs());
   CHECK_EQUAL(ql.size(), 1)("ambiguous coupling of y and z")(y.TransformsAs())(z.TransformsAs());
   return TriangularThreeSite(x,y,z, ql[0]);
}

MpOpTriangular TriangularFourSite(SiteOperator const& w, SiteOperator const& x, 
                                  SiteOperator const& y, SiteOperator const& z, 
                                  QuantumNumber const& xyz_trans,
                                  QuantumNumber const& yz_trans)
{
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());

   BasisList b(x.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(z.TransformsAs());
   b.push_back(yz_trans);
   b.push_back(xyz_trans);
   b.push_back(Ident);

   SimpleOperator I = SimpleOperator::make_identity(x.Basis1().Basis());

   OperatorComponent data(x.Basis1().Basis(), b, b);
   data.set_operator(0,0,I);
   data.set_operator(4,4,I);
   data.set_operator(1,0,z);
   data.set_operator(2,1,y);
   data.set_operator(3,2,x);
   data.set_operator(4,3,w);

   return MpOpTriangular(data);
}

MpOpTriangular TriangularFourSite(SiteOperator const& w, SiteOperator const& x, 
                                  SiteOperator const& y, SiteOperator const& z)
{
   QuantumNumbers::QuantumNumberList ql = transform_targets(y.TransformsAs(), z.TransformsAs());
   CHECK_EQUAL(ql.size(), 1)("ambiguous coupling of y and z")(y.TransformsAs())(z.TransformsAs());
   QuantumNumber yz = ql[0];

   ql = transform_targets(x.TransformsAs(), yz);
   CHECK_EQUAL(ql.size(), 1)("ambiguous coupling of x and yz")(x.TransformsAs())
      (y.TransformsAs())(z.TransformsAs());
   QuantumNumber xyz = ql[0];

   return TriangularFourSite(w, x, y, z, xyz, yz);
}


MpOpTriangular ZigZagChain(SiteOperator const& S, SiteOperator const& T, double J1, double J2)
{
   QuantumNumbers::QuantumNumber Ident(S.GetSymmetryList());

   BasisList b(S.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(S.TransformsAs());
   b.push_back(S.TransformsAs());
   b.push_back(Ident);

   SimpleOperator I = SimpleOperator::make_identity(S.Basis1().Basis());

   OperatorComponent data(S.Basis1().Basis(), b, b);
   data.set_operator(0,0,I);
   data.set_operator(3,3,I);
   data.set_operator(1,0,S);
   data.set_operator(2,1,I);
   data.set_operator(3,2,J2*T);
   data.set_operator(3,1,J1*T);

   return MpOpTriangular(data);
}

// temporary hacks for specific hamiltonians

MpOpTriangular TriangularTwoSitePBC(SiteOperator const& x, SiteOperator const& y, QuantumNumber const& Trans)
{
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());
   BasisList b(x.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(y.TransformsAs());
   b.push_back(y.TransformsAs());
   b.push_back(Trans);

   SimpleOperator I = SimpleOperator::make_identity(x.Basis1().Basis());

   OperatorComponent data(x.Basis1().Basis(), b, b);
   data.set_operator(0,0,I);
   data.set_operator(2,2,I);
   data.set_operator(3,3,I);
   data.set_operator(1,0,y);
   data.set_operator(3,1,x);

   return MpOpTriangular(data);
}

MpOpTriangular TriangularTwoSitePBC(SiteOperator const& x, SiteOperator const& y)
{
   CHECK_EQUAL(num_transform_targets(x.TransformsAs(), y.TransformsAs()), 1);
   return TriangularTwoSitePBC(x, y, transform_targets(x.TransformsAs(), y.TransformsAs())[0]);
}

MpOpTriangular TriangularTwoSitePBC_Boundary(SiteOperator const& x, SiteOperator const& y, QuantumNumber const& Trans)
{
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());
   BasisList b(x.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(y.TransformsAs());
   b.push_back(y.TransformsAs());
   b.push_back(Trans);

   SimpleOperator I = SimpleOperator::make_identity(x.Basis1().Basis());

   OperatorComponent data(x.Basis1().Basis(), b, b);
   data.set_operator(0,0,I);
   data.set_operator(3,3,I);
   data.set_operator(1,0,y);
   data.set_operator(2,0,y);
   data.set_operator(3,1,x);
   data.set_operator(3,2,x);

   return MpOpTriangular(data);
}

MpOpTriangular TriangularTwoSitePBC_Boundary(SiteOperator const& x, SiteOperator const& y)
{
   CHECK_EQUAL(num_transform_targets(x.TransformsAs(), y.TransformsAs()), 1);
   return TriangularTwoSitePBC_Boundary(x, y, transform_targets(x.TransformsAs(), y.TransformsAs())[0]);
}

// arithmetic

MpOpTriangular& operator*=(MpOpTriangular& Op, double x)
{
   // Multiply the first column by x
   SimpleOperator T = SimpleOperator::make_identity(Op.data().Basis1());
   T(0,0) = x;
   Op.data() = prod(Op.data(), T);
   // except for the (0,0) entry which is always the identity
   Op.data().set_operator(0,0, SimpleOperator::make_identity(Op.data().LocalBasis1()));
   return Op;
}

MpOpTriangular& operator*=(MpOpTriangular& Op, std::complex<double> x)
{
   // Multiply the first column by x
   SimpleOperator T = SimpleOperator::make_identity(Op.data().Basis1());
   T(0,0) = x;
   Op.data() = prod(Op.data(), T);
   // except for the (0,0) entry which is always the identity
   Op.data().set_operator(0,0, SimpleOperator::make_identity(Op.data().LocalBasis1()));
   return Op;
}

MpOpTriangular operator*(MpOpTriangular const& Op, double x)
{
   MpOpTriangular Result(Op);
   Result *= x;
   return Result;
}

MpOpTriangular operator*(double x, MpOpTriangular const& Op)
{
   MpOpTriangular Result(Op);
   Result *= x;
   return Result;
}

MpOpTriangular operator*(MpOpTriangular const& Op, std::complex<double> x)
{
   MpOpTriangular Result(Op);
   Result *= x;
   return Result;
}

MpOpTriangular operator*(std::complex<double> x, MpOpTriangular const& Op)
{
   MpOpTriangular Result(Op);
   Result *= x;
   return Result;
}

MpOpTriangular operator+(MpOpTriangular const& x, MpOpTriangular const& y)
{
   PRECONDITION_EQUAL(x.GetSymmetryList(), y.GetSymmetryList());
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());

   int xSz = x.data().Basis1().size();
   int ySz = y.data().Basis1().size();

   // This is a somewhat simplistic approach.  We just combine the first and last
   // rows/columns, and join the remainder.

   BasisList b(x.GetSymmetryList());
   for (int i = 0; i < xSz-1; ++i)
      b.push_back(x.data().Basis1()[i]);
   for (int i = 1; i < ySz; ++i)
      b.push_back(y.data().Basis1()[i]);

   SimpleOperator xProjector(b, x.data().Basis1());
   for (int i = 0; i < xSz-1; ++i)         // include the first row in x only
      xProjector(i,i) = 1.0;
   xProjector(xSz+ySz-3, xSz-1) = 1.0;   // last row is special

   SimpleOperator yProjector(b, y.data().Basis1());
   yProjector(0, 0) = 1.0;
   for (int i = 1; i < ySz; ++i)
      yProjector(i+xSz-2,i) = 1.0;

   OperatorComponent Result = ::prod(::prod(xProjector, x.data()), herm(xProjector))
      + ::prod(::prod(yProjector, y.data()), herm(yProjector));

   SimpleOperator I = SimpleOperator::make_identity(x.data().LocalBasis1());
   Result.set_operator(0,0, I);
   Result.set_operator(xSz+ySz-3,xSz+ySz-3, I);

   return MpOpTriangular(Result);
}

MpOpTriangular& operator+=(MpOpTriangular& x, MpOpTriangular const& y)
{
   x = x+y;
   return x;
}

MpOpTriangular operator-(MpOpTriangular const& x, MpOpTriangular const& y)
{
   return x + (-1.0)*y;
}

MpOpTriangular operator*(MpOpTriangular const& x, MpOpTriangular const& y)
{
   Tensor::ProductBasis<BasisList, BasisList> B(x.Basis1(), y.Basis1());
   OperatorComponent Result(x.LocalBasis1(), y.LocalBasis2(), B.Basis(), B.Basis());
   Result.data() = direct_product(x.data().data(), y.data().data(), 
                               LinearAlgebra::Multiplication<SimpleRedOperator, SimpleRedOperator>());
   return MpOpTriangular(Result);
}
			 
// initial operators

StateComponent Initial_E(MpOpTriangular const& m, VectorBasis const& Vac)
{
   StateComponent Result(m.data().Basis1(), Vac, Vac);
   Result[m.data().Basis1().size()-1] = MatrixOperator::make_identity(Vac);
   return Result;
}

StateComponent Initial_F(MpOpTriangular const& m, VectorBasis const& Vac)
{
   StateComponent Result(m.data().Basis2(), Vac, Vac);
   Result[0] = MatrixOperator::make_identity(Vac);
   return Result;
}

StateComponent Initial_E(MpOpTriangular const& m)
{
   VectorBasis Vac(make_vacuum_basis(m.GetSymmetryList()));
   return Initial_E(m, Vac);
}

StateComponent Initial_F(MpOpTriangular const& m)
{
   VectorBasis Vac(make_vacuum_basis(m.GetSymmetryList()));
   return Initial_F(m, Vac);
}

MpOpTriangular Section(MpOpTriangular const& x, int Offset, int n)
{
   BasisList NewB(x.GetSymmetryList());
   for (int i = 0; i < n; ++i)
      NewB.push_back(x.Basis1()[i+Offset]);

   SimpleOperator P(NewB, x.Basis1());
   for (int i = 0; i < n; ++i)
      P(i,i+Offset) = 1.0;

   return MpOpTriangular(prod(P, prod(x.data(), herm(P))));
}

MpOpTriangular SectionTopLeft(MpOpTriangular const& x, int n)
{
   return Section(x, 0, n);
}

MpOpTriangular SectionBottomRight(MpOpTriangular const& x, int n)
{
   return Section(x, x.Basis1().size() - n, n);
}

StateComponent Initial_E(OperatorComponent const& m, VectorBasis const& Vac)
{
   StateComponent Result(m.Basis1(), Vac, Vac);
   Result[m.Basis1().size()-1] = MatrixOperator::make_identity(Vac);
   return Result;
}

StateComponent Initial_F(OperatorComponent const& m, VectorBasis const& Vac)
{
   StateComponent Result(m.Basis2(), Vac, Vac);
   Result[0] = MatrixOperator::make_identity(Vac);
   return Result;
}

std::vector<BasisList>
ExtractLocalBasis(SimpleMPOperator const& Op)
{
   std::vector<BasisList> Result;
   Result.reserve(Op.size());
   for (SimpleMPOperator::const_iterator I = Op.begin(); I != Op.end(); ++I)
   {
      Result.push_back(I->LocalBasis1());
   }
   return Result;
}


