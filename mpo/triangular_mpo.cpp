// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mpo/triangular_mpo.cpp
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

#include "triangular_mpo.h"
#include "common/statistics.h"

void
TriangularMPO::check_structure() const
{
   for (int i = 0; i < this->size(); ++i)
   {
      Data_[i].check_structure();
   }
}

PStream::opstream&
 operator<<(PStream::opstream& out, TriangularMPO const& Op)
{
   out << Op.Data_;
   return out;
}

PStream::ipstream&
 operator>>(PStream::ipstream& in, TriangularMPO& Op)
{
   in >> Op.Data_;
   return in;
}

GenericMPO extract_column(TriangularMPO const& Op, int Col)
{
   GenericMPO MPOp(Op.begin(), Op.end());

   std::set<int> Cols;
   Cols.insert(Col);
   SimpleOperator ColP = make_projector_onto(MPOp.Basis2(), Cols);
   MPOp.back() = prod(MPOp.back(), herm(ColP));

   cull_unused_elements(MPOp);
   return MPOp;
}

GenericMPO extract_lower_column(TriangularMPO const& Op, int Col)
{
   GenericMPO MPOp(Op.begin(), Op.end());

   // remove the diagonal element (if it exists), by setting the elements of the first operator
   // at (Col, x) to zero.
   for (unsigned k = 0; k < MPOp.front().Basis2().size(); ++k)
   {
      zero_element(MPOp.front().data(), Col, k);
   }

   std::set<int> Cols;
   Cols.insert(Col);
   MPOp.back() = project_columns(MPOp.back(), Cols);

   cull_unused_elements(MPOp);
   return MPOp;
}

FiniteMPO 
TriangularMPO::operator()(int Row, int Col) const
{
   GenericMPO MPOp(Data_);

   std::set<int> Rows;
   Rows.insert(Row);
   MPOp.front() = project_rows(MPOp.front(), Rows);

   std::set<int> Cols;
   Cols.insert(Col);
   MPOp.back() = project_columns(MPOp.back(), Cols);
   //   TRACE(MPOp.back());

   cull_unused_elements(MPOp);
   return FiniteMPO(MPOp);
}

std::ostream&
operator<<(std::ostream& out, TriangularMPO const& op)
{
   return out << op.data();
}

#if 1
void optimize(TriangularMPO& Op)
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
      for (int i = 1; i < Op.size(); ++i)
      {
         Op[i] = T * Op[i];
         T = TruncateBasis2(Op[i]);
         if (T.size1() != T.size2())
            Reduced = true;
      }
      Op.front() = T * Op.front();
      
      // Working right to left, optimize Basis1
      T = TruncateBasis1(Op.back());
      if (T.size1() != T.size2())
         Reduced = true;
      for (int i = Op.size()-2; i >= 0; --i)
      {
         Op[i] = Op[i] * T;
	 T = TruncateBasis1(Op[i]);
         if (T.size1() != T.size2())
            Reduced = true;
      }
      Op.back() = Op.back() * T;
   }
}

#else

void optimize(TriangularMPO& Op)
{
   if (Op.size() < 2)
      return;

   double const Eps = 1E-13;

   TRACE(Op);

   bool Reduced = true; // flag to indicate that we reduced a dimension
   // loop until we do a complete sweep with no reduction in dimensions
   bool First = true;
   bool Second = true;
   while (Reduced || Second)
   {
      Reduced = false;

      OperatorComponent Op2 = Op.front();
      if (!First && Second)
      {
	 TRACE("XXXXX");
      }
      SimpleOperator T2 = TruncateBasis2MkII(Op2, First ? 0.0 : Eps);
      TRACE(norm_frob(Op.front() - Op2*T2));

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

   TRACE(Op);

}
#endif

void balance(TriangularMPO& Op)
{
   // 
}

void print_structure(TriangularMPO const& Op, std::ostream& out, double UnityEpsilon)
{
   out << "TriangularMPO has " << Op.size() << " sites\n";
   for (int i = 0; i < Op.size(); ++i)
   {
      out << "Site " << i << " dimension " << Op[i].size1() << " x " << Op[i].size2() << '\n';
      print_structure(Op[i], out, UnityEpsilon);
   }
}

// arithmetic

TriangularMPO& operator*=(TriangularMPO& Op, double x)
{
   for (int i = 0; i < Op.size(); ++i)
   {
      for (unsigned j = 1; j < Op[i].Basis2().size(); ++j)
      {
         if (iterate_at(Op[i].data(), 0, j))
            set_element(Op[i].data(), 0, j, get_element(Op[i].data(),0,j) * x);
      }
   }
   return Op;
}

TriangularMPO& operator*=(TriangularMPO& Op, std::complex<double> x)
{
   for (int i = 0; i < Op.size(); ++i)
   {
      for (unsigned j = 1; j < Op[i].Basis2().size(); ++j)
      {
         if (iterate_at(Op[i].data(), 0, j))
            set_element(Op[i].data(), 0, j, get_element(Op[i].data(),0,j) * x);
      }
   }
   return Op;
}

TriangularMPO operator*(TriangularMPO const& Op, double x)
{
   TriangularMPO Result(Op);
   Result *= x;
   return Result;
}

TriangularMPO operator*(double x, TriangularMPO const& Op)
{
   TriangularMPO Result(Op);
   Result *= x;
   return Result;
}

TriangularMPO operator*(TriangularMPO const& Op, std::complex<double> x)
{
   TriangularMPO Result(Op);
   Result *= x;
   return Result;
}

TriangularMPO operator*(std::complex<double> x, TriangularMPO const& Op)
{
   TriangularMPO Result(Op);
   Result *= x;
   return Result;
}

TriangularMPO operator+(TriangularMPO const& x, TriangularMPO const& y)
{
   if (x.size() == 0)
      return y;
   if (y.size() == 0)
      return x;

   if (x.size() != y.size())
   {
      int NewSize = statistics::lcm(x.size(), y.size());
      return repeat(x, NewSize/x.size()) + repeat(y, NewSize/y.size());
   }

   PRECONDITION_EQUAL(x.size(), y.size());
   PRECONDITION_EQUAL(x.GetSymmetryList(), y.GetSymmetryList());
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());

   TriangularMPO Result(x.size());

   for (int Here = 0; Here < x.size(); ++Here)
   {
      int x_rows = x[Here].Basis1().size();
      int y_rows = y[Here].Basis1().size();

      int x_cols = x[Here].Basis2().size();
      int y_cols = y[Here].Basis2().size();

      CHECK(norm_frob(x[Here](0,0) - y[Here](0,0)) < 1E-10)(x[Here](0,0))(y[Here](0,0));
      CHECK(norm_frob(x[Here](x_rows-1, x_cols-1) - y[Here](y_rows-1, y_cols-1)) < 1E-10)
            (x[Here](x_rows-1, x_cols-1))(y[Here](y_rows-1, y_cols-1));

      // This is a somewhat simplistic approach.  We just combine the first and last
      // rows/columns, and join the remainder.

      BasisList NewBasis1(x.GetSymmetryList());
      for (int i = 0; i < x_rows-1; ++i)
         NewBasis1.push_back(x[Here].Basis1()[i]);
      for (int i = 1; i < y_rows; ++i)
         NewBasis1.push_back(y[Here].Basis1()[i]);

      SimpleOperator xProjector1(NewBasis1, x[Here].Basis1());
      for (int i = 0; i < x_rows-1; ++i)         // include the first row in x only
         xProjector1(i,i) = 1.0;
      xProjector1(x_rows+y_rows-3, x_rows-1) = 1.0;   // last row is special

      SimpleOperator yProjector1(NewBasis1, y[Here].Basis1());
      yProjector1(0, 0) = 1.0;  // first row is special
      for (int i = 1; i < y_rows; ++i)
         yProjector1(i+x_rows-2,i) = 1.0;

      BasisList NewBasis2(x.GetSymmetryList());
      for (int i = 0; i < x_cols-1; ++i)
         NewBasis2.push_back(x[Here].Basis2()[i]);
      for (int i = 1; i < y_cols; ++i)
         NewBasis2.push_back(y[Here].Basis2()[i]);

      SimpleOperator xProjector2(NewBasis2, x[Here].Basis2());
      for (int i = 0; i < x_cols-1; ++i)         // include the first row in x only
         xProjector2(i,i) = 1.0;
      xProjector2(x_cols+y_cols-3, x_cols-1) = 1.0;   // last row is special

      SimpleOperator yProjector2(NewBasis2, y[Here].Basis2());
      yProjector2(0, 0) = 1.0;  // first row is special
      for (int i = 1; i < y_cols; ++i)
         yProjector2(i+x_cols-2,i) = 1.0;

      OperatorComponent Next = triple_prod(xProjector1, x[Here], herm(xProjector2)) 
         + triple_prod(yProjector1, y[Here], herm(yProjector2));

      Next(0,0) = x[Here](0,0);
      Next(Next.size1()-1, Next.size2()-1) = x[Here](x[Here].size1()-1, x[Here].size2()-1);

      Result[Here] = Next;
   }

   //   TRACE(x)(y)(Result);

   optimize(Result);

   return Result;
}

TriangularMPO& operator+=(TriangularMPO& x, TriangularMPO const& y)
{
   x = x+y;
   return x;
}

TriangularMPO operator-(TriangularMPO const& x, TriangularMPO const& y)
{
   return x + (-1.0)*y;
}

TriangularMPO operator*(TriangularMPO const& x, TriangularMPO const& y)
{
   QuantumNumber qx = x.TransformsAs();
   QuantumNumber qy = y.TransformsAs();
   QuantumNumbers::QuantumNumberList ql = transform_targets(x.TransformsAs(), y.TransformsAs());
   CHECK_EQUAL(ql.size(), 1)("Operator product is not irreducible - must specify a quantum number")(ql);
   return prod(x, y, ql.front());
}

TriangularMPO prod(TriangularMPO const& x, TriangularMPO const& y, QuantumNumbers::QuantumNumber const& q)
{
   if (x.size() != y.size())
   {
      int NewSize = statistics::lcm(x.size(), y.size());
      return prod(repeat(x, NewSize/x.size()), repeat(y, NewSize/y.size()), q);
   }

   PRECONDITION(is_transform_target(y.TransformsAs(), x.TransformsAs(), q))(x.TransformsAs())(y.TransformsAs())(q)
      (x.Basis())(x.Basis().front());
   PRECONDITION_EQUAL(x.size(), y.size());

   typedef Tensor::ProductBasis<BasisList, BasisList> PBasisType;

   TriangularMPO Result(x.size());

   // The basis that wraps around gets the final element projected onto component q only
   PBasisType ProjectedBasis = PBasisType::MakeTriangularProjected(x.front().Basis1(), y.front().Basis1(), q);

   for (int Here = 0; Here < x.size(); ++Here)
   {
      Tensor::ProductBasis<BasisList, BasisList> B1 = (Here == 0) ? ProjectedBasis : PBasisType(x[Here].Basis1(), y[Here].Basis1());
      Tensor::ProductBasis<BasisList, BasisList> B2 = (Here == x.size()-1) ? ProjectedBasis : PBasisType(x[Here].Basis2(), y[Here].Basis2());

      OperatorComponent Op(x[Here].LocalBasis1(), y[Here].LocalBasis2(), B1.Basis(), B2.Basis());

      for (OperatorComponent::const_iterator I1 = iterate(x[Here].data()); I1; ++I1)
      {
	 for (OperatorComponent::const_inner_iterator I2 = iterate(I1); I2; ++I2)
	 {

	    for (OperatorComponent::const_iterator J1 = iterate(y[Here].data()); J1; ++J1)
	    {
	       for (OperatorComponent::const_inner_iterator J2 = iterate(J1); J2; ++J2)
	       {

		  ProductBasis<BasisList, BasisList>::const_iterator B1End = B1.end(I2.index1(), J2.index1());
		  ProductBasis<BasisList, BasisList>::const_iterator B2End = B2.end(I2.index2(), J2.index2());

		  ProductBasis<BasisList, BasisList>::const_iterator B1Iter = B1.begin(I2.index1(), J2.index1());

		  while (B1Iter != B1End)
		  {
		     ProductBasis<BasisList, BasisList>::const_iterator B2Iter = B2.begin(I2.index2(), J2.index2());
		     while (B2Iter != B2End)
		     {
			SimpleRedOperator ToInsert(x[Here].LocalBasis1(), y[Here].LocalBasis2());

			// iterate over components of the reducible operators *I1 and *I2
			for (SimpleRedOperator::const_iterator IComponent = I2->begin(); IComponent != I2->end(); ++IComponent)
			{
			   for (SimpleRedOperator::const_iterator JComponent = J2->begin(); JComponent != J2->end(); ++JComponent)
			   {

			      QuantumNumbers::QuantumNumberList ql = transform_targets(IComponent->TransformsAs(),
										       JComponent->TransformsAs());

			      for (QuantumNumbers::QuantumNumberList::const_iterator q = ql.begin(); q != ql.end(); ++q)
			      {
				 if (is_transform_target(B2[*B2Iter], *q, B1[*B1Iter]))
				 {
				    SimpleOperator Next = 
				       tensor_coefficient(B1, B2, 
							  IComponent->TransformsAs(), JComponent->TransformsAs(), *q,
							  I2.index1(), J2.index1(), *B1Iter,
							  I2.index2(), J2.index2(), *B2Iter) *
				       prod(*IComponent, *JComponent, *q);
				    if (!is_zero(Next.data()))
				       ToInsert += Next;

				 }
			      }
			   }
			}

			if (!ToInsert.empty())
			   Op.data()(*B1Iter, *B2Iter) = ToInsert;

			++B2Iter;
		     }
		     ++B1Iter;
		  }
	       }
	    }
	 }
      }

      Result[Here] = Op;
   }
   optimize(Result);
   return Result;
}

TriangularMPO
repeat(TriangularMPO const& x, int count)
{
   TriangularMPO Result(x.size() * count);
   for (int i = 0; i < count; ++i)
      for (int j = 0; j < x.size(); ++j)
         Result[i*x.size()+j] = x[j];
   return Result;
}

TriangularMPO dot(TriangularMPO const& x, TriangularMPO const& y)
{
   if (x.empty())
      return x;
   if (y.empty())
      return y;
   // If the MPO's are reducible then sum over every combination that
   // leads to a scalar.  I'm not sure if this is physically useful?
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());
   CHECK(is_transform_target(x.TransformsAs(), y.TransformsAs(), Ident));
   return std::sqrt(double(degree(x.TransformsAs()))) * prod(x, y, Ident);
}

TriangularMPO inner(TriangularMPO const& x, TriangularMPO const& y)
{
   return dot(adjoint(x), y);
}

TriangularMPO cross(TriangularMPO const& x, TriangularMPO const& y)
{
   CHECK(cross_product_exists(x.TransformsAs(), y.TransformsAs()))
      ("Cross product does not exist for these operators")
      (x.TransformsAs())(y.TransformsAs());

   return cross_product_factor(x.TransformsAs(), y.TransformsAs())
      * prod(x, y, cross_product_transforms_as(x.TransformsAs(), y.TransformsAs()));
}

TriangularMPO outer(TriangularMPO const& x, TriangularMPO const& y)
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

   return outer_coefficient(dx, dy, dq) *  prod(x,y,q);
}

TriangularMPO
pow(TriangularMPO const& x, int n)
{
   if (n == 0)
   {
      PANIC("TriangularOperator to power zero is not yet implemented!");
      return x;
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

TriangularMPO
operator-(TriangularMPO const& x)
{
   return x * -1.0;
}

TriangularMPO coarse_grain(TriangularMPO const& x, int N)
{
   int MinSize = statistics::lcm(N, x.size());
   return TriangularMPO(coarse_grain(repeat(x, MinSize/x.size()).data(), N).data());
}
			 
// initial operators

StateComponent Initial_E(TriangularMPO const& m, VectorBasis const& Vac)
{
   StateComponent Result(m.data().Basis1(), Vac, Vac);
   Result[0] = MatrixOperator::make_identity(Vac);
   return Result;
}

StateComponent Initial_F(TriangularMPO const& m, VectorBasis const& Vac)
{
   StateComponent Result(m.data().Basis2(), Vac, Vac);
   Result[m.data().Basis2().size()-1] = MatrixOperator::make_identity(Vac);
   return Result;
}

StateComponent Initial_E(TriangularMPO const& m)
{
   VectorBasis Vac(make_vacuum_basis(m.GetSymmetryList()));
   return Initial_E(m, Vac);
}

StateComponent Initial_F(TriangularMPO const& m)
{
   VectorBasis Vac(make_vacuum_basis(m.GetSymmetryList()));
   return Initial_F(m, Vac);
}

StateComponent Initial_E(OperatorComponent const& m, VectorBasis const& Vac)
{
   StateComponent Result(m.Basis1(), Vac, Vac);
   Result[0] = MatrixOperator::make_identity(Vac);
   return Result;
}

StateComponent Initial_F(OperatorComponent const& m, VectorBasis const& Vac)
{
   StateComponent Result(m.Basis2(), Vac, Vac);
   Result[m.Basis2().size()-1] = MatrixOperator::make_identity(Vac);
   return Result;
}



bool remove_redundant_by_row(OperatorComponent& Op)
{
   int Size = Op.Basis1().size();
   int i = 1;
   OperatorComponent const& ConstOp = Op;  // constant reference to the same object
   std::set<int> ToDelete;
   // This map says that row ToCollapseOnto[j] is equal to ToCollapseOnto[j].second * row ToCollapseOnto[j].first
   std::map<int, std::pair<int, std::complex<double> > > ToCollapseOnto;
   while (i < Size-1)
   {
      if (ToDelete.count(i))
      {
         ++i;
         continue;
      }

      // Firstly, check to see if the row is zero
      double iNormSq = 0;
      for (int j = 0; j < i; ++j)
         iNormSq += norm_frob_sq(ConstOp(i,j));

      if (iNormSq < std::numeric_limits<double>::epsilon() * 1000)
      {
         //         TRACE("Null row")(i);
         ToDelete.insert(i);
      }
      else
      {
         // not zero.  Check the other rows for something parallel, with the same diagonal
         SimpleRedOperator Diag = ConstOp(i,i);
         int j = i+1;
         while (j < Size-1)
         { 
            if (ToDelete.count(j) || norm_frob_sq(ConstOp(j,j) - Diag) >= std::numeric_limits<double>::epsilon() * 1000)
            {
               ++j;
               continue;
            }
            else
            {
               if ((norm_frob_sq(ConstOp(j,j) - Diag)) > std::numeric_limits<double>::epsilon())
                  std::cout << __LINE__ << " " << (norm_frob_sq(ConstOp(j,j) - Diag)) << '\n';
            }

            //            TRACE(i)(j);
            // see if this row is identical, up to constant factor, up to the diagonal component
            double jNormSq = 0;
            std::complex<double> ijInner;
            for (int k = 0; k < j; ++k)
            {
               jNormSq += norm_frob_sq(ConstOp(j,k));
               ijInner += inner_prod(ConstOp(i,k), ConstOp(j,k));
            }

            // if row j is zero, skip over it and continue
            if (jNormSq < std::numeric_limits<double>::epsilon()*1000)
            {
               ++j;
               continue;
            }
            //            TRACE(ijInner)(iNormSq)(jNormSq);
            // are rows i,j proportional to each other?
            if (LinearAlgebra::norm_frob(LinearAlgebra::norm_frob_sq(ijInner) - iNormSq*jNormSq) 
                < std::numeric_limits<double>::epsilon() * iNormSq*jNormSq*100)
            {
               //DEBUG_TRACE((LinearAlgebra::norm_frob_sq(ijInner) - iNormSq*jNormSq)/(iNormSq*jNormSq))(iNormSq)(jNormSq);
               // they are proportional: get the constant of proportionality
               std::complex<double> Factor = ijInner / iNormSq;
               ToCollapseOnto[j] = std::make_pair(i, Factor);
               ToDelete.insert(j);
            }
            else if (LinearAlgebra::norm_frob(LinearAlgebra::norm_frob_sq(ijInner) - iNormSq*jNormSq) 
                < std::numeric_limits<double>::epsilon() * iNormSq*jNormSq*10000000)
            {
               TRACE("smallish norm")(LinearAlgebra::norm_frob_sq(ijInner))(iNormSq)(jNormSq);
            }
            ++j;
         }
      }
      ++i;
   }

   if (ToDelete.empty())
      return false;

   // Now we assemble a new operator
   std::vector<int> NewBasisMapping(Size, -1);
   BasisList NewBasis(Op.GetSymmetryList());
   for (int i = 0; i < Size; ++i)
   {
      if (ToDelete.count(i) == 0)
      {
         NewBasisMapping[i] = NewBasis.size();
         NewBasis.push_back(Op.Basis1()[i]);
      }
   }

   OperatorComponent Result(Op.LocalBasis1(), Op.LocalBasis2(), NewBasis, NewBasis);
   for (int i = 0; i < Size; ++i)
   {
      if (NewBasisMapping[i] < 0)
         continue;

      for (int j = 0; j < i; ++j)
      {
         if (!iterate_at(ConstOp.data(), i, j))
            continue;

         if (ToCollapseOnto.find(j) != ToCollapseOnto.end())
         {
            Result(NewBasisMapping[i], NewBasisMapping[ToCollapseOnto[j].first]) += ToCollapseOnto[j].second * ConstOp(i,j);
         }
         else if (NewBasisMapping[j] >= 0)
         {
            Result(NewBasisMapping[i], NewBasisMapping[j]) += ConstOp(i,j);
         }
      }
      // the diagonal part
      if (iterate_at(ConstOp.data(), i, i))
         Result(NewBasisMapping[i], NewBasisMapping[i]) = ConstOp(i,i);
   }

   Op = Result;
   return true;
}

bool remove_redundant_by_column(OperatorComponent& Op)
{
   int Size = Op.Basis1().size();
   int i = Size-2;
   OperatorComponent const& ConstOp = Op;  // constant reference to the same object
   std::set<int> ToDelete;
   // This map says that row ToCollapseOnto[j] is equal to ToCollapseOnto[j].second * row ToCollapseOnto[j].first
   std::map<int, std::pair<int, std::complex<double> > > ToCollapseOnto;
   while (i > 0)
   {
      if (ToDelete.count(i))
      {
         --i;
         continue;
      }

      // Firstly, check to see if the column is zero
      double iNormSq = 0;
      for (int k = i+1; k < Size; ++k)
         iNormSq += norm_frob_sq(ConstOp(k,i));

      if (iNormSq < std::numeric_limits<double>::epsilon() * 1000)
      {
         ToDelete.insert(i);
      }
      else
      {
         if (iNormSq < 0.001 || i == 132)
         {
            TRACE(iNormSq)(i);
         }

         // not zero.  Check the other columns for something parallel, with the same diagonal
         SimpleRedOperator Diag = ConstOp(i,i);
         int j = i-1;
         while (j > 0)
         { 
            if (ToDelete.count(j) || norm_frob_sq(ConstOp(j,j) - Diag) >= std::numeric_limits<double>::epsilon() * 1000)
            {
               --j;
               continue;
            }

            //            TRACE(i)(j);
            // see if this row is identical, up to constant factor, up to the diagonal component
            double jNormSq = 0;
            std::complex<double> ijInner;
            for (int k = j+1; k < Size; ++k)
	    {
               jNormSq += norm_frob_sq(ConstOp(k,j));
	    }
            for (int k = j+1; k < Size; ++k)
            {
               //               TRACE(i)(j)(k)(ConstOp(k,i))(ConstOp(k,j));
               ijInner += inner_prod(ConstOp(k,i), ConstOp(k,j));
            }
            if (jNormSq < std::numeric_limits<double>::epsilon()*1000)
            {
               --j;
               continue;
            }
            //            TRACE(ijInner)(iNormSq)(jNormSq);
            // are rows i,j proportional to each other?
            if (LinearAlgebra::norm_frob(LinearAlgebra::norm_frob_sq(ijInner) - iNormSq*jNormSq) 
                < std::numeric_limits<double>::epsilon()*iNormSq*jNormSq*100)
            {
               // yes, get the constant of proportionality
               std::complex<double> Factor = ijInner / iNormSq;
               ToCollapseOnto[j] = std::make_pair(i, Factor);
               ToDelete.insert(j);
               //               TRACE("Collapsible")(i)(j)(Factor);
            }
            else if (LinearAlgebra::norm_frob(LinearAlgebra::norm_frob_sq(ijInner) 
					      - iNormSq*jNormSq) 
                < std::numeric_limits<double>::epsilon() * iNormSq*jNormSq*10000000)
            {
               TRACE("smallish norm")(LinearAlgebra::norm_frob_sq(ijInner))(iNormSq)(jNormSq);
            }

            --j;
         }
      }
      --i;
   }

   if (ToDelete.empty())
      return false;

   // Now we assemble a new operator
   std::vector<int> NewBasisMapping(Size, -1);
   BasisList NewBasis(Op.GetSymmetryList());
   for (int i = 0; i < Size; ++i)
   {
      if (ToDelete.count(i) == 0)
      {
         NewBasisMapping[i] = NewBasis.size();
         NewBasis.push_back(Op.Basis1()[i]);
      }
   }

   OperatorComponent Result(Op.LocalBasis1(), Op.LocalBasis2(), NewBasis, NewBasis);
   for (int i = 0; i < Size; ++i)
   {
      if (NewBasisMapping[i] < 0)
         continue;

      for (int j = i+1; j < Size; ++j)
      {
         if (!iterate_at(ConstOp.data(), j, i))
            continue;

         if (ToCollapseOnto.find(j) != ToCollapseOnto.end())
         {
	    DEBUG_CHECK(NewBasisMapping[ToCollapseOnto[j].first] > NewBasisMapping[i]);
            Result(NewBasisMapping[ToCollapseOnto[j].first], NewBasisMapping[i]) 
	       += ToCollapseOnto[j].second * ConstOp(j,i);
         }
         else if (NewBasisMapping[j] >= 0)
         {
	    DEBUG_CHECK(NewBasisMapping[j] > NewBasisMapping[i]);
            Result(NewBasisMapping[j], NewBasisMapping[i]) += ConstOp(j,i);
         }
      }
      // the diagonal part
      if (iterate_at(ConstOp.data(), i, i))
         Result(NewBasisMapping[i], NewBasisMapping[i]) = ConstOp(i,i);
   }

   Op = Result;
   return true;
}

void remove_redundant(OperatorComponent& Op)
{
   bool Done = false;
   while (!Done)
   {
      Done = !remove_redundant_by_row(Op) && !remove_redundant_by_column(Op);
   }
}

// smod: sensible mod operation, that is guaranteed to return a +ve result
// if k is positive.
inline
int smod(int n, int k)
{
   int r = n%k;
   if (r < 0 && k > 0)
   {
      r += k;
   }
   return r;
}

std::vector<SimpleOperator>
MakeIdentityUnitCell(std::vector<BasisList> const& Sites)
{
   std::vector<SimpleOperator> Result;
   for (unsigned i = 0; i < Sites.size(); ++i)
   {
      Result.push_back(SimpleOperator::make_identity(Sites[i]));
   }
   return Result;
}

