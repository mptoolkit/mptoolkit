// -*- C++ -*- $Id$

#include "triangular_operator.h"

MPOperator extract_column(TriangularOperator const& Op, int Col)
{
   MPOperator MPOp(Op.begin(), Op.end());

   std::set<int> Cols;
   Cols.insert(Col);
   SimpleOperator ColP = make_projector_onto(MPOp.Basis2(), Cols);
   MPOp.back() = prod(MPOp.back(), herm(ColP));

   cull_unused_elements(MPOp);
   return MPOp;
}

MPOperator extract_lower_column(TriangularOperator const& Op, int Col)
{
   MPOperator MPOp(Op.begin(), Op.end());

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

void mask_lower_column(TriangularOperator const& Op, int Col, std::vector<std::vector<int> >& Mask)
{
   initialize_mask(Op.data(), Mask);
   std::fill(Mask.back().begin(), Mask.back().end(), false);
   Mask.back()[Col] = true;
   mask_unused_elements(Op.data(), Mask);
}

MPOperator 
TriangularOperator::operator()(int Row, int Col) const
{
   MPOperator MPOp(Data_);

   std::set<int> Rows;
   Rows.insert(Row);
   MPOp.front() = project_rows(MPOp.front(), Rows);

   std::set<int> Cols;
   Cols.insert(Col);
   MPOp.back() = project_columns(MPOp.back(), Cols);

   cull_unused_elements(MPOp);
   return MPOp;
}

std::ostream&
operator<<(std::ostream& out, TriangularOperator const& op)
{
   return out << op.data();
}

TriangularOperator TriangularOneSite(SimpleOperator const& x)
{
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());
   BasisList b(x.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(x.TransformsAs());

   SimpleOperator I = SimpleOperator::make_identity(x.Basis1());

   OperatorComponent data(x.Basis1(), b, b);
   data.set_operator(0,0,I);
   data.set_operator(1,1,I);
   data.set_operator(1,0,x);

   return TriangularOperator(data);
}

TriangularOperator TriangularOneSite(SimpleOperator const& x, double Momentum)
{
   std::complex<double> r = std::exp(std::complex<double>(0.0, 1.0) * Momentum);

   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());
   BasisList b(x.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(x.TransformsAs());

   SimpleOperator I = SimpleOperator::make_identity(x.Basis1());

   OperatorComponent data(x.Basis1(), b, b);
   data.set_operator(0,0,r*I);
   data.set_operator(1,1,I);
   data.set_operator(1,0,x);

   return TriangularOperator(data);
}

TriangularOperator TriangularTwoSite(SimpleOperator const& x, SimpleOperator const& y, QuantumNumber const& Trans)
{
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());
   BasisList b(x.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(y.TransformsAs());
   b.push_back(Trans);

   SimpleOperator I = SimpleOperator::make_identity(x.Basis1());

   OperatorComponent data(x.Basis1(), b, b);
   data.set_operator(0,0,I);
   data.set_operator(2,2,I);
   data.set_operator(1,0,y);
   data.set_operator(2,1,x);

   return TriangularOperator(data);
}

TriangularOperator TriangularTwoSite(SimpleOperator const& x, SimpleOperator const& y)
{
   QuantumNumber q;
   if (num_transform_targets(x.TransformsAs(), y.TransformsAs()) == 1)
      q = transform_targets(x.TransformsAs(), y.TransformsAs())[0];
   else
      q = QuantumNumber(x.GetSymmetryList());
   return TriangularTwoSite(x, y, q);
}

TriangularOperator TriangularThreeSite(SimpleOperator const& x, SimpleOperator const& y, SimpleOperator const& z, 
				   QuantumNumber const& yz_trans)
{
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());

   BasisList b(x.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(z.TransformsAs());
   b.push_back(yz_trans);
   b.push_back(Ident);

   SimpleOperator I = SimpleOperator::make_identity(x.Basis1());

   OperatorComponent data(x.Basis1(), b, b);
   data.set_operator(0,0,I);
   data.set_operator(3,3,I);
   data.set_operator(1,0,z);
   data.set_operator(2,1,y);
   data.set_operator(3,2,x);

   return TriangularOperator(data);
}

TriangularOperator TriangularStretchedTwoSite(SimpleOperator const& x, int NumNeighbors,
					  SimpleOperator const& y)
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

   SimpleOperator I = SimpleOperator::make_identity(x.Basis1());

   OperatorComponent data(x.Basis1(), b, b);
   data.set_operator(0,0,I);
   data.set_operator(1,0,y);
   for (int i = 0; i < NumNeighbors; ++i)
   {
      data.set_operator(i+2, i+1, I);
   }
   data.set_operator(2+NumNeighbors,1+NumNeighbors,x);
   data.set_operator(2+NumNeighbors,2+NumNeighbors,I);

   return TriangularOperator(data);
}


TriangularOperator TriangularThreeSite(SimpleOperator const& x, SimpleOperator const& y, SimpleOperator const& z)
{
   QuantumNumbers::QuantumNumberList ql = transform_targets(y.TransformsAs(), z.TransformsAs());
   CHECK_EQUAL(ql.size(), 1)("ambiguous coupling of y and z")(y.TransformsAs())(z.TransformsAs());
   return TriangularThreeSite(x,y,z, ql[0]);
}

TriangularOperator TriangularFourSite(SimpleOperator const& w, SimpleOperator const& x, 
                                  SimpleOperator const& y, SimpleOperator const& z, 
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

   SimpleOperator I = SimpleOperator::make_identity(x.Basis1());

   OperatorComponent data(x.Basis1(), b, b);
   data.set_operator(0,0,I);
   data.set_operator(4,4,I);
   data.set_operator(1,0,z);
   data.set_operator(2,1,y);
   data.set_operator(3,2,x);
   data.set_operator(4,3,w);

   return TriangularOperator(data);
}

TriangularOperator TriangularFourSite(SimpleOperator const& w, SimpleOperator const& x, 
                                  SimpleOperator const& y, SimpleOperator const& z)
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


TriangularOperator ZigZagChain(SimpleOperator const& S, SimpleOperator const& T, double J1, double J2)
{
   QuantumNumbers::QuantumNumber Ident(S.GetSymmetryList());

   BasisList b(S.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(S.TransformsAs());
   b.push_back(S.TransformsAs());
   b.push_back(Ident);

   SimpleOperator I = SimpleOperator::make_identity(S.Basis1());

   OperatorComponent data(S.Basis1(), b, b);
   data.set_operator(0,0,I);
   data.set_operator(3,3,I);
   data.set_operator(1,0,S);
   data.set_operator(2,1,I);
   data.set_operator(3,2,J2*T);
   data.set_operator(3,1,J1*T);

   return TriangularOperator(data);
}

// temporary hacks for specific hamiltonians

TriangularOperator TriangularTwoSitePBC(SimpleOperator const& x, SimpleOperator const& y, QuantumNumber const& Trans)
{
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());
   BasisList b(x.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(y.TransformsAs());
   b.push_back(y.TransformsAs());
   b.push_back(Trans);

   SimpleOperator I = SimpleOperator::make_identity(x.Basis1());

   OperatorComponent data(x.Basis1(), b, b);
   data.set_operator(0,0,I);
   data.set_operator(2,2,I);
   data.set_operator(3,3,I);
   data.set_operator(1,0,y);
   data.set_operator(3,1,x);

   return TriangularOperator(data);
}

TriangularOperator TriangularTwoSitePBC(SimpleOperator const& x, SimpleOperator const& y)
{
   CHECK_EQUAL(num_transform_targets(x.TransformsAs(), y.TransformsAs()), 1);
   return TriangularTwoSitePBC(x, y, transform_targets(x.TransformsAs(), y.TransformsAs())[0]);
}

TriangularOperator TriangularTwoSitePBC_Boundary(SimpleOperator const& x, SimpleOperator const& y, QuantumNumber const& Trans)
{
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());
   BasisList b(x.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(y.TransformsAs());
   b.push_back(y.TransformsAs());
   b.push_back(Trans);

   SimpleOperator I = SimpleOperator::make_identity(x.Basis1());

   OperatorComponent data(x.Basis1(), b, b);
   data.set_operator(0,0,I);
   data.set_operator(3,3,I);
   data.set_operator(1,0,y);
   data.set_operator(2,0,y);
   data.set_operator(3,1,x);
   data.set_operator(3,2,x);

   return TriangularOperator(data);
}

TriangularOperator TriangularTwoSitePBC_Boundary(SimpleOperator const& x, SimpleOperator const& y)
{
   CHECK_EQUAL(num_transform_targets(x.TransformsAs(), y.TransformsAs()), 1);
   return TriangularTwoSitePBC_Boundary(x, y, transform_targets(x.TransformsAs(), y.TransformsAs())[0]);
}

TriangularOperator repeat(TriangularOperator const& x, int Count)
{
   std::vector<OperatorComponent> Op;
   for (int i = 0; i < Count; ++i)
   {
      for (unsigned j = 0; j < x.size(); ++j)
      {
         Op.push_back(x[j]);
      }
   }
   return TriangularOperator(Op);
}

// arithmetic

TriangularOperator& operator*=(TriangularOperator& Op, double x)
{
   for (unsigned i = 0; i < Op.size(); ++i)
   {
      for (unsigned j = 1; j < Op[i].Basis1().size(); ++j)
      {
         if (iterate_at(Op[i].data(), j, 0))
            set_element(Op[i].data(), j, 0, get_element(Op[i].data(),j,0) * x);
      }
   }
   return Op;
}

TriangularOperator& operator*=(TriangularOperator& Op, std::complex<double> x)
{
   for (unsigned i = 0; i < Op.size(); ++i)
   {
      for (unsigned j = 1; j < Op[i].Basis1().size(); ++j)
      {
         if (iterate_at(Op[i].data(), j, 0))
            set_element(Op[i].data(), j, 0, get_element(Op[i].data(),j,0) * x);
      }
   }
   return Op;
}

TriangularOperator operator*(TriangularOperator const& Op, double x)
{
   TriangularOperator Result(Op);
   Result *= x;
   return Result;
}

TriangularOperator operator*(double x, TriangularOperator const& Op)
{
   TriangularOperator Result(Op);
   Result *= x;
   return Result;
}

TriangularOperator operator*(TriangularOperator const& Op, std::complex<double> x)
{
   TriangularOperator Result(Op);
   Result *= x;
   return Result;
}

TriangularOperator operator*(std::complex<double> x, TriangularOperator const& Op)
{
   TriangularOperator Result(Op);
   Result *= x;
   return Result;
}

TriangularOperator operator+(TriangularOperator const& x, TriangularOperator const& y)
{
   if (x.size() == 0)
      return y;
   if (y.size() == 0)
      return x;

   PRECONDITION_EQUAL(x.size(), y.size());
   PRECONDITION_EQUAL(x.GetSymmetryList(), y.GetSymmetryList());
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());

   TriangularOperator Result(x.size());

   for (unsigned Here = 0; Here < x.size(); ++Here)
   {

      int xSz1 = x[Here].Basis1().size();
      int ySz1 = y[Here].Basis1().size();

      int xSz2 = x[Here].Basis2().size();
      int ySz2 = y[Here].Basis2().size();

      // This is a somewhat simplistic approach.  We just combine the first and last
      // rows/columns, and join the remainder.

      BasisList NewBasis1(x.GetSymmetryList());
      for (int i = 0; i < xSz1-1; ++i)
         NewBasis1.push_back(x[Here].Basis1()[i]);
      for (int i = 1; i < ySz1; ++i)
         NewBasis1.push_back(y[Here].Basis1()[i]);

      SimpleOperator xProjector1(NewBasis1, x[Here].Basis1());
      for (int i = 0; i < xSz1-1; ++i)         // include the first row in x only
         xProjector1(i,i) = 1.0;
      xProjector1(xSz1+ySz1-3, xSz1-1) = 1.0;   // last row is special

      SimpleOperator yProjector1(NewBasis1, y[Here].Basis1());
      yProjector1(0, 0) = 1.0;  // first row is special
      for (int i = 1; i < ySz1; ++i)
         yProjector1(i+xSz1-2,i) = 1.0;



      BasisList NewBasis2(x.GetSymmetryList());
      for (int i = 0; i < xSz2-1; ++i)
         NewBasis2.push_back(x[Here].Basis2()[i]);
      for (int i = 1; i < ySz2; ++i)
         NewBasis2.push_back(y[Here].Basis2()[i]);

      SimpleOperator xProjector2(NewBasis2, x[Here].Basis2());
      for (int i = 0; i < xSz2-1; ++i)         // include the first row in x only
         xProjector2(i,i) = 1.0;
      xProjector2(xSz2+ySz2-3, xSz2-1) = 1.0;   // last row is special

      SimpleOperator yProjector2(NewBasis2, y[Here].Basis2());
      yProjector2(0, 0) = 1.0;  // first row is special
      for (int i = 1; i < ySz2; ++i)
         yProjector2(i+xSz2-2,i) = 1.0;

      OperatorComponent Next = triple_prod(xProjector1, x[Here], herm(xProjector2)) 
         + triple_prod(yProjector1, y[Here], herm(yProjector2));

      Next(0,0) = x[Here](0,0);
      Next(Next.size1()-1, Next.size2()-1) = x[Here](x[Here].size1()-1, x[Here].size2()-1);

      Result[Here] = Next;
   }

   return Result;
}

TriangularOperator& operator+=(TriangularOperator& x, TriangularOperator const& y)
{
   x = x+y;
   return x;
}

TriangularOperator operator-(TriangularOperator const& x, TriangularOperator const& y)
{
   return x + (-1.0)*y;
}

TriangularOperator operator*(TriangularOperator const& x, TriangularOperator const& y)
{
   QuantumNumbers::QuantumNumberList ql = transform_targets(x.TransformsAs(), y.TransformsAs());
   CHECK_EQUAL(ql.size(), 1)("Operator product is not irreducible - must specify a quantum number")(ql);
   return prod(x, y, ql.front());
}

TriangularOperator prod(TriangularOperator const& x, TriangularOperator const& y, QuantumNumbers::QuantumNumber const& q)
{
   PRECONDITION(is_transform_target(y.TransformsAs(), x.TransformsAs(), q))(x.TransformsAs())(y.TransformsAs())(q);
   PRECONDITION_EQUAL(x.size(), y.size());

   typedef Tensor::ProductBasis<BasisList, BasisList> PBasisType;

   TriangularOperator Result(x.size());

   // The basis that wraps around gets the final element projected onto component q only
   PBasisType ProjectedBasis = PBasisType::MakeTriangularProjected(x.front().Basis1(), y.front().Basis1(), q);

   for (unsigned Here = 0; Here < x.size(); ++Here)
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
   return Result;
}

TriangularOperator extend(TriangularOperator const& x, int count)
{
   TriangularOperator Result(x.size() * count);
   for (int i = 0; i < count; ++i)
      for (unsigned j = 0; j < x.size(); ++j)
         Result[i*x.size()+j] = x[j];
   return Result;
}
			 
// initial operators

StateComponent Initial_E(TriangularOperator const& m, VectorBasis const& Vac)
{
   StateComponent Result(m.data().Basis1(), Vac, Vac);
   Result[m.data().Basis1().size()-1] = MatrixOperator::make_identity(Vac);
   return Result;
}

StateComponent Initial_F(TriangularOperator const& m, VectorBasis const& Vac)
{
   StateComponent Result(m.data().Basis2(), Vac, Vac);
   Result[0] = MatrixOperator::make_identity(Vac);
   return Result;
}

StateComponent Initial_E(TriangularOperator const& m)
{
   VectorBasis Vac(make_vacuum_basis(m.GetSymmetryList()));
   return Initial_E(m, Vac);
}

StateComponent Initial_F(TriangularOperator const& m)
{
   VectorBasis Vac(make_vacuum_basis(m.GetSymmetryList()));
   return Initial_F(m, Vac);
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

TriangularOperator OnePointOperator(std::vector<BasisList> const& Sites, 
                                    int n, SimpleOperator const& x)
{
   int const Size = Sites.size();

   // construct a list of added quantum numbers that are needed for the bond bases
   std::vector<BasisList> BondBasis(Size, BasisList(Sites[0].GetSymmetryList()));

   // Add the identity component to the bond basis
   for (int i = 0; i < Size; ++i)
   {
      BondBasis[i].push_back(QuantumNumber(Sites[0].GetSymmetryList()));
   }

   // the actual operator
   for (int i = 0; i < Size; ++i)
   {
      BondBasis[i].push_back(QuantumNumber(Sites[0].GetSymmetryList()));
   }

   // Assemble the operators
   MPOperator Result(Size);
   for (int i = 0; i < Size; ++i)
   {
      Result[i] = OperatorComponent(Sites[i], Sites[i], BondBasis[i], BondBasis[smod(i+1,Size)]);
      Result[i](0,0) = SimpleOperator::make_identity(Sites[i]);
      Result[i](BondBasis[i].size()-1, BondBasis[smod(i+1,Size)].size()-1) 
         = SimpleOperator::make_identity(Sites[i]);
   }

   // now the operators.  Keep track of which component we insert them into
   std::vector<int> Loc(Size, 0);
   Result[smod(n,Size)](1,0) = x;
   return TriangularOperator(Result.data());
}

TriangularOperator TwoPointOperator(std::vector<BasisList> const& Sites, 
                                    int n1, SimpleOperator const& x1,
                                    int n2, SimpleOperator const& x2)
{
   // Normal order the sites.
   if (n1 > n2)
      return TwoPointOperator(Sites, n2, x2, n1, x1);

   int const Size = Sites.size();

   // construct a list of added quantum numbers that are needed for the bond bases
   std::vector<BasisList> BondBasis(Size, BasisList(Sites[0].GetSymmetryList()));

   // Add the identity component to the bond basis
   for (int i = 0; i < Size; ++i)
   {
      BondBasis[i].push_back(QuantumNumber(Sites[0].GetSymmetryList()));
   }

   // the actual operator
   for (int i = n2; i > n1; --i)
   {
      BondBasis[smod(i,Size)].push_back(x2.TransformsAs());
   }
   
   // Finally, the Hamiltonian component
   for (int i = 0; i < Size; ++i)
   {
      BondBasis[i].push_back(QuantumNumber(Sites[0].GetSymmetryList()));
   }

   // Assemble the operators
   MPOperator Result(Size);
   for (int i = 0; i < Size; ++i)
   {
      Result[i] = OperatorComponent(Sites[i], Sites[i], BondBasis[i], BondBasis[smod(i+1,Size)]);
      Result[i](0,0) = SimpleOperator::make_identity(Sites[i]);
      Result[i](BondBasis[i].size()-1, BondBasis[smod(i+1,Size)].size()-1) 
         = SimpleOperator::make_identity(Sites[i]);
   }

   // now the operators.  Keep track of which component we insert them into
   std::vector<int> Loc(Size, 0);
   Result[smod(n2,Size)](1,0) = x2;
   ++Loc[smod(n2-1,Size)];
   for (int i = n2-1; i > n1; --i)
   {
      ++Loc[smod(i-1,Size)];
      Result[smod(i,Size)](Loc[smod(i-1,Size)], Loc[smod(i,Size)]) 
         = SimpleOperator::make_identity(Sites[smod(i,Size)]);
   }
   Result[smod(n1,Size)](Loc[smod(n1-1,Size)]+1,Loc[smod(n1,Size)]) = x1;
   return TriangularOperator(Result.data());
}

TriangularOperator TwoPointStringOperator(std::vector<BasisList> const& Sites, 
					  int n1, SimpleOperator const& x1,
					  SimpleOperator const& String,
					  int n2, SimpleOperator const& x2)
{
   CHECK(n1 < n2)("TwoPointStringOperator: error: sites must be normal ordered.")(n1)(n2);

   int const Size = Sites.size();

   // construct a list of added quantum numbers that are needed for the bond bases
   std::vector<BasisList> BondBasis(Size, BasisList(Sites[0].GetSymmetryList()));

   // Add the identity component to the bond basis
   for (int i = 0; i < Size; ++i)
   {
      BondBasis[i].push_back(QuantumNumber(Sites[0].GetSymmetryList()));
   }

   // the actual operator
   for (int i = n2; i > n1; --i)
   {
      BondBasis[smod(i,Size)].push_back(x2.TransformsAs());
   }
   
   // Finally, the Hamiltonian component
   for (int i = 0; i < Size; ++i)
   {
      BondBasis[i].push_back(QuantumNumber(Sites[0].GetSymmetryList()));
   }

   // Assemble the operators
   MPOperator Result(Size);
   for (int i = 0; i < Size; ++i)
   {
      Result[i] = OperatorComponent(Sites[i], Sites[i], BondBasis[i], BondBasis[smod(i+1,Size)]);
      Result[i](0,0) = SimpleOperator::make_identity(Sites[i]);
      Result[i](BondBasis[i].size()-1, BondBasis[smod(i+1,Size)].size()-1) 
         = SimpleOperator::make_identity(Sites[i]);
   }

   // now the operators.  Keep track of which component we insert them into
   std::vector<int> Loc(Size, 0);
   Result[smod(n2,Size)](1,0) = x2;
   ++Loc[smod(n2-1,Size)];
   for (int i = n2-1; i > n1; --i)
   {
      ++Loc[smod(i-1,Size)];
      Result[smod(i,Size)](Loc[smod(i-1,Size)], Loc[smod(i,Size)]) = String;
   }
   Result[smod(n1,Size)](Loc[smod(n1-1,Size)]+1,Loc[smod(n1,Size)]) = x1;
   return TriangularOperator(Result.data());
}
