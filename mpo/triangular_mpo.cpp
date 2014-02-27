// -*- C++ -*- $Id$

#include "triangular_mpo.h"

void
TriangularMPO::check_structure() const
{
   for (unsigned i = 0; i < this->size(); ++i)
   {
      Data_[i].check_structure();
   }
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
   TRACE(MPOp.back());

   cull_unused_elements(MPOp);
   return FiniteMPO(MPOp);
}

std::ostream&
operator<<(std::ostream& out, TriangularMPO const& op)
{
   return out << op.data();
}

TriangularMPO TriangularOneSite(SimpleOperator const& x)
{
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());
   BasisList b(x.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(x.TransformsAs());

   SimpleOperator I = SimpleOperator::make_identity(x.Basis1());

   OperatorComponent data(x.Basis1(), b, b);

   data.set_operator(0,0,I);
   data.set_operator(1,1,I);
   data.set_operator(0,1,x);

   return TriangularMPO(data);
}

TriangularMPO TriangularOneSite(SimpleOperator const& x, double Momentum)
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
   data.set_operator(0,1,x);

   return TriangularMPO(data);
}

TriangularMPO TriangularTwoSite(SimpleOperator const& x, SimpleOperator const& y, QuantumNumber const& Trans)
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
   data.set_operator(0,1,x);
   data.set_operator(1,2,y);

   return TriangularMPO(data);
}

TriangularMPO TriangularTwoSite(SimpleOperator const& x, SimpleOperator const& y)
{
   QuantumNumber q;
   if (num_transform_targets(x.TransformsAs(), y.TransformsAs()) == 1)
      q = transform_targets(x.TransformsAs(), y.TransformsAs())[0];
   else
      q = QuantumNumber(x.GetSymmetryList());
   return TriangularTwoSite(x, y, q);
}

TriangularMPO TriangularTwoSiteExponential(SimpleOperator const& x, SimpleOperator const& y, std::complex<double> Factor, QuantumNumber const& Trans)
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
   data.set_operator(1,1,Factor*I);
   data.set_operator(0,1,x);
   data.set_operator(1,2,y);

   return TriangularMPO(data);
}

TriangularMPO TriangularTwoSiteExponential(SimpleOperator const& x, SimpleOperator const& y, std::complex<double> Factor)
{
   QuantumNumber q;
   if (num_transform_targets(x.TransformsAs(), y.TransformsAs()) == 1)
      q = transform_targets(x.TransformsAs(), y.TransformsAs())[0];
   else
      q = QuantumNumber(x.GetSymmetryList());
   return TriangularTwoSiteExponential(x, y, Factor, q);
}



TriangularMPO TriangularThreeSite(SimpleOperator const& x, SimpleOperator const& y, SimpleOperator const& z, 
				   QuantumNumber const& yz_trans)
{
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());

   BasisList b(x.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(yz_trans);
   b.push_back(z.TransformsAs());
   b.push_back(Ident);

   SimpleOperator I = SimpleOperator::make_identity(x.Basis1());

   OperatorComponent data(x.Basis1(), b, b);
   data.set_operator(0,0,I);
   data.set_operator(3,3,I);
   data.set_operator(0,1,x);
   data.set_operator(1,2,y);
   data.set_operator(2,3,z);

   return TriangularMPO(data);
}

TriangularMPO TriangularStretchedTwoSite(SimpleOperator const& x, int NumNeighbors,
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
   data.set_operator(0,1,x);
   for (int i = 0; i < NumNeighbors; ++i)
   {
      data.set_operator(i+1, i+2, I);
   }
   data.set_operator(1+NumNeighbors,2+NumNeighbors,y);
   data.set_operator(2+NumNeighbors,2+NumNeighbors,I);

   return TriangularMPO(data);
}


TriangularMPO TriangularThreeSite(SimpleOperator const& x, SimpleOperator const& y, SimpleOperator const& z)
{
   QuantumNumbers::QuantumNumberList ql = transform_targets(y.TransformsAs(), z.TransformsAs());
   CHECK_EQUAL(ql.size(), 1)("ambiguous coupling of y and z")(y.TransformsAs())(z.TransformsAs());
   return TriangularThreeSite(x,y,z, ql[0]);
}

TriangularMPO TriangularFourSite(SimpleOperator const& w, SimpleOperator const& x, QuantumNumber const& wx_trans,
                                 SimpleOperator const& y, QuantumNumber const& wxy_trans,
                                 SimpleOperator const& z)
{
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());

   BasisList b(w.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(adjoint(w.TransformsAs()));
   b.push_back(adjoint(wx_trans));
   b.push_back(adjoint(wxy_trans));
   b.push_back(Ident);

   SimpleOperator I = SimpleOperator::make_identity(x.Basis1());

   OperatorComponent data(x.Basis1(), b, b);
   data.set_operator(0,0,I);
   data.set_operator(4,4,I);
   data.set_operator(0,1,w);
   data.set_operator(1,2,x);
   data.set_operator(2,3,y);
   data.set_operator(3,4,z);

   return TriangularMPO(data);
}

TriangularMPO TriangularFourSite(SimpleOperator const& w, SimpleOperator const& x, 
                                  SimpleOperator const& y, SimpleOperator const& z)
{
   QuantumNumbers::QuantumNumberList ql = transform_targets(w.TransformsAs(), x.TransformsAs());
   CHECK_EQUAL(ql.size(), 1)("ambiguous coupling of w and x")(w.TransformsAs())(x.TransformsAs());
   QuantumNumber wx = ql[0];

   ql = transform_targets(wx, y.TransformsAs());
   CHECK_EQUAL(ql.size(), 1)("ambiguous coupling of wx and y")(y.TransformsAs())
      (y.TransformsAs())(z.TransformsAs());
   QuantumNumber wxy = ql[0];

   return TriangularFourSite(w, x, wx, y, wxy, z);
}

#if 0
TriangularMPO ZigZagChain(SimpleOperator const& S, SimpleOperator const& T, double J1, double J2)
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

   return TriangularMPO(data);
}

// temporary hacks for specific hamiltonians

TriangularMPO TriangularTwoSitePBC(SimpleOperator const& x, SimpleOperator const& y, QuantumNumber const& Trans)
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

   return TriangularMPO(data);
}

TriangularMPO TriangularTwoSitePBC(SimpleOperator const& x, SimpleOperator const& y)
{
   CHECK_EQUAL(num_transform_targets(x.TransformsAs(), y.TransformsAs()), 1);
   return TriangularTwoSitePBC(x, y, transform_targets(x.TransformsAs(), y.TransformsAs())[0]);
}

TriangularMPO TriangularTwoSitePBC_Boundary(SimpleOperator const& x, SimpleOperator const& y, QuantumNumber const& Trans)
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

   return TriangularMPO(data);
}

TriangularMPO TriangularTwoSitePBC_Boundary(SimpleOperator const& x, SimpleOperator const& y)
{
   CHECK_EQUAL(num_transform_targets(x.TransformsAs(), y.TransformsAs()), 1);
   return TriangularTwoSitePBC_Boundary(x, y, transform_targets(x.TransformsAs(), y.TransformsAs())[0]);
}
#endif

// arithmetic

TriangularMPO& operator*=(TriangularMPO& Op, double x)
{
   for (unsigned i = 0; i < Op.size(); ++i)
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
   for (unsigned i = 0; i < Op.size(); ++i)
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

   PRECONDITION_EQUAL(x.size(), y.size());
   PRECONDITION_EQUAL(x.GetSymmetryList(), y.GetSymmetryList());
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());

   TriangularMPO Result(x.size());

   for (unsigned Here = 0; Here < x.size(); ++Here)
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
   PRECONDITION(is_transform_target(y.TransformsAs(), x.TransformsAs(), q))(x.TransformsAs())(y.TransformsAs())(q);
   PRECONDITION_EQUAL(x.size(), y.size());

   typedef Tensor::ProductBasis<BasisList, BasisList> PBasisType;

   TriangularMPO Result(x.size());

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

TriangularMPO
repeat(TriangularMPO const& x, int count)
{
   TriangularMPO Result(x.size() * count);
   for (int i = 0; i < count; ++i)
      for (unsigned j = 0; j < x.size(); ++j)
         Result[i*x.size()+j] = x[j];
   return Result;
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

TriangularMPO OnePointStringOperator(std::vector<BasisList> const& Sites, 
					  std::vector<SimpleOperator> const& String,
					  int n, SimpleOperator const& x, double Momentum)
{
   CHECK(n >= 0 && n < int(Sites.size()))(n)(Sites.size());
   CHECK_EQUAL(Sites.size(), String.size());
   CHECK_EQUAL(x.Basis1(), Sites[n]);
   CHECK_EQUAL(x.Basis2(), Sites[n]);
   int const Size = Sites.size();
   std::complex<double> r = std::exp(std::complex<double>(0.0, 1.0) * (Momentum / Size));

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
      BondBasis[i].push_back(x.TransformsAs());
   }

   // Assemble the operators
   GenericMPO Result(Size);
   for (int i = 0; i < Size; ++i)
   {
      CHECK_EQUAL(String[i].Basis1(), Sites[i]);
      CHECK_EQUAL(String[i].Basis2(), Sites[i]);

      Result[i] = OperatorComponent(Sites[i], Sites[i], BondBasis[i], BondBasis[smod(i+1,Size)]);
      Result[i](0,0) = i == 0 ? (r * String[i]) : r * String[i];
      Result[i](1,1) = SimpleOperator::make_identity(Sites[i]);
   }

   Result[n](0,1) = x;
   return TriangularMPO(Result.data());
}

TriangularMPO OnePointOperator(std::vector<BasisList> const& Sites, 
                                    int n, SimpleOperator const& x, double Momentum)
{
   return OnePointStringOperator(Sites, MakeIdentityUnitCell(Sites), n, x, Momentum);
}

TriangularMPO TwoPointOperator(std::vector<BasisList> const& Sites, 
                                    int n1, SimpleOperator const& x1,
                                    int n2, SimpleOperator const& x2)
{
   DEBUG_TRACE(n1)(x1)(n2)(x2);
   // Normal order the sites.
   if (n1 > n2)
      return TwoPointOperator(Sites, n2, x2, n1, x1);
   else if (n1 == n2)
      return OnePointOperator(Sites, n1, x1*x2);

   int const Size = Sites.size();

   // construct a list of added quantum numbers that are needed for the bond bases
   std::vector<BasisList> BondBasis(Size, BasisList(Sites[0].GetSymmetryList()));

   // Add the identity component to the bond basis
   for (int i = 0; i < Size; ++i)
   {
      BondBasis[i].push_back(QuantumNumber(Sites[0].GetSymmetryList()));
   }

   // the actual operator
   for (int i = n1+1; i <= n2; ++i)
   {
      BondBasis[smod(i,Size)].push_back(x2.TransformsAs());
   }
   
   // Finally, the Hamiltonian component
   for (int i = 0; i < Size; ++i)
   {
      BondBasis[i].push_back(QuantumNumber(Sites[0].GetSymmetryList()));
   }

   // Assemble the operators
   GenericMPO Result(Size);
   for (int i = 0; i < Size; ++i)
   {
      Result[i] = OperatorComponent(Sites[i], Sites[i], BondBasis[i], BondBasis[smod(i+1,Size)]);
      Result[i](0,0) = SimpleOperator::make_identity(Sites[i]);
      Result[i](BondBasis[i].size()-1, BondBasis[smod(i+1,Size)].size()-1) 
         = SimpleOperator::make_identity(Sites[i]);
   }

   // now the operators.  Keep track of which component we insert them into
   std::vector<int> Loc(Size, 0);
   Result[smod(n1,Size)](0,1) = x1;
   ++Loc[smod(n1+1,Size)];
   for (int i = n1+1; i < n2; ++i)
   {
      Result[smod(i,Size)](Loc[smod(i,Size)], Loc[smod(i+1,Size)]+1) 
         = SimpleOperator::make_identity(Sites[smod(i,Size)]);
      ++Loc[smod(i+1,Size)];
   }
   //   TRACE(n1)(n2)(LinearAlgebra::Vector<int>(Loc));
   Result[smod(n2,Size)](Loc[smod(n2,Size)],Loc[smod(n2+1,Size)]+1) = x2;
   TriangularMPO TriResult(Result.data());
   //TRACE(TriResult);
   TriResult.debug_check_structure();
   return TriResult;
}

TriangularMPO TwoPointStringOperator(std::vector<BasisList> const& Sites, 
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
   for (int i = n1+1; i <= n2; ++i)
   {
      BondBasis[smod(i,Size)].push_back(x2.TransformsAs());
   }
   
   // Finally, the Hamiltonian component
   for (int i = 0; i < Size; ++i)
   {
      BondBasis[i].push_back(QuantumNumber(Sites[0].GetSymmetryList()));
   }

   // Assemble the operators
   GenericMPO Result(Size);
   for (int i = 0; i < Size; ++i)
   {
      Result[i] = OperatorComponent(Sites[i], Sites[i], BondBasis[i], BondBasis[smod(i+1,Size)]);
      Result[i](0,0) = SimpleOperator::make_identity(Sites[i]);
      Result[i](BondBasis[i].size()-1, BondBasis[smod(i+1,Size)].size()-1) 
         = SimpleOperator::make_identity(Sites[i]);
   }

   // now the operators.  Keep track of which component we insert them into
   std::vector<int> Loc(Size, 0);
   Result[smod(n1,Size)](0,1) = x1;
   ++Loc[smod(n1+1,Size)];
   for (int i = n1+1; i < n2; ++i)
   {
      Result[smod(i,Size)](Loc[smod(i,Size)], Loc[smod(i+1,Size)]+1) = String;
      ++Loc[smod(i+1,Size)];
   }
   Result[smod(n2,Size)](Loc[smod(n2,Size)],Loc[smod(n2+1,Size)]+1) = x2;
   return TriangularMPO(Result.data());
}
