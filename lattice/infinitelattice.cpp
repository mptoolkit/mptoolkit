// -*- C++ -*- $Id$

#include "infinitelattice.h"

InfiniteLattice::InfiniteLattice()
{
}

InfiniteLattice::InfiniteLattice(UnitCell const& uc)
   : UnitCell_(uc)
{
}

bool
InfiniteLattice::triangular_operator_exists(std::string const& s) const
{
   triangular_map_type::const_iterator I = OperatorMap_.find(s);
   return I != OperatorMap_.end();
}

TriangularMPO const&
InfiniteLattice::TriangularOperator(std::string const& Op) const
{
   triangular_map_type::const_iterator I = OperatorMap_.find(Op);
   CHECK(I != OperatorMap_.end())("Operator does not exist in the lattice!")(Op);
   return I->second;
}

TriangularMPO const&
InfiniteLattice::operator[](std::string const& Op) const
{
   triangular_map_type::const_iterator I = OperatorMap_.find(Op);
   CHECK(I != OperatorMap_.end())("Operator does not exist in the lattice!")(Op);
   return I->second;
}

TriangularMPO&
InfiniteLattice::operator[](std::string const& Op)
{
   triangular_map_type::iterator I = OperatorMap_.find(Op);
   if (I == OperatorMap_.end())
   {
      OperatorMap_[Op] = make_zero(this->GetUnitCell());
      I = OperatorMap_.find(Op);
   }
   return I->second;
}

TriangularMPO
InfiniteLattice::TriangularOperatorFunction(std::string const& Op,
					    std::vector<std::complex<double> > const& Params) const
{
   PANIC("Triangular operator functions are not yet implemented");
   return TriangularMPO();
}

void
InfiniteLattice::add_triangular(std::string const& Name, InfiniteMPO const& Op)
{
   OperatorMap_[Name] = Op;
}

PStream::opstream&
operator<<(PStream::opstream& out, InfiniteLattice const& L)
{
   out << L.UnitCell_;
   out << L.OperatorMap_;
   return out;
}

PStream::ipstream&
operator>>(PStream::ipstream& in, InfiniteLattice& L)
{
   in >> L.UnitCell_;
   in >> L.OperatorMap_;
   return in;
}

// utility to join another basis list to the end of an existing basis
void
JoinBasis(BasisList& b, BasisList const& Other)
{
   for (BasisList::const_iterator I = Other.begin(); I != Other.end(); ++I)
   {
      b.push_back(*I);
   }
}

// utility to set a subset of components of C to the same as Op.
// Since we can't do sections of a sparse matrix.
void
SetComponents(OperatorComponent& C, OperatorComponent const& Op, int xStart, int yStart)
{
   for (OperatorComponent::const_iterator I = iterate(Op); I; ++I)
   {
      for (OperatorComponent::const_inner_iterator J = iterate(I); J; ++J)
      {
	 C.data()(J.index1()+xStart, J.index2()+yStart) = *J;
      }
   }
}

TriangularMPO sum_unit(UnitCellMPO const& Op)
{
   return sum_unit(Op.GetUnitCell(), Op.MPO(), Op.Commute());
}

TriangularMPO sum_unit(UnitCell const& Cell, FiniteMPO const& Op, LatticeCommute Com)
{
   CHECK(Op.is_irreducible());
   CHECK(Op.size() % Cell.size() == 0)
      ("Operator for sum_unit() must be a multiple of the unit cell");

   // Suppose that the unit cell size is 1.  Then if Op is A \times B \otimes C, then
   // the resulting TriangularMPO is
   // ( X A 0 0 )
   // ( 0 0 B 0 )
   // ( 0 0 0 C )
   // ( 0 0 0 I )
   // where X is the commutation operator.

   // If the unit cell is larger than 1 site, we can follow a more general procedure,
   // eg 6-site operator ABCDEF on a 3-site unit cell, construct operators
   // (X 0 0 0) (X 0 0 0) (X 0 0 0)
   // (0 A 0 0) (0 B 0 0) (0 C 0 0)
   // (0 0 D 0) (0 0 E 0) (0 0 F 0)
   // (0 0 0 I) (0 0 0 I) (0 0 0 I)
   // and merge the first two rows of the first site, and the
   // last two columns of the last site, to give:
   // (X A 0 0) (X 0 0 0) (X 0 0)
   // (0 0 D 0) (0 B 0 0) (0 C 0)
   // (0 0 0 I) (0 0 E 0) (0 0 F)
   //           (0 0 0 I) (0 0 I)
   //
   // Upon coarse-graining this gives:
   // (XXX ABC 000)
   // (000 000 DEF)
   // (000 000 III)
   // which is the desired operator.  This is quite straightforward,
   // since X,A have 1 row, and F,I have 1 column.

   TRACE(Op.size())(Cell.size());

   // To construct this operator, we firstly split Op into unit cells
   std::vector<std::vector<OperatorComponent> > SplitOp;
   int n = 0;
   FiniteMPO::const_iterator I = Op.begin();
   FiniteMPO::const_iterator J = I;
   while (J != Op.end())
   {
      std::advance(J, Cell.size());
      SplitOp.push_back(std::vector<OperatorComponent>(I,J));
      I = J;
   }

   // And we need th JW string
   FiniteMPO JW = Cell.string_mpo(Com.SignOperator(), Op.qn1());
   
   // and finally the identity operator
   FiniteMPO Ident = Cell.identity_mpo(Op.qn2());

   TriangularMPO Result(Cell.size());
   for (int i = 0; i < Cell.size(); ++i)
   {
      // Construct the basis
      BasisList Basis1(Cell.GetSymmetryList());
      BasisList Basis2(Cell.GetSymmetryList());
      if (i != 0)
	 JoinBasis(Basis1, JW[i].Basis1());
      JoinBasis(Basis2, JW[i].Basis2());

      for (unsigned n = 0; n < SplitOp.size(); ++n)
      {
	 JoinBasis(Basis1, SplitOp[n][i].Basis1());
	 JoinBasis(Basis2, SplitOp[n][i].Basis2());
      }
      JoinBasis(Basis1, Ident[i].Basis1());
      if (i != int(Cell.size())-1)
	 JoinBasis(Basis2, Ident[i].Basis2());

      TRACE(Basis1)(Basis2)(SplitOp.size());

      // Construct the OperatorComponent
      OperatorComponent C(Op[i].LocalBasis1(), Op[i].LocalBasis2(), Basis1, Basis2);

      // The JW goes in the top left
      int r = 0;
      int c = 0;
      SetComponents(C, JW[i], r, c);
      if (i != 0)
	 r += JW[i].Basis1().size();
      c += JW[i].Basis2().size();

      // the finite MPO components go along the diagonal
      for (unsigned n = 0; n < SplitOp.size()-1; ++n)
      {
	 SetComponents(C, SplitOp[n][i], r, c);
	 r += SplitOp[n][i].Basis1().size();
	 c += SplitOp[n][i].Basis2().size();
      }
      SetComponents(C, SplitOp.back()[i], r, c);
      r += SplitOp.back()[i].Basis1().size();
      if (i != int(Cell.size())-1)
	 c += SplitOp.back()[i].Basis2().size();
      // The identity goes in the bottom right
      SetComponents(C, Ident[i], r, c);

      // check that we're at the end
      CHECK_EQUAL(r+Ident[i].Basis1().size(), Basis1.size());
      CHECK_EQUAL(c+Ident[i].Basis2().size(), Basis2.size());

      TRACE(C);
      
      C.debug_check_structure();

      Result[i] = C;
   }

   Result.check_structure();
   return Result;
}

TriangularMPO make_zero(UnitCell const& Cell)
{
   // The auxiliary basis is the same at every site
   // Construct the basis
   BasisList b(Cell.GetSymmetryList());
   QuantumNumbers::QuantumNumber Ident(Cell.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(Ident);

   TriangularMPO Result(Cell.size());
   for (int i = 0; i < Cell.size(); ++i)
   {
      
      OperatorComponent C(Cell[i].Basis1(), Cell[i].Basis2(), b, b);
      C(0,0) = SimpleOperator::make_identity(Cell[i].Basis1());
      C(1,1) = SimpleOperator::make_identity(Cell[i].Basis1());
      Result[i] = C;
   }

   Result.check_structure();
   return Result;
}
