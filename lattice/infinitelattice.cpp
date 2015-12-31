// -*- C++ -*- $Id$

#include "infinitelattice.h"
#include "lattice/infinite-parser.h"
#include "mp/copyright.h" // for EscapeArgument

PStream::VersionTag LatticeVersion(2);

InfiniteLattice::InfiniteLattice()
{
}

InfiniteLattice::InfiniteLattice(UnitCell const& uc)
   : UnitCell_(uc)
{
}

InfiniteLattice::InfiniteLattice(std::string const& Description, UnitCell const& uc)
   : Description_(Description), UnitCell_(uc)
{
}

void
InfiniteLattice::set_command_line(int argc, char** argv)
{
   std::ostringstream Out;
   Out << EscapeArgument(argv[0]);
   for (int i = 1; i < argc; ++i)
      Out << ' ' << EscapeArgument(argv[i]);
   CommandLine_ = Out.str();

   // timestamp
   time_t now = time(NULL);
   char s[200];
   int const max = 200;
   strftime(s, max, "%a, %d %b %Y %T %z", localtime(&now));
   Timestamp_ = s;
}

// operators

void
InfiniteLattice::set_operator_descriptions(OperatorDescriptions const& Desc)
{
   // iterate through the descriptions
   for (OperatorDescriptions::const_iterator I = Desc.begin(); I != Desc.end(); ++I)
   {
      if (this->operator_exists(I->first))
      {
	 Operators_[I->first].set_description(I->second);
      }
      else
      {
	 std::cerr << "warning: operator " << I->first 
		   << " has a description but is not defined in the lattice.\n";
      }
   }

   // Now go through the operators and check for any that don't have a description
   for (const_operator_iterator I = this->begin_operator(); I != this->end_operator(); ++I)
   {
      if (I->second.description().empty())
      {
	 std::cerr << "warning: lattice operator " << I->first << " has no description.\n";
      }
   }

   // Same for the functions
   for (OperatorDescriptions::const_function_iterator I = Desc.begin_function();
	I != Desc.end_function(); ++I)
   {
      if (this->function_exists(I->first))
      {
	 Functions_[I->first].set_description(I->second);
      }
      else
      {
	 std::cerr << "warning: function " << I->first
		   << " has a descrption but is not defined in the lattice.\n";
      }
   }

   // Check that all functions have a definition
   for (const_function_iterator I = this->begin_function(); I != this->end_function(); ++I)
   {
      if (I->second.description().empty())
      {
	 std::cerr << "warning: lattice function " << I->first << " has no description.\n";
      }
   }
}

bool
InfiniteLattice::operator_exists(std::string const& s) const
{
   return (Operators_.find(s) != Operators_.end());
}

bool
InfiniteLattice::triangular_operator_exists(std::string const& s) const
{
   OperatorListType::const_iterator I = Operators_.find(s);
   return I != Operators_.end() && I->second.is_triangular();
}

bool
InfiniteLattice::product_operator_exists(std::string const& s) const
{
   OperatorListType::const_iterator I = Operators_.find(s);
   return I != Operators_.end() && I->second.is_product();
}

InfiniteMPO const&
InfiniteLattice::operator[](std::string const& Op) const
{
   OperatorListType::const_iterator I = Operators_.find(Op);
   CHECK(I != Operators_.end())("Operator does not exist in the lattice!")(Op);
   return I->second;
}

InfiniteMPO&
InfiniteLattice::operator[](std::string const& Op)
{
   return Operators_[Op];
}

TriangularMPO const&
InfiniteLattice::as_triangular_mpo(std::string const& Op) const
{
   OperatorListType::const_iterator I = Operators_.find(Op);
   CHECK(I != Operators_.end())("Operator does not exist in the lattice!")(Op);
   return I->second.as_triangular_mpo();
}

ProductMPO const&
InfiniteLattice::as_product_mpo(std::string const& Op) const
{
   OperatorListType::const_iterator I = Operators_.find(Op);
   CHECK(I != Operators_.end())("Operator does not exist in the lattice!")(Op);
   return I->second.as_product_mpo();
}

// arguments

std::complex<double>
InfiniteLattice::arg(std::string const& a) const
{
   const_argument_iterator I = this->find_arg(a);
   if (I != this->end_arg())
      return I->second;
   return 0.0;
}

// functions

InfiniteLattice::function_type
InfiniteLattice::func(std::string const& s) const
{
   const_function_iterator I = this->find_function(s);
   CHECK(I != this->end_function())("Function not found")(s);
   return I->second;
}

// functor to parse default arguments of a UnitCell operator
struct ParseInfiniteArgument
{
   ParseInfiniteArgument(InfiniteLattice const& Lattice_) : Lattice(Lattice_) {}

   std::complex<double> operator()(Function::ArgumentList const& Args,
				   std::string const& Str) const
   {
      return ParseInfiniteNumber(Lattice, Str, Args);
   }

   InfiniteLattice const& Lattice;
};

InfiniteMPO
InfiniteLattice::eval_function(Function::OperatorFunction const& Func,
			       Function::ParameterList const& Params) const
{
   Function::ArgumentList Args = GetArguments(Func.Args, Params, 
					      ParseInfiniteArgument(*this));
   return ParseInfiniteOperator(*this, Func.Def, Args);
}

InfiniteMPO
InfiniteLattice::eval_function(std::string const& Func,
			       Function::ParameterList const& Params) const
{
   return this->eval_function(this->func(Func), Params);
}


PStream::opstream&
operator<<(PStream::opstream& out, InfiniteLattice const& L)
{
   PStream::VersionSentry Sentry(out, LatticeVersion, LatticeVersion.default_version());
   out << Sentry.version();
   out << L.Description_;
   out << L.CommandLine_;
   out << L.Timestamp_;
   out << L.UnitCell_;
   out << L.Operators_;
   out << L.Arguments_;
   out << L.Functions_;
   return out;
}

PStream::ipstream&
operator>>(PStream::ipstream& in, InfiniteLattice& L)
{
   PStream::VersionSentry Sentry(in, LatticeVersion, in.read<int>());

   if (Sentry.version() > 2)
   {
      PANIC("This program is too old to read this lattice file format,"
	    "  Maximum readable version number is 2")(Sentry.version());
   } else if (Sentry.version() < 0)
   {
      PANIC("Lattice file is too old, please reconstruct the lattice."
	    "  Expected Version >= 0")(Sentry.version());
   }

   in >> L.Description_;
   in >> L.CommandLine_;
   in >> L.Timestamp_;
   in >> L.UnitCell_;
   in >> L.Operators_;
   in >> L.Arguments_;
   in >> L.Functions_;
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
	 CHECK_EQUAL(C.Basis1()[J.index1()+xStart], Op.Basis1()[J.index1()])(J.index1())(xStart);
	 CHECK_EQUAL(C.Basis2()[J.index2()+yStart], Op.Basis2()[J.index2()])(J.index2())(yStart);
	 C.data()(J.index1()+xStart, J.index2()+yStart) = *J;
      }
   }
}

TriangularMPO sum_unit(UnitCellMPO const& Op)
{
   if (Op.is_null())
      return TriangularMPO();
   return sum_unit(*Op.GetSiteList(), Op.MPO(), Op.Commute(), Op.GetSiteList()->size());
}

TriangularMPO sum_unit(UnitCellMPO const& Op, int UnitCellSize)
{
   if (Op.is_null())
      return TriangularMPO();
   return sum_unit(*Op.GetSiteList(), ExtendToCoverUnitCell(Op, UnitCellSize).MPO(), Op.Commute(), UnitCellSize);
}

TriangularMPO sum_unit(UnitCell const& Cell, FiniteMPO const& Op, LatticeCommute Com)
{
   return sum_unit(*Cell.GetSiteList(), Op, Com, Cell.size());
}

TriangularMPO sum_unit(SiteListType const& SiteList, FiniteMPO const& Op, LatticeCommute Com, int UnitCellSize)
{
   return sum_unit(SiteList, string_mpo(SiteList, Com.SignOperator(), Op.qn1()), Op, UnitCellSize);
}

std::vector<std::vector<OperatorComponent>> 
   SplitOperator(FiniteMPO const& Op, int UnitCellSize)
{
   DEBUG_CHECK_EQUAL(Op.size() % UnitCellSize, 0);
   std::vector<std::vector<OperatorComponent> > Result;
   FiniteMPO::const_iterator I = Op.begin();
   FiniteMPO::const_iterator J = I;
   while (J != Op.end())
   {
      std::advance(J, UnitCellSize);
      Result.push_back(std::vector<OperatorComponent>(I,J));
      I = J;
   }
   return Result;
}

TriangularMPO sum_unit(SiteListType const& SiteList, FiniteMPO const& JW2, FiniteMPO const& Op2, int UnitCellSize)
{
   FiniteMPO Op = Op2;
   optimize(Op);

   FiniteMPO JW = JW2;
   optimize(JW);

   if (Op.is_null())
      return TriangularMPO();
   CHECK(Op.is_irreducible());
   CHECK(UnitCellSize % SiteList.size() == 0)
      ("UnitCellSize for sum_unit() must be a multiple of UnitCell.size()");
   CHECK(Op.size() % UnitCellSize == 0)
      ("Operator for sum_unit() must be a multiple of the unit cell")
      (Op.size())(UnitCellSize);
   // THe JW string operator must have a unit cell that is a divisor of UnitCellSize
   CHECK(UnitCellSize % JW.size() == 0)
      ("JW string for sum_unit() must divide UnitCell.size()");

   // Suppose that the unit cell size is 1.  Then if Op is A \times B \otimes C, then
   // the resulting TriangularMPO is
   // ( X A 0 0 )
   // ( 0 0 B 0 )
   // ( 0 0 0 C )
   // ( 0 0 0 I )
   // where X is the JW string operator.

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

   DEBUG_TRACE(Op.size())(UnitCellSize)(SiteList.size());

   // To construct this operator, we firstly split Op into UnitCellSize pieces
   std::vector<std::vector<OperatorComponent>> SplitOp = SplitOperator(Op, UnitCellSize);

   // and finally the identity operator
   FiniteMPO Ident = identity_mpo(SiteList, Op.qn2());

   TriangularMPO Result(UnitCellSize);
   for (int i = 0; i < UnitCellSize; ++i)
   {
      // Construct the basis
      BasisList Basis1(Op.GetSymmetryList());
      BasisList Basis2(Op.GetSymmetryList());
      if (i != 0)
	 JoinBasis(Basis1, JW[i%JW.size()].Basis1());
      JoinBasis(Basis2, JW[i%JW.size()].Basis2());

      for (unsigned n = 0; n < SplitOp.size(); ++n)
      {
	 JoinBasis(Basis1, SplitOp[n][i].Basis1());
	 JoinBasis(Basis2, SplitOp[n][i].Basis2());
      }
      JoinBasis(Basis1, Ident[i%SiteList.size()].Basis1());
      if (i != int(UnitCellSize)-1)
	 JoinBasis(Basis2, Ident[i%SiteList.size()].Basis2());

      DEBUG_TRACE(Basis1)(Basis2)(SplitOp.size());

      // Construct the OperatorComponent
      OperatorComponent C(Op[i].LocalBasis1(), Op[i].LocalBasis2(), Basis1, Basis2);

      // The JW goes in the top left
      int r = 0;
      int c = 0;
      SetComponents(C, JW[i%JW.size()], r, c);
      if (i != 0)
	 r += JW[i%JW.size()].Basis1().size();
      c += JW[i%JW.size()].Basis2().size();

      // the finite MPO components go along the diagonal
      for (unsigned n = 0; n < SplitOp.size()-1; ++n)
      {
	 SetComponents(C, SplitOp[n][i], r, c);
	 r += SplitOp[n][i].Basis1().size();
	 c += SplitOp[n][i].Basis2().size();
      }
      SetComponents(C, SplitOp.back()[i], r, c);
      r += SplitOp.back()[i].Basis1().size();
      if (i != int(UnitCellSize)-1)
	 c += SplitOp.back()[i].Basis2().size();
      // The identity goes in the bottom right
      SetComponents(C, Ident[i%SiteList.size()], r, c);

      // check that we're at the end
      CHECK_EQUAL(r+Ident[i%SiteList.size()].Basis1().size(), Basis1.size());
      CHECK_EQUAL(c+Ident[i%SiteList.size()].Basis2().size(), Basis2.size());

      DEBUG_TRACE(C);
      
      C.debug_check_structure();

      Result[i] = C;
   }

   Result.check_structure();
   optimize(Result);
   return Result;
}


TriangularMPO sum_kink(UnitCellMPO const& Kink, UnitCellMPO const& Op, int UnitCellSize)
{
   if (Op.is_null() || Kink.is_null())
      return TriangularMPO();
   return sum_kink(*Op.GetSiteList(), 
		   ExtendToCoverUnitCell(Kink, UnitCellSize).MPO(), 
		   ExtendToCoverUnitCell(Op, UnitCellSize).MPO(), 
		   Op.Commute(), 
		   UnitCellSize);
}

TriangularMPO sum_kink(UnitCellMPO const& Kink, UnitCellMPO const& Op)
{
   if (Op.is_null() || Kink.is_null())
      return TriangularMPO();
   CHECK_EQUAL(Kink.GetSiteList()->size(), Op.GetSiteList()->size())
      ("Operators for sum_kink must have the same unit cell!");
   return sum_kink(*Op.GetSiteList(), Kink.MPO(), Op.MPO(), Op.Commute(), Op.GetSiteList()->size());
}

TriangularMPO sum_kink(SiteListType const& SiteList, FiniteMPO const& Kink, FiniteMPO const& Op, LatticeCommute Com, int UnitCellSize)
{
   FiniteMPO Ident = repeat(string_mpo(SiteList, Com.SignOperator(), Op.qn1()), Kink.size() / SiteList.size());
   return sum_unit(SiteList, Kink*Ident, Op, UnitCellSize);
}

TriangularMPO sum_k(SiteListType const& SiteList, std::complex<double> const& k,
		       FiniteMPO const& Op, LatticeCommute Com, int UnitCellSize)
{
   return sum_unit(SiteList, exp(std::complex<double>(0,1)*k)
		   *string_mpo(SiteList, Com.SignOperator(), Op.qn1()), Op, UnitCellSize);
}
   
TriangularMPO sum_k(std::complex<double> const& k, UnitCellMPO const& Op, int UnitCellSize)
{
   if (Op.is_null())
      return TriangularMPO();
   return sum_unit(*Op.GetSiteList(), exp(std::complex<double>(0,1)*k)
		   *string_mpo(*Op.GetSiteList(), Op.Commute().SignOperator(), Op.qn1()), Op.MPO(), UnitCellSize);
}   

TriangularMPO sum_k(std::complex<double> const& k, UnitCellMPO const& Op)
{
   if (Op.is_null())
      return TriangularMPO();
   return sum_unit(*Op.GetSiteList(), exp(std::complex<double>(0,1)*k)
		   *string_mpo(*Op.GetSiteList(), Op.Commute().SignOperator(), Op.qn1()), Op.MPO(), 
		   Op.GetSiteList()->size());
}

TriangularMPO make_zero(SiteListType const& SiteList)
{
   if (SiteList.empty())
      return TriangularMPO();

   // The auxiliary basis is the same at every site
   // Construct the basis
   BasisList b(SiteList[0].GetSymmetryList());
   QuantumNumbers::QuantumNumber Ident(SiteList[0].GetSymmetryList());
   b.push_back(Ident);
   b.push_back(Ident);

   TriangularMPO Result(SiteList.size());
   for (unsigned i = 0; i < SiteList.size(); ++i)
   {
      
      OperatorComponent C(SiteList[i].Basis1(), SiteList[i].Basis2(), b, b);
      C(0,0) = SimpleOperator::make_identity(SiteList[i].Basis1());
      C(1,1) = SimpleOperator::make_identity(SiteList[i].Basis1());
      Result[i] = C;
   }

   Result.check_structure();
   optimize(Result);
   return Result;
}

//
// sum_string
//

TriangularMPO sum_string(SiteListType const& SiteList, FiniteMPO const& JW, FiniteMPO const& Op1,
			 FiniteMPO const& String, FiniteMPO const& Op2, int UnitCellSize,
			 QuantumNumbers::QuantumNumber q)
{
   FiniteMPO X = JW;
   FiniteMPO A = Op1;
   FiniteMPO B = String;
   FiniteMPO C = Op2;
   // Suppose the unit cell is 1 site, and all operators are 1-site local.  
   // We want to generate the string
   // X * ...  * X * A * B * ... * B * C * I * ....
   // This is the MPO 
   // ( X A 0 )
   // ( 0 B C )
   // ( 0 0 I )

   // If A and C are more then one unit cell, then set A = A1 * A2 * A3
   // C = C1 * C2 * C3
   // and the MPO is
   // ( X  A1 0  0  0  0  0  )
   // ( 0  0  A2 0  0  0  0  )
   // ( 0  0  0  A3 0  0  0  )
   // ( 0  0  0  B  C0 0  0  )
   // ( 0  0  0  0  0  C1 0  )
   // ( 0  0  0  0  0  0  C2 )
   // ( 0  0  0  0  0  0  I  )

   // Suppose now that we have 6-site operators A,C with a 3-site unit cell.  
   // Let A = A0 A1 A2 A3 A4 A5
   // Let C = C0 C1 C2 C3 C4 C5
   // ( X  A0 0  0  0  0  0  ) ( X  0  0  0  0  0  0 ) ( X  0  0  0  0  )
   // ( 0  0  A3 0  0  0  0  ) ( 0  A1 0  0  0  0  0 ) ( 0  A2 0  0  0  )
   // ( 0  0  0  B  C0 0  0  ) ( 0  0  A4 0  0  0  0 ) ( 0  0  A5 0  0  )
   // ( 0  0  0  0  0  C3 0  ) ( 0  0  0  B  0  0  0 ) ( 0  0  B  0  0  )
   // ( 0  0  0  0  0  0  I  ) ( 0  0  0  0  C1 0  0 ) ( 0  0  0  C2 0  )
   //                          ( 0  0  0  0  0  C4 0 ) ( 0  0  0  0  C5 )
   //                          ( 0  0  0  0  0  0  I ) ( 0  0  0  0  I  )
   // 
   // Upon coarse-graining this gives
   // ( XXX A0A1A2 0      0             )
   // ( 0   0      A3A4A5 0             )
   // ( 0   0      BBB    C0C1C2        )
   // ( 0   0      0      0      C3C4C5 )
   // ( 0   0      0      0      I      )

   CHECK(UnitCellSize % SiteList.size() == 0)
      ("UnitCellSize for sum_string() must be a multiple of UnitCell.size()");
   CHECK(A.size() % UnitCellSize == 0)

      ("Operator for sum_string() must be a multiple of the unit cell")
      (A.size())(UnitCellSize);
   CHECK(C.size() % UnitCellSize == 0)
      ("Operator for sum_string() must be a multiple of the unit cell")
      (C.size())(UnitCellSize);
   // THe JW string operator must have a unit cell that is a divisor of UnitCellSize
   CHECK(UnitCellSize % X.size() == 0)
      ("JW string for sum_string() must divide UnitCell.size()");
   CHECK(UnitCellSize % B.size() == 0)
      ("middle string for sum_string() must divide UnitCell.size()");

   // The aux basis for the operators is already fixed correctly
   CHECK_EQUAL(X.qn2(), A.qn1());
   CHECK_EQUAL(A.qn2(), B.qn1());
   CHECK_EQUAL(B.qn2(), C.qn1());
   CHECK(is_scalar(C.qn2()));

   // Split A and C into UnitCellSize pieces
   std::vector<std::vector<OperatorComponent>> SplitA = SplitOperator(A, UnitCellSize);
   std::vector<std::vector<OperatorComponent>> SplitC = SplitOperator(C, UnitCellSize);

   // and make the identity operator; C.qn2() is always the identity
   FiniteMPO Ident = identity_mpo(SiteList, C.qn2());

   TriangularMPO Result(UnitCellSize);
   for (int i = 0; i < UnitCellSize; ++i)
   {
      // Construct the basis
      BasisList Basis1(A.GetSymmetryList());
      BasisList Basis2(A.GetSymmetryList());
      if (i != 0)
	 JoinBasis(Basis1, JW[i%JW.size()].Basis1());
      JoinBasis(Basis2, JW[i%JW.size()].Basis2());

      for (unsigned n = 0; n < SplitA.size(); ++n)
      {
	 JoinBasis(Basis1, SplitA[n][i].Basis1());
	 JoinBasis(Basis2, SplitA[n][i].Basis2());
      }

      if (i != 0)
	 JoinBasis(Basis1, B[i%B.size()].Basis1());
      if (i != int(UnitCellSize-1))
	 JoinBasis(Basis2, B[i%B.size()].Basis2());

      for (unsigned n = 0; n < SplitC.size(); ++n)
      {
	 JoinBasis(Basis1, SplitC[n][i].Basis1());
	 JoinBasis(Basis2, SplitC[n][i].Basis2());
      }

      JoinBasis(Basis1, Ident[i%SiteList.size()].Basis1());
      if (i != int(UnitCellSize)-1)
	 JoinBasis(Basis2, Ident[i%SiteList.size()].Basis2());

      // Construct the OperatorComponent
      OperatorComponent Comp(A[i].LocalBasis1(), A[i].LocalBasis2(), Basis1, Basis2);

      // The JW goes in the top left
      int r = 0;

      int c = 0;
      SetComponents(Comp, X[i%X.size()], r, c);
      if (i != 0)
	 r += X[i%X.size()].Basis1().size();
      c += X[i%X.size()].Basis2().size();

      // the A components go along the diagonal
      for (unsigned n = 0; n < SplitA.size(); ++n)
      {
	 SetComponents(Comp, SplitA[n][i], r, c);
	 r += SplitA[n][i].Basis1().size();
	 c += SplitA[n][i].Basis2().size();
      }

      if (i == int(UnitCellSize-1))
	 --c;  // operator A is guaranteed to have only one column for the final entry

      // the B operator
      SetComponents(Comp, B[i%B.size()], r, c);
      if (i != 0)
	 r += B[i%B.size()].Basis1().size();
      //      if (i != int(UnitCellSize-1))
      c += B[i%B.size()].Basis2().size();

      // the C components continue on the diagonal
      for (unsigned n = 0; n < SplitC.size()-1; ++n)
      {
	 SetComponents(Comp, SplitC[n][i], r, c);
	 r += SplitC[n][i].Basis1().size();
	 c += SplitC[n][i].Basis2().size();
      }

      SetComponents(Comp, SplitC.back()[i], r, c);
      r += SplitC.back()[i].Basis1().size();
      if (i != int(UnitCellSize)-1)
	 c += SplitC.back()[i].Basis2().size();

      // The identity goes in the bottom right
      SetComponents(Comp, Ident[i%SiteList.size()], r, c);

      // check that we're at the end
      CHECK_EQUAL(r+Ident[i%SiteList.size()].Basis1().size(), Basis1.size());
      CHECK_EQUAL(c+Ident[i%SiteList.size()].Basis2().size(), Basis2.size());

      DEBUG_TRACE(Comp);
      
      Comp.debug_check_structure();

      Result[i] = Comp;
   }

   Result.check_structure();
   optimize(Result);
   return Result;
}

// This version of sum_string takes UnitCellMPO's for the operator arguments.  The String term
// must be a scalar with bosonic commutation, and cannot be any longer than UnitCellSize.
TriangularMPO sum_string(UnitCellMPO const& Op1_, UnitCellMPO const& String_, UnitCellMPO const& Op2_, 
			 int UnitCellSize,
			 QuantumNumbers::QuantumNumber q)
{
   CHECK(is_transform_target(Op1_.TransformsAs(), Op2_.TransformsAs(), q));
   CHECK(String_.size() <= UnitCellSize)("String operator cannot exceed the UnitCellSize in sum_string()");
   CHECK(is_scalar(String_.TransformsAs()));

   FiniteMPO Op1 = ExtendToCoverUnitCell(Op1_, UnitCellSize).MPO();
   FiniteMPO String = ExtendToCoverUnitCell(String_, UnitCellSize).MPO();
   FiniteMPO Op2 = ExtendToCoverUnitCell(Op2_, UnitCellSize).MPO();

   SiteListType SiteList = *Op1_.GetSiteList();

   CHECK(UnitCellSize % SiteList.size() == 0)
      ("UnitCellSize for sum_string must be a multiple of the lattice unit cell size");

   // shifting the quantum numbers in the aux basis
   FiniteMPO IdentShiftAB = identity_mpo(SiteList, Op2.qn1());
   String = String * repeat(IdentShiftAB, String.size() / SiteList.size());
   Op1 = project(Op1 * repeat(IdentShiftAB, Op1.size() / SiteList.size()), q);

   // Multiply through the JW string of Op2
   FiniteMPO JW2 = repeat(string_mpo(SiteList, Op2_.Commute().SignOperator(), QuantumNumber(Op2.GetSymmetryList())), 
			  UnitCellSize / SiteList.size());
   String = String * JW2;
   Op1 = Op1 * repeat(JW2, Op1.size() / JW2.size());

   // compute the overall JW string
   FiniteMPO JW = repeat(string_mpo(SiteList, Op1_.Commute().SignOperator(), Op1.qn1()),
			 UnitCellSize / SiteList.size()) * JW2;

   return sum_string(SiteList, JW, Op1, String, Op2, UnitCellSize, q);
}

TriangularMPO sum_string_dot(UnitCellMPO const& Op1_, UnitCellMPO const& String_, UnitCellMPO const& Op2_, 
			     int UnitCellSize)
{
   
   return sum_string(std::sqrt(degree(Op1_.TransformsAs())) * Op1_, String_, Op2_, UnitCellSize,
			 QuantumNumber(Op1_.GetSymmetryList()));
}

TriangularMPO sum_string_dot(UnitCellMPO const& Op1_, UnitCellMPO const& String_, UnitCellMPO const& Op2_)
{
   return sum_string_dot(Op1_, String_, Op2_, Op1_.GetSiteList()->size());
}

TriangularMPO sum_string_inner(UnitCellMPO const& Op1_, UnitCellMPO const& String_, UnitCellMPO const& Op2_, 
			     int UnitCellSize)
{
   
   return sum_string_dot(adjoint(Op1_), String_, Op2_, UnitCellSize);
}

TriangularMPO sum_string_inner(UnitCellMPO const& Op1_, UnitCellMPO const& String_, UnitCellMPO const& Op2_)
{
   return sum_string_dot(adjoint(Op1_), String_, Op2_, Op1_.GetSiteList()->size());
}
