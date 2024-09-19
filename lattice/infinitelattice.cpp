// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// lattice/infinitelattice.cpp
//
// Copyright (C) 2015-2023 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2024 Jesse Osborne <j.osborne@uqconnect.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#include "infinitelattice.h"
#include "lattice/infinite-parser.h"
#include "mp/copyright.h" // for EscapeArgument

PStream::VersionTag LatticeVersion(6);
// version 5 doesn't change anything in the InfiniteLattice directly, but
// adds versioning support for the InfiniteMPO structure
// version 6 doesn't change anything in the InfiniteLattice directly, but
// adds versioning support for the ProductMPO structure


InfiniteLattice::InfiniteLattice() : UnitCell_(new UnitCell()), OwnUnitCell_(true)
{
}

InfiniteLattice::InfiniteLattice(UnitCell* uc)
   : UnitCell_(uc), OwnUnitCell_(false)
{
}

InfiniteLattice::InfiniteLattice(UnitCell&& uc)
   : UnitCell_(new UnitCell(std::move(uc))), OwnUnitCell_(true)
{
}

InfiniteLattice::InfiniteLattice(std::string const& Description, UnitCell& uc)
   : Description_(Description), UnitCell_(&uc), OwnUnitCell_(false)
{
}


InfiniteLattice::InfiniteLattice(std::string const& Description, UnitCell&& uc)
   : Description_(Description), UnitCell_(new UnitCell(std::move(uc))), OwnUnitCell_(true)
{
}

InfiniteLattice::InfiniteLattice(InfiniteLattice const& Other)
   : Description_(Other.Description_),
     Authors_(Other.Authors_),
     CommandLine_(Other.CommandLine_),
     Timestamp_(Other.Timestamp_),
     UnitCell_(Other.OwnUnitCell_ ? new UnitCell(*Other.UnitCell_) : Other.UnitCell_),
     OwnUnitCell_(Other.OwnUnitCell_),
     Operators_(Other.Operators_),
     Arguments_(Other.Arguments_),
     Functions_(Other.Functions_)
{
}

InfiniteLattice::InfiniteLattice(InfiniteLattice&& Other)
   : Description_(std::move(Other.Description_)),
     Authors_(std::move(Other.Authors_)),
     CommandLine_(std::move(Other.CommandLine_)),
     Timestamp_(std::move(Other.Timestamp_)),
     UnitCell_(std::move(Other.UnitCell_)),
     OwnUnitCell_(std::move(Other.OwnUnitCell_)),
     Operators_(std::move(Other.Operators_)),
     Arguments_(std::move(Other.Arguments_)),
     Functions_(std::move(Other.Functions_))
{
   Other.UnitCell_ = nullptr;
   Other.OwnUnitCell_ = false;
}

InfiniteLattice&
InfiniteLattice::operator=(InfiniteLattice const& Other)
{
   Description_ = Other.Description_;
   Authors_ = Other.Authors_;
   CommandLine_ = Other.CommandLine_;
   Timestamp_ = Other.Timestamp_;
   UnitCell_ = Other.OwnUnitCell_ ? new UnitCell(*Other.UnitCell_) : Other.UnitCell_;
   OwnUnitCell_ = Other.OwnUnitCell_;
   Operators_ = Other.Operators_;
   Arguments_ = Other.Arguments_;
   Functions_ = Other.Functions_;
   return *this;
}

InfiniteLattice&
InfiniteLattice::operator=(InfiniteLattice&& Other)
{
   Description_ = std::move(Other.Description_);
   Authors_ = std::move(Other.Authors_);
   CommandLine_ = std::move(Other.CommandLine_);
   Timestamp_ = std::move(Other.Timestamp_);
   UnitCell_ = std::move(Other.UnitCell_);
   OwnUnitCell_ = std::move(Other.OwnUnitCell_);
   Operators_ = std::move(Other.Operators_);
   Arguments_ = std::move(Other.Arguments_);
   Functions_ = std::move(Other.Functions_);
   Other.UnitCell_ = nullptr;
   Other.OwnUnitCell_ = false;
   return *this;
}

InfiniteLattice::~InfiniteLattice()
{
   if (OwnUnitCell_)
      delete UnitCell_;
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
   // set the main description
   if (!Desc.description().empty())
      this->set_description(Desc.description());

   // Author information
   this->authors().insert(this->authors().end(), Desc.authors().begin(), Desc.authors().end());

   // iterate through the descriptions
   for (OperatorDescriptions::const_iterator I = Desc.begin(); I != Desc.end(); ++I)
   {
      if (this->operator_exists(std::get<0>(*I)))
      {
         // see if the operator is conditional
         if (std::get<3>(*I) && (!(*std::get<3>(*I))()))
         {
            std::cerr << "warning: conditional lattice operator " << std::get<0>(*I)
                      << " (conditional on: " << std::get<2>(*I) << ") should not be defined, but is!\n";
         }
         Operators_[std::get<0>(*I)].set_description(std::get<1>(*I));
      }
      else
      {
         // is the operator optional?
         if (!std::get<2>(*I).empty() || std::get<3>(*I))
         {
            // yes, check and see that we satisfy the condition
            if (std::get<3>(*I))
            {
               // invoke the function
               if (((*std::get<3>(*I))()))
               {
                  std::cerr << "warning: conditional lattice operator "  << std::get<0>(*I)
                            << " should be defined but is not.\n";
               }
            }
         }
         else
         {
            std::cerr << "warning: operator " << std::get<0>(*I)
                      << " has a description but is not defined in the lattice.\n";
         }
      }
   }

   // Now go through the operators and check for any that don't have a description, or are optional but
   // should not have been defined
   for (const_operator_iterator I = this->begin_operator(); I != this->end_operator(); ++I)
   {
      if (I->second.description().empty())
      {
         std::cerr << "warning: lattice operator " << I->first << " has no description.\n";
      }
   }

   // Same for the functions
   for (OperatorDescriptions::const_iterator I = Desc.begin_function(); I != Desc.end_function(); ++I)
   {
      if (this->function_exists(std::get<0>(*I)))
      {
         // see if the function is conditional
         if (std::get<3>(*I) && (!(*std::get<3>(*I))()))
         {
            std::cerr << "warning: conditional lattice function " << std::get<0>(*I)
                      << " (conditional on: " << std::get<2>(*I) << ") should not be defined, but is!\n";
         }
         Functions_[std::get<0>(*I)].set_description(std::get<1>(*I));
      }
      else
      {
         // is the operator optional?
         if (!std::get<2>(*I).empty() || std::get<3>(*I))
         {
            // yes, check and see that we satisfy the condition
            if (std::get<3>(*I))
            {
               // invoke the function
               if (((*std::get<3>(*I))()))
               {
                  std::cerr << "warning: conditional lattice function "  << std::get<0>(*I)
                            << " should be defined but is not.\n";
               }
            }
         }
         else
         {
            std::cerr << "warning: function " << std::get<0>(*I)
                      << " has a description but is not defined in the lattice.\n";
         }
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

   // set the main description
   if (!Desc.description().empty())
   {
      this->set_description(Desc.description());
   }

   // cell operator descriptions
   for (OperatorDescriptions::const_iterator I = Desc.cell_begin(); I != Desc.cell_end(); ++I)
   {
      if (this->GetUnitCell().operator_exists(std::get<0>(*I)))
      {
         // see if the operator is conditional
         if (std::get<3>(*I) && (!(*std::get<3>(*I))()))
         {
            std::cerr << "warning: conditional unit cell operator " << std::get<0>(*I)
                      << " (conditional on: " << std::get<2>(*I) << ") should not be defined, but is!\n";
         }
         this->GetUnitCell()[std::get<0>(*I)].set_description(std::get<1>(*I));
      }
      else
      {
         // is the operator optional?
         if (!std::get<2>(*I).empty() || std::get<3>(*I))
         {
            // yes, check and see that we satisfy the condition
            if (std::get<3>(*I))
            {
               // invoke the function
               if (((*std::get<3>(*I))()))
               {
                  std::cerr << "warning: conditional unit cell operator "  << std::get<0>(*I)
                            << "should be defined but is not.\n";
               }
            }
         }
         else
         {
            std::cerr << "warning: cell operator " << std::get<0>(*I)
                      << " has a description but is not defined in the unit cell.\n";
         }
      }
   }

   // Now go through the unit cell operators and check for any that don't have a description
   for (UnitCell::const_operator_iterator I = this->GetUnitCell().begin_operator();
        I != this->GetUnitCell().end_operator(); ++I)
   {
      if (std::get<1>(*I).description().empty())
      {
         std::cerr << "warning: cell operator " << I->first << " has no description.\n";
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

BasicTriangularMPO const&
InfiniteLattice::as_basic_triangular_mpo(std::string const& Op) const
{
   OperatorListType::const_iterator I = Operators_.find(Op);
   CHECK(I != Operators_.end())("Operator does not exist in the lattice!")(Op);
   return I->second.as_basic_triangular_mpo();
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
   Function::ArgumentList Args = GetArguments(Func.args(), Params,
                                              ParseInfiniteArgument(*this));
   return ParseInfiniteOperator(*this, Func.definition(), Args);
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
   out << L.Authors_;
   out << L.CommandLine_;
   out << L.Timestamp_;
   out << *L.UnitCell_;
   out << L.Operators_;
   out << L.Arguments_;
   out << L.Functions_;
   return out;
}

PStream::ipstream&
operator>>(PStream::ipstream& in, InfiniteLattice& L)
{
   PStream::VersionSentry Sentry(in, LatticeVersion, in.read<int>());

   if (Sentry.version() > 6)
   {
      PANIC("This program is too old to read this lattice file format,"
            "  Maximum readable version number is 6")(Sentry.version());
   } else if (Sentry.version() < 0)
   {
      PANIC("Lattice file is too old, please reconstruct the lattice."
            "  Expected Version >= 0")(Sentry.version());
   }

   in >> L.Description_;
   if (Sentry.version() >= 3)
   {
      in >> L.Authors_;
   }
   else
   {
      L.Authors_.clear();
   }
   in >> L.CommandLine_;
   in >> L.Timestamp_;
   if (L.OwnUnitCell_)
   {
      delete L.UnitCell_;
   }
   L.OwnUnitCell_ = true;
   L.UnitCell_ = new UnitCell();
   in >> (*L.UnitCell_);
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

BasicTriangularMPO make_triangular(UnitCellMPO const& Op)
{
   if (Op.is_null())
      return BasicTriangularMPO();

   CHECK_EQUAL(Op.Basis1().size(), 1)("Operator Basis1 is not one dimensional.  Basis must contain only one quantum number for sum_unit()");
   CHECK_EQUAL(Op.Basis2().size(), 1)("Operator Basis2 is not one dimensional.  Basis must contain only one quantum number for sum_unit()");

   return sum_unit(Op.GetJWStringUnit(), Op.MPO(), Op.size());
}

BasicTriangularMPO sum_unit(UnitCellMPO const& Op)
{
   if (Op.is_null())
      return BasicTriangularMPO();
   return sum_unit(Op, Op.GetSiteList()->size());
}

BasicTriangularMPO sum_unit(UnitCellMPO const& Op, int UnitCellSize)
{
   if (Op.is_null())
      return BasicTriangularMPO();

   CHECK_EQUAL(Op.Basis1().size(), 1)("Operator Basis1 is not one dimensional.  Basis must contain only one quantum number for sum_unit()");
   CHECK_EQUAL(Op.Basis2().size(), 1)("Operator Basis2 is not one dimensional.  Basis must contain only one quantum number for sum_unit()");

   CHECK(UnitCellSize % Op.coarse_grain_factor() == 0);

   return sum_unit(Op.GetJWStringUnit(), ExtendToCoverUnitCell(Op, UnitCellSize / Op.coarse_grain_factor()).MPO(),
                   UnitCellSize / Op.coarse_grain_factor());
}

std::vector<std::vector<OperatorComponent>>
   SplitOperator(BasicFiniteMPO const& Op, int UnitCellSize)
{
   DEBUG_CHECK_EQUAL(Op.size() % UnitCellSize, 0);
   std::vector<std::vector<OperatorComponent> > Result;
   BasicFiniteMPO::const_iterator I = Op.begin();
   BasicFiniteMPO::const_iterator J = I;
   while (J != Op.end())
   {
      std::advance(J, UnitCellSize);
      Result.push_back(std::vector<OperatorComponent>(I,J));
      I = J;
   }
   return Result;
}

BasicTriangularMPO sum_unit(BasicFiniteMPO const& JW2, BasicFiniteMPO const& Op2, int UnitCellSize)
{
   BasicFiniteMPO Op = Op2;
   optimize(Op);

   BasicFiniteMPO JW = JW2;
   optimize(JW);

   if (Op.is_null())
      return BasicTriangularMPO();
   CHECK(Op.is_irreducible());
   CHECK(Op.size() % UnitCellSize == 0)
      ("Operator for sum_unit() must be a multiple of the unit cell")
      (Op.size())(UnitCellSize);
   // THe JW string operator must have a unit cell that is a divisor of UnitCellSize
   CHECK(UnitCellSize % JW.size() == 0)
      ("JW string for sum_unit() must divide UnitCellSize");

   // Suppose that the unit cell size is 1.  Then if Op is A \times B \otimes C, then
   // the resulting BasicTriangularMPO is
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

   // To construct this operator, we firstly split Op into UnitCellSize pieces
   std::vector<std::vector<OperatorComponent>> SplitOp = SplitOperator(Op, UnitCellSize);

   BasicTriangularMPO Result(UnitCellSize);
   for (int i = 0; i < UnitCellSize; ++i)
   {
      OperatorComponent Ident = OperatorComponent::make_identity(Op[i].LocalBasis1());

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
      JoinBasis(Basis1, Ident.Basis1());
      if (i != int(UnitCellSize)-1)
         JoinBasis(Basis2, Ident.Basis2());

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
      SetComponents(C, Ident, r, c);

      // check that we're at the end
      CHECK_EQUAL(r+Ident.Basis1().size(), Basis1.size());
      CHECK_EQUAL(c+Ident.Basis2().size(), Basis2.size());

      C.debug_check_structure();

      Result[i] = C;
   }

   Result.check_structure();
   optimize(Result);
   return Result;
}

BasicTriangularMPO sum_kink(UnitCellMPO const& Kink, UnitCellMPO const& Op)
{
   if (Op.is_null() || Kink.is_null())
      return BasicTriangularMPO();

   CHECK_EQUAL(Kink.GetSiteList()->size(), Op.GetSiteList()->size())
      ("Operators for sum_kink must have the same unit cell!");

   return sum_kink(Kink, Op, Op.GetSiteList()->size());
}

BasicTriangularMPO sum_kink(UnitCellMPO const& Kink, UnitCellMPO const& Op, int UnitCellSize)
{
   if (Op.is_null() || Kink.is_null())
      return BasicTriangularMPO();

   CHECK(UnitCellSize % Op.coarse_grain_factor() == 0);

   UnitCellMPO Kink_ = ExtendToCoverUnitCell(Kink, UnitCellSize / Op.coarse_grain_factor());
   UnitCellMPO Op_ = ExtendToCoverUnitCell(Op, UnitCellSize / Op.coarse_grain_factor());

   BasicFiniteMPO Ident = repeat(Op_.GetJWStringUnit(), Kink_.size() / (Op_.GetSiteList()->size() / Op_.coarse_grain_factor()));
   return sum_unit(Kink_.MPO()*Ident, Op_.MPO(), UnitCellSize / Op.coarse_grain_factor());
}

BasicTriangularMPO sum_k(std::complex<double> const& k, UnitCellMPO const& Op)
{
   if (Op.is_null())
      return BasicTriangularMPO();

   return sum_k(k, Op, Op.GetSiteList()->size());
}

BasicTriangularMPO sum_k(std::complex<double> const& k, UnitCellMPO const& Op, int UnitCellSize)
{
   if (Op.is_null())
      return BasicTriangularMPO();

   CHECK(UnitCellSize % Op.coarse_grain_factor() == 0);

   return sum_unit(exp(std::complex<double>(0,1)*k) * Op.GetJWStringUnit(),
                   ExtendToCoverUnitCell(Op, UnitCellSize / Op.coarse_grain_factor()).MPO(),
                   UnitCellSize / Op.coarse_grain_factor());
}

BasicTriangularMPO make_zero(SiteListType const& SiteList)
{
   if (SiteList.empty())
      return BasicTriangularMPO();

   // The auxiliary basis is the same at every site
   // Construct the basis
   BasisList b(SiteList[0].GetSymmetryList());
   QuantumNumbers::QuantumNumber Ident(SiteList[0].GetSymmetryList());
   b.push_back(Ident);
   b.push_back(Ident);

   BasicTriangularMPO Result(SiteList.size());
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

BasicTriangularMPO sum_string(SiteListType const& SiteList, BasicFiniteMPO const& JW, BasicFiniteMPO const& Op1,
                         BasicFiniteMPO const& String, BasicFiniteMPO const& Op2, int UnitCellSize,
                         QuantumNumbers::QuantumNumber q)
{
   BasicFiniteMPO X = JW;
   BasicFiniteMPO A = Op1;
   BasicFiniteMPO B = String;
   BasicFiniteMPO C = Op2;
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
   BasicFiniteMPO Ident = identity_mpo(SiteList, C.qn2());

   BasicTriangularMPO Result(UnitCellSize);
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
BasicTriangularMPO sum_string(UnitCellMPO const& Op1_, UnitCellMPO const& String_, UnitCellMPO const& Op2_,
                         int UnitCellSize,
                         QuantumNumbers::QuantumNumber q)
{
   CHECK(is_transform_target(Op1_.TransformsAs(), Op2_.TransformsAs(), q));
   CHECK(String_.size() <= UnitCellSize)("String operator cannot exceed the UnitCellSize in sum_string()");
   CHECK(is_scalar(String_.TransformsAs()));

   BasicFiniteMPO Op1 = ExtendToCoverUnitCell(Op1_, UnitCellSize).MPO();
   BasicFiniteMPO String = ExtendToCoverUnitCell(String_, UnitCellSize).MPO();
   BasicFiniteMPO Op2 = ExtendToCoverUnitCell(Op2_, UnitCellSize).MPO();

   SiteListType SiteList = *Op1_.GetSiteList();

   CHECK(UnitCellSize % SiteList.size() == 0)
      ("UnitCellSize for sum_string must be a multiple of the lattice unit cell size");

   // shifting the quantum numbers in the aux basis
   BasicFiniteMPO IdentShiftAB = identity_mpo(SiteList, Op2.qn1());
   String = String * repeat(IdentShiftAB, String.size() / SiteList.size());
   Op1 = project(Op1 * repeat(IdentShiftAB, Op1.size() / SiteList.size()), q);

   // Multiply through the JW string of Op2
   BasicFiniteMPO JW2 = repeat(string_mpo(SiteList, Op2_.Commute().SignOperator(), QuantumNumber(Op2.GetSymmetryList())),
                          UnitCellSize / SiteList.size());
   String = String * JW2;
   Op1 = Op1 * repeat(JW2, Op1.size() / JW2.size());

   // compute the overall JW string
   BasicFiniteMPO JW = repeat(string_mpo(SiteList, Op1_.Commute().SignOperator(), Op1.qn1()),
                         UnitCellSize / SiteList.size()) * JW2;

   return sum_string(SiteList, JW, Op1, String, Op2, UnitCellSize, q);
}

BasicTriangularMPO sum_string(UnitCellMPO const& Op1_, UnitCellMPO const& String_, UnitCellMPO const& Op2_)
{
   return sum_string(Op1_, String_, Op2_, Op1_.GetSiteList()->size(), QuantumNumber(Op1_.GetSymmetryList()));
}

BasicTriangularMPO sum_string_dot(UnitCellMPO const& Op1_, UnitCellMPO const& String_, UnitCellMPO const& Op2_,
                             int UnitCellSize)
{

   return sum_string(std::sqrt(degree(Op1_.TransformsAs())) * Op1_, String_, Op2_, UnitCellSize,
                         QuantumNumber(Op1_.GetSymmetryList()));
}

BasicTriangularMPO sum_string_dot(UnitCellMPO const& Op1_, UnitCellMPO const& String_, UnitCellMPO const& Op2_)
{
   return sum_string_dot(Op1_, String_, Op2_, Op1_.GetSiteList()->size());
}

BasicTriangularMPO sum_string_inner(UnitCellMPO const& Op1_, UnitCellMPO const& String_, UnitCellMPO const& Op2_,
                             int UnitCellSize)
{

   return sum_string_dot(adjoint(Op1_), String_, Op2_, UnitCellSize);
}

BasicTriangularMPO sum_string_inner(UnitCellMPO const& Op1_, UnitCellMPO const& String_, UnitCellMPO const& Op2_)
{
   return sum_string_dot(adjoint(Op1_), String_, Op2_, Op1_.GetSiteList()->size());
}

/* **DEPRECATED**
ProductMPO prod_unit(UnitCellMPO const& Op_)
{
   return prod_unit_left_to_right(Op_.MPO(), Op_.GetSiteList()->size());
}

ProductMPO prod_unit(UnitCellMPO const& Op_, int UnitCellSize)
{
   return prod_unit_left_to_right(Op_.MPO(), UnitCellSize);
}
*/

// 5-term version of sum_string
BasicTriangularMPO sum_string(SiteListType const& SiteList, BasicFiniteMPO const& JW, BasicFiniteMPO const& Op1,
                         BasicFiniteMPO const& String1, BasicFiniteMPO const& Op2,
                         BasicFiniteMPO const& String2, BasicFiniteMPO const& Op3,
                         int UnitCellSize, QuantumNumbers::QuantumNumber q)
{
   BasicFiniteMPO X = JW;
   BasicFiniteMPO A = Op1;
   BasicFiniteMPO B = String1;
   BasicFiniteMPO C = Op2;
   BasicFiniteMPO D = String2;
   BasicFiniteMPO E = Op3;
   // Suppose the unit cell is 1 site, and all operators are 1-site local.
   // We want to generate the string
   // X * ...  * X * A * B * ... * B * C * D * ... * D * E * I * ....
   // This is the MPO
   // ( X A 0 0 )
   // ( 0 B C 0 )
   // ( 0 0 D E )
   // ( 0 0 0 I )

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
   CHECK(D.size() % UnitCellSize == 0)
      ("Operator for sum_string() must be a multiple of the unit cell")
      (D.size())(UnitCellSize);
   CHECK(E.size() % UnitCellSize == 0)
      ("Operator for sum_string() must be a multiple of the unit cell")
      (E.size())(UnitCellSize);

   // The aux basis for the operators is already fixed correctly
   CHECK_EQUAL(X.qn2(), A.qn1());
   CHECK_EQUAL(A.qn2(), B.qn1());
   CHECK_EQUAL(B.qn2(), C.qn1());
   CHECK_EQUAL(C.qn2(), D.qn1());
   CHECK_EQUAL(D.qn2(), E.qn1());
   CHECK(is_scalar(E.qn2()));

   // Split A, C, E into UnitCellSize pieces
   std::vector<std::vector<OperatorComponent>> SplitA = SplitOperator(A, UnitCellSize);
   std::vector<std::vector<OperatorComponent>> SplitC = SplitOperator(C, UnitCellSize);
   std::vector<std::vector<OperatorComponent>> SplitE = SplitOperator(E, UnitCellSize);

   // and make the identity operator; C.qn2() is always the identity
   BasicFiniteMPO Ident = identity_mpo(SiteList, C.qn2());

   BasicTriangularMPO Result(UnitCellSize);
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

      if (i != 0)
         JoinBasis(Basis1, D[i%D.size()].Basis1());
      if (i != int(UnitCellSize-1))
         JoinBasis(Basis2, D[i%D.size()].Basis2());

      for (unsigned n = 0; n < SplitE.size(); ++n)
      {
         JoinBasis(Basis1, SplitE[n][i].Basis1());
         JoinBasis(Basis2, SplitE[n][i].Basis2());
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

      // The D operator

      SetComponents(Comp, D[i%D.size()], r, c);
      if (i != 0)
         r += D[i%D.size()].Basis1().size();
      //      if (i != int(UnitCellSize-1))
      c += D[i%D.size()].Basis2().size();

      // the E components continue on the diagonal
      for (unsigned n = 0; n < SplitE.size()-1; ++n)
      {
         SetComponents(Comp, SplitE[n][i], r, c);
         r += SplitE[n][i].Basis1().size();
         c += SplitE[n][i].Basis2().size();
      }

      SetComponents(Comp, SplitE.back()[i], r, c);
      r += SplitE.back()[i].Basis1().size();
      if (i != int(UnitCellSize)-1)
         c += SplitE.back()[i].Basis2().size();

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
BasicTriangularMPO sum_string(UnitCellMPO const& Op1_, UnitCellMPO const& String1_, UnitCellMPO const& Op2_,
                              UnitCellMPO const& String2_, UnitCellMPO const& Op3_,
                              int UnitCellSize,
                              QuantumNumbers::QuantumNumber q)
{
   CHECK(is_transform_target(Op1_.TransformsAs(), Op2_.TransformsAs(), q));
   CHECK(String1_.size() <= UnitCellSize)("String operator 1 cannot exceed the UnitCellSize in sum_string()");
   CHECK(String2_.size() <= UnitCellSize)("String operator 2 cannot exceed the UnitCellSize in sum_string()");
   CHECK(is_scalar(String1_.TransformsAs()));
   CHECK(is_scalar(String2_.TransformsAs()));

   BasicFiniteMPO Op1 = ExtendToCoverUnitCell(Op1_, UnitCellSize).MPO();
   BasicFiniteMPO String1 = ExtendToCoverUnitCell(String1_, UnitCellSize).MPO();
   BasicFiniteMPO Op2 = ExtendToCoverUnitCell(Op2_, UnitCellSize).MPO();
   BasicFiniteMPO String2 = ExtendToCoverUnitCell(String2_, UnitCellSize).MPO();
   BasicFiniteMPO Op3 = ExtendToCoverUnitCell(Op3_, UnitCellSize).MPO();

   SiteListType SiteList = *Op1_.GetSiteList();

   CHECK(UnitCellSize % SiteList.size() == 0)
      ("UnitCellSize for sum_string must be a multiple of the lattice unit cell size");

   // shifting the quantum numbers in the aux basis for Op3
   BasicFiniteMPO IdentShift = identity_mpo(SiteList, Op3.qn1());
   String2 = String2 * repeat(IdentShift, String2.size() / SiteList.size());
   Op2 = Op2 * repeat(IdentShift, Op2.size() / SiteList.size());
   IdentShift = IdentShift * identity_mpo(SiteList, Op2.qn1());
   String1 = String1 * repeat(IdentShift, String2.size() / SiteList.size());
   Op1 = project(Op1 * repeat(IdentShift, Op1.size() / SiteList.size()), q);

   // Multiply through the JW string of Op3
   BasicFiniteMPO JW3 = repeat(string_mpo(SiteList, Op3_.Commute().SignOperator(), QuantumNumber(Op3.GetSymmetryList())),
                          UnitCellSize / SiteList.size());
   String2 = String2 * JW3;
   Op2 = Op2 * repeat(JW3, Op2.size() / JW3.size());
   String2 = String2 * JW3;
   Op1 = Op1 * repeat(JW3, Op1.size() / JW3.size());

   // Multiply through the JW string of Op2
   BasicFiniteMPO JW2 = repeat(string_mpo(SiteList, Op2_.Commute().SignOperator(), QuantumNumber(Op2.GetSymmetryList())),
                          UnitCellSize / SiteList.size());
   String1 = String1 * JW2;
   Op1 = Op1 * repeat(JW2, Op1.size() / JW2.size());

   // compute the overall JW string
   BasicFiniteMPO JW = repeat(string_mpo(SiteList, Op1_.Commute().SignOperator(), Op1.qn1()),
                         UnitCellSize / SiteList.size()) * JW2 * JW3;

   return sum_string(SiteList, JW, Op1, String1, Op2, String2, Op3, UnitCellSize, q);
}

BasicTriangularMPO sum_string(UnitCellMPO const& Op1_, UnitCellMPO const& String1_, UnitCellMPO const& Op2_,
                              UnitCellMPO const& String2_, UnitCellMPO const& Op3_)
{
   return sum_string(Op1_, String1_, Op2_, String2_, Op3_, Op1_.GetSiteList()->size(), QuantumNumber(Op1_.GetSymmetryList()));
}

BasicTriangularMPO sum_partial(SiteListType const& SiteList, BasicTriangularMPO const& Op, BasicFiniteMPO const& Pivot, /* BasicFiniteMPO const& JW2, */ int UnitCellSize)
{
   //   CHECK_EQUAL(Pivot.size(), UnitCellSize);

   //BasicFiniteMPO JW = JW2;

   // Suppose that the unit cell size is 1.  Then we take the triangular MPO and form the sum
   // of all partial sums of it.
   // Eg, suppose the original triangular operator is 3x3. Then the resulting BasicTriangularMPO is
   // 4xr, and
   // ( X X X 0 )
   // ( 0 X X 0 )
   // ( 0 0 X P )
   // ( 0 0 0 I )
   // Now what happens if the JW string of P is non-trivial?  For the time being, we assume that
   // there is no JW string associated with P.

   // If the unit cell of P is bigger then 1 site, then we are summing over a larger unit cell.  In this case,
   // we take a unit cell of X, followed by a unit cell of P.  Eg, suppose we have a unit cell of 3 sites.
   // Then
   // (X X X 0  0) (X X X 0  0) (X X X 0)
   // (0 X X 0  0) (0 X X 0  0) (0 X X 0)
   // (0 0 X PP 0) (0 0 X 0  0) (0 0 X 0)
   // (0 0 0 0  I) (0 0 0 PP 0) (0 0 0 P)
   //              (0 0 0 PP 0) (0 0 0 P)
   //              (0 0 0 0  I) (0 0 0 I)
   // where we show the shape of P if it has bond dimensions 1x2x1
   //
   // Now it can also be the case that the pivot operator is multiple unit cells.  eg suppose it is
   // PPPQQQ.  Then we have
   // (X X X 0  0  0) (X X X 0  0  0) (X X X 0  0)
   // (0 X X 0  0  0) (0 X X 0  0  0) (0 X X 0  0)
   // (0 0 X PP 0  0) (0 0 X 0  0  0) (0 0 X 0  0)
   // (0 0 0 0  QQ 0) (0 0 0 PP 0  0) (0 0 0 PP 0)
   // (0 0 0 0  QQ 0) (0 0 0 PP 0  0) (0 0 0 PP 0)
   // (0 0 0 0  0  I) (0 0 0 0  QQ 0) (0 0 0 0  Q)
   //                 (0 0 0 0  QQ 0) (0 0 0 0  Q)
   //                 (0 0 0 0  0  I) (0 0 0 0  I)

   DEBUG_TRACE(Op.size())(UnitCellSize)(SiteList.size());

   // To construct this operator, we firstly split Pivot into UnitCellSize pieces
   std::vector<std::vector<OperatorComponent>> SplitPivot = SplitOperator(Pivot, UnitCellSize);

   // and we need the identity operator
   BasicFiniteMPO Ident = identity_mpo(SiteList, Pivot.qn2());

   BasicTriangularMPO Result(UnitCellSize);
   for (int i = 0; i < UnitCellSize; ++i)
   {
      // Construct the basis.  Start from the triangular operator
      BasisList Basis1(Op.Basis1());
      BasisList Basis2(Op.Basis2());

      // Now the pivot operator
      for (unsigned n = 0; n < SplitPivot.size(); ++n)
      {
         if (i != 0 || n != 0)
            JoinBasis(Basis1, SplitPivot[n][i].Basis1());
         JoinBasis(Basis2, SplitPivot[n][i].Basis2());
      }

      // And finally the identity
      JoinBasis(Basis1, Ident[i%SiteList.size()].Basis1());
      if (i != int(UnitCellSize)-1)
         JoinBasis(Basis2, Ident[i%SiteList.size()].Basis2());

      DEBUG_TRACE(Basis1)(Basis2)(SplitPivot.size());

      // Construct the OperatorComponent
      OperatorComponent C(Op[i%Op.size()].LocalBasis1(), Op[i%Op.size()].LocalBasis2(), Basis1, Basis2);

      int r = 0, c = 0;
      // The triangular MPO goes in the top left
      SetComponents(C, Op[i%Op.size()], r, c);
      r += Op[i%Op.size()].Basis1().size();
      c += Op[i%Op.size()].Basis2().size();

      // at the first entry, move back a row.
      if (i == 0)
	     --r;

      // Now the pivot operator
      for (unsigned n = 0; n < SplitPivot.size(); ++n)
      {
         SetComponents(C, SplitPivot[n][i], r, c);
         r += SplitPivot[n][i].Basis1().size();
         c += SplitPivot[n][i].Basis2().size();
      }

      // At the last entry, move back a column
      if (i == int(UnitCellSize)-1)
	     --c;

      // Finally, the identity in the bottom right
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

BasicTriangularMPO sum_partial(BasicTriangularMPO const& Op, UnitCellMPO const& Pivot, int UnitCellSize)
{
   return sum_partial(*Pivot.GetSiteList(), Op,
		      ExtendToCoverUnitCell(Pivot, UnitCellSize).MPO(),
		      UnitCellSize);
}

BasicTriangularMPO sum_partial(BasicTriangularMPO const& Op, UnitCellMPO const& Pivot)
{
   return sum_partial(Op, Pivot, Op.size());
}
