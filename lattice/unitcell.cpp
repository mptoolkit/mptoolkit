// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// lattice/unitcell.cpp
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

#include "unitcell.h"
#include "unitcell-parser.h"
#include "siteoperator-parser.h"

//
// UnitCell members
//

UnitCell::UnitCell()
{
}

UnitCell::UnitCell(UnitCell const& Other)
   : Sites(Other.Sites), Operators(Other.Operators),
     Arguments(Other.Arguments), Functions(Other.Functions)
{
}

UnitCell::UnitCell(LatticeSite const& s)
   : Sites(new SiteListType(1, s)),
     Arguments(s.begin_arg(), s.end_arg())
{
   this->SetDefaultOperators();
   // convert the operators to the UnitCell equivalents
   for (LatticeSite::const_operator_iterator I = s.begin_operator(); I != s.end_operator(); ++I)
   {
      Operators[I->first] = this->local_operator(I->first, 0, 0);
   }
   // Note, we cannot import the functions from the LatticeSite, because they generally won't
   // be valid UnitCell operator expressions (ie, they won't include the cell index)
}

UnitCell::UnitCell(LatticeSite const& s, LatticeSite const& t)
   : Sites(new SiteListType(1, s))
{
   {
      pvalue_lock<SiteListType> Lock(Sites);
      Lock->push_back(t);
   }
   this->SetDefaultOperators();
   this->check_structure();
}

UnitCell::UnitCell(LatticeSite const& s, LatticeSite const& t, LatticeSite const& u)
   : Sites(new SiteListType(1, s))
{
   {
      pvalue_lock<SiteListType> Lock(Sites);
      Lock->push_back(t);
      Lock->push_back(u);
   }
   this->SetDefaultOperators();
   this->check_structure();
}

UnitCell::UnitCell(LatticeSite const& s, LatticeSite const& t, LatticeSite const& u, LatticeSite const& v)
   : Sites(new SiteListType(1, s))
{
   {
      pvalue_lock<SiteListType> Lock(Sites);
      Lock->push_back(t);
      Lock->push_back(u);
      Lock->push_back(v);
   }
   this->SetDefaultOperators();
   this->check_structure();
}

UnitCell::UnitCell(std::vector<LatticeSite> const& s)
   : Sites(new SiteListType(s))
{
   this->SetDefaultOperators();
   this->check_structure();
}

UnitCell::UnitCell(SymmetryList const& sl, LatticeSite const& s)
   : Sites(new SiteListType(1, CoerceSymmetryList(s,sl)))
{
   this->SetDefaultOperators();
   this->check_structure();
}

UnitCell::UnitCell(SymmetryList const& sl, LatticeSite const& s, LatticeSite const& t)
   : Sites(new SiteListType(1, CoerceSymmetryList(s, sl)))
{
   {
      pvalue_lock<SiteListType> Lock(Sites);
      Lock->push_back(CoerceSymmetryList(t, sl));
   }
   this->SetDefaultOperators();
   this->check_structure();
}

UnitCell::UnitCell(SymmetryList const& sl, LatticeSite const& s, LatticeSite const& t, LatticeSite const& u)
   : Sites(new SiteListType(1, CoerceSymmetryList(s, sl)))
{
   {
      pvalue_lock<SiteListType> Lock(Sites);
      Lock->push_back(CoerceSymmetryList(t, sl));
      Lock->push_back(CoerceSymmetryList(u, sl));
   }
   this->SetDefaultOperators();
   this->check_structure();
}

UnitCell::UnitCell(SymmetryList const& sl, LatticeSite const& s, LatticeSite const& t,
                   LatticeSite const& u, LatticeSite const& v)
   : Sites(new SiteListType(1, CoerceSymmetryList(s, sl)))
{
   {
      pvalue_lock<SiteListType> Lock(Sites);
      Lock->push_back(CoerceSymmetryList(t, sl));
      Lock->push_back(CoerceSymmetryList(u, sl));
      Lock->push_back(CoerceSymmetryList(v, sl));
   }
   this->SetDefaultOperators();
   this->check_structure();
}

UnitCell::UnitCell(int RepeatCount, UnitCell const& l)
   : Sites(new SiteListType())
{
   {
      pvalue_lock<SiteListType> Lock(Sites);
      Lock->reserve(l.size()*RepeatCount);
      for (int i = 0; i < RepeatCount; ++i)
      {
         Lock->insert(Lock->end(), l.Sites->begin(), l.Sites->end());
      }
   }
   if (RepeatCount == 1)
   {
      Operators = l.Operators;
      Arguments = l.Arguments;
      Functions = l.Functions;
   }
   else
   {
      this->SetDefaultOperators();
   }
   this->check_structure();
}

UnitCell::UnitCell(UnitCell const& x1, UnitCell const& x2)
   : Sites(x1.Sites)
{
   {
      pvalue_lock<SiteListType> Lock(Sites);
      Lock->insert(Lock->end(), x2.Sites->begin(), x2.Sites->end());
   }
   this->SetDefaultOperators();
   this->check_structure();
}

UnitCell::UnitCell(UnitCell const& x1, UnitCell const& x2, UnitCell const& x3)
   : Sites(x1.Sites)
{
   {
      pvalue_lock<SiteListType> Lock(Sites);
      Lock->insert(Lock->end(), x2.Sites->begin(), x2.Sites->end());
      Lock->insert(Lock->end(), x3.Sites->begin(), x3.Sites->end());
   }
   this->SetDefaultOperators();
   this->check_structure();
}

UnitCell::UnitCell(UnitCell const& x1, UnitCell const& x2, UnitCell const& x3, UnitCell const& x4)
   : Sites(x1.Sites)
{
   {
      pvalue_lock<SiteListType> Lock(Sites);
      Lock->insert(Lock->end(), x2.Sites->begin(), x2.Sites->end());
      Lock->insert(Lock->end(), x3.Sites->begin(), x3.Sites->end());
      Lock->insert(Lock->end(), x4.Sites->begin(), x4.Sites->end());
   }
   this->SetDefaultOperators();
   this->check_structure();
}

UnitCell::UnitCell(int Size, LatticeSite const& s)
   : Sites(new SiteListType(Size, s))
{
   this->SetDefaultOperators();
   this->check_structure();
}

UnitCell&
UnitCell::operator=(UnitCell const& Other)
{
   Sites = Other.Sites;
   Operators = Other.Operators;
   Arguments = Other.Arguments;
   Functions = Other.Functions;
   return *this;
}

LatticeSite const&
UnitCell::operator[](int n) const
{
   return (*Sites)[n];
}

void
UnitCell::assign_operator(std::string const& Name, operator_type Op, int Offset)
{
   Operators[Name] = translate(std::move(Op), -Offset * int(this->size()));
}

bool
UnitCell::operator_exists(std::string const& s) const
{
   if (Operators.find(s) != Operators.end())
      return true;
   if (Sites->size() == 1 && Sites->front().operator_exists(s))
      return true;
   return false;
}

UnitCell::operator_type
UnitCell::operator[](std::string const& Op) const
{
   OperatorListType::const_iterator I = Operators.find(Op);
   if (I != Operators.end())
      return I->second;

   if (Sites->size() != 1)
   {
      PANIC("Operator not found in unit cell")(Op);
   }

   return this->local_operator(Op, 0);
}

UnitCell::operator_type
UnitCell::identity() const
{
   return this->operator[]("I");
}

UnitCell::operator_type&
UnitCell::operator[](std::string const& Op)
{
   OperatorListType::iterator I = Operators.find(Op);
   if (I != Operators.end())
      return I->second;
   // else if the unit cell is 1 site, look for a local operator
   // and copy it into our Operators.  We need to do this so that
   // we have a local copy and can return by reference (and it can be subsequently modified)
   if (Sites->size() == 1 && (*Sites)[0].operator_exists(Op))
   {
      Operators[Op] = this->local_operator(Op, 0);
   }
   return Operators[Op];
}

UnitCell::operator_type
UnitCell::operator()(std::string const& Op, int n) const
{
   OperatorListType::const_iterator I = Operators.find(Op);
   if (I != Operators.end())
   {
      return UnitCellMPO(Sites, I->second.MPO(), I->second.Commute(), I->second.offset()+n*this->size());
   }
   // else
   if (Sites->size() != 1)
   {
      PANIC("Operator not found in unit cell")(Op);
   }

   return this->local_operator(Op, n, 0);
}


bool
UnitCell::local_operator_exists(std::string const& Op, int n) const
{
   if (n < 0 || n >= int(Sites->size()))
      return false;
   return (*Sites)[n].operator_exists(Op);
}

LatticeCommute
UnitCell::Commute(std::string const& Op, int n) const
{
   CHECK(0 <= n && n < int(Sites->size()))("Site index is out of range")(n)(Sites->size());
   return (*Sites)[n][Op].Commute();
}

UnitCell::operator_type
UnitCell::local_operator(std::string const& Op, int Cell, int n) const
{
   CHECK(0 <= n && n < int(Sites->size()))("Site index is out of range")(n)(Sites->size());
   CHECK(this->local_operator_exists(Op, n));
   SiteOperator Operator = (*Sites)[n][Op];
   return this->map_local_operator(Operator, Cell, n);
}

UnitCell::operator_type
UnitCell::map_local_operator(SiteOperator const& Operator, int Cell, int n) const
{
   std::string SignOperator = Operator.Commute().SignOperator();

   FiniteMPO Result(Sites->size());

   BasisList Vacuum = make_vacuum_basis(Operator.GetSymmetryList());
   BasisList Basis = make_single_basis(Operator.TransformsAs());

   // Assemble the JW-string.
   for (int i = 0; i < n; ++i)
   {
      if (!(*Sites)[i].operator_exists(SignOperator))
      {
         WARNING("JW-string operator doesn't exist at a lattice site, using the identity")
            (i)(SignOperator);
      }
      SimpleOperator Op = (*Sites)[i].operator_exists(SignOperator) ?
         (*Sites)[i][SignOperator] : (*Sites)[i].identity();
      Result[i] = OperatorComponent(Op.Basis1(), Op.Basis2(), Basis, Basis);
      Result[i](0,0) = Op;
   }
   Result[n] = OperatorComponent(Operator.Basis1(), Operator.Basis2(), Basis, Vacuum);
   Result[n](0,0) = Operator;
   for (int i = n+1; i < this->size(); ++i)
   {
      SimpleOperator I = (*Sites)[i].identity();
      Result[i] = OperatorComponent(I.Basis1(), I.Basis2(), Vacuum, Vacuum);
      Result[i](0,0) = I;
   }

   return UnitCellMPO(Sites, Result, Operator.Commute(), Cell*this->size(),
                      Operator.description());
}

UnitCell::operator_type
UnitCell::local_operator(std::string const& Op, int n) const
{
   return this->local_operator(Op, 0, n);
}

UnitCell::function_type
UnitCell::func(std::string const& s) const
{
   const_function_iterator I = this->find_function(s);
   if (I != this->end_function())
      return I->second;
   // else
   return function_type();
}

// functor to parse default arguments of a UnitCell operator
struct ParseUnitCellExpression
{
   ParseUnitCellExpression(UnitCell const& Cell_) : Cell(Cell_) {}

   std::complex<double> operator()(Function::ArgumentList const& Args,
                                   std::string const& Str) const
   {
      return ParseUnitCellNumber(Cell, 0, Str, Args);
   }

   UnitCell const& Cell;
};



boost::variant<UnitCell::operator_type, std::complex<double> >
UnitCell::eval_function(Function::OperatorFunction const& Func, int Cell,
                        Function::ParameterList const& Params) const
{
   Function::ArgumentList Args = GetArguments(Func.args(), Params, ParseUnitCellExpression(*this));
   boost::variant<UnitCell::operator_type, std::complex<double> > Result
      = ParseUnitCellElement(*this, 0, Func.definition(), Args);
   if (boost::get<UnitCell::operator_type>(&Result))
   {
      boost::get<UnitCell::operator_type>(&Result)->translate(Cell*this->size());
   }
   return Result;
}

boost::variant<UnitCell::operator_type, std::complex<double> >
UnitCell::eval_function(std::string const& Func, int Cell,
                        Function::ParameterList const& Params) const
{
   if (this->function_exists(Func))
      return this->eval_function(this->func(Func), Cell, Params);
   // else
   return this->eval_local_function(Func, Cell, 0, Params);
}

boost::variant<UnitCell::operator_type, std::complex<double> >
UnitCell::eval_local_function(std::string const& Func, int Cell, int Site,
                              Function::ParameterList const& Params) const
{
   SiteElementType Element = Sites->operator[](Site).eval_function(Func, Params);

   // if the result is a c-number, we can return it without further processing
   std::complex<double>* c = boost::get<std::complex<double> >(&Element);
   if (c)
   {
      return *c;
   }
   // else we have to convert the SiteOperator into a UnitCellMPO
   SiteOperator Op = boost::get<SiteOperator>(Element);
   return this->map_local_operator(Op, Cell, Site);
}

UnitCellMPO
UnitCell::swap_gate(int Cell_i, int i, int Cell_j, int j) const
{
   // count the number of fermion swaps.  This is the number of fermions in
   // site i + fermions in site j, multiplied by the number of fermions between them.
   // And of course we need to swap site j with site i, so add N_i * N_j.
   // So in summary, the number of fermion swaps is
   // (N_i + N_j)*(sum_{k=i+1}^{j-1} N_k) + N_i*N_j

   // Now we need to turn that into something that can use the local P operators only.
   // Note that P_i = (-1)^N_i.  Therefore, the number of fermions on site i (mod 2) is
   // 0.5*(1-P_i).  Call this F_i.  We can replace all of the N_i above with F_i.
   // The F_i on each site commute, so we can expand out the exponential.
   // Hence we can write this as (-1)^{....} =
   //   (-1)^{F_i F_{i+1}                    }
   // * (-1)^{F_i         F_{i+2}            }
   // * ......
   // * (-1)^{F_i                 F_{j-1}    }

   // * (-1)^{    F_{i+1}                 F_j}
   // * ......
   // * (-1)^{                    F_{j-1} F_j}

   // * (-1)^{F_i                         F_j}

   // Now F_k is diagonal and eigenvalues 0,1, so the product is also eigenvalues 0,1.

   // maybe easier to do products of nearest-neighbor swaps?
   // If the swap is nearest-neighbor then the number of fermions is F_i * F_j.

   // Need to write it as a function of the individual basis states in F_i and F_j,
   // and write it as a modifier on the sign of the swap gate between those states.
   // if they are both fermions or neither are fermions then there is no change.
   // If there is one fermion, then insert the product of parity operators.
   // This only works if F_i and F_j are diagonal.

   CHECK(i >= 0 && i < this->size())("Site index is outside the unit cell!")(i)(this->size());
   CHECK(j >= 0 && j < this->size())("Site index is outside the unit cell!")(j)(this->size());
   // normal-order the sites
   if (Cell_i > Cell_j || (Cell_i == Cell_j && i > j))
   {
      std::swap(Cell_i, Cell_j);
      std::swap(i,j);
   }

   if (Cell_i == Cell_j && i == j)
   {
      return UnitCellMPO(Sites, identity_mpo(*Sites), LatticeCommute::Bosonic, Cell_i*this->size());
   }

   BasisList Basis_i = this->operator[](i).Basis1();
   BasisList Basis_j = this->operator[](j).Basis1();

   // Construct the parity operators
   // TODO: we should verify that the P operator is diagonal.  We are also
   // assuming that the diagonal matrix elements are purely real.
   LinearAlgebra::Vector<double> Parity_i(Basis_i.size(), 1.0);
   for (unsigned n = 0; n < Parity_i.size(); ++n)
   {
      Parity_i[n] = this->operator[](i)["P"](n,n).real();
   }

   LinearAlgebra::Vector<double> Parity_j(Basis_j.size(), 1.0);
   for (unsigned n = 0; n < Parity_j.size(); ++n)
   {
      Parity_j[n] = this->operator[](j)["P"](n,n).real();
   }

   ProductBasis<BasisList, BasisList> Basis_ij(Basis_i, Basis_j);
   ProductBasis<BasisList, BasisList> Basis_ji(Basis_j, Basis_i);

   // The actual gate operator
   SimpleOperator Op = ::swap_gate_fermion(Basis_i, Parity_i,
                                           Basis_j, Parity_j,
                                           Basis_ji, Basis_ij);

   // decompose it back into sites
   OperatorComponent R1, R2;
   std::tie(R1, R2) = decompose_local_tensor_prod(Op, Basis_ji, Basis_ij);

   // now turn this into a FiniteMPO
   FiniteMPO Result(this->size() * (Cell_j-Cell_i+1));

   // components up to site i are identity
   for (int n = 0; n <i; ++n)
   {
      Result[n] = OperatorComponent::make_identity(this->LocalBasis(n % this->size()));
   }

   // site i
   Result[i] = R1;

   // sites between i and j need to get the parity treatment
   for (int n = i + 1; n < (Cell_j-Cell_i)*this->size() + j; ++n)
   {
      // Construct the operator component manually
      // if Parity_i * Parity_j is +1, then the component is the identity
      // if -1 then its the local P operator
      Result[n] = OperatorComponent(this->LocalBasis(n % this->size()), R1.Basis2(), R1.Basis2());
      for (unsigned p = 0; p < Basis_ji.size(); ++p)
      {
         std::pair<int,int> x = Basis_ji.rmap(p);
         // cast to int to make sure we're comparing properly
         if (int(Parity_i[x.second]) == int(Parity_j[x.first]))
         {
            // identity
            Result[n](p,p) = this->operator[](n % this->size())["I"];
         }
         else
         {
            // parity
            Result[n](p,p) = this->operator[](n % this->size())["P"];
         }
      }
   }

   // site j
   Result[(Cell_j-Cell_i)*this->size() + j] = R2;

   // remaining sites to the end of the unit cell are identity
   for (int n = (Cell_j-Cell_i)*this->size() + j + 1; n < (Cell_j-Cell_i+1)*this->size(); ++n)
   {
      Result[n] = OperatorComponent::make_identity(this->LocalBasis(n % this->size()));
   }

   Result.debug_check_structure();
   return UnitCellMPO(Sites, Result, LatticeCommute::Bosonic, Cell_i*this->size());
}

UnitCellMPO
UnitCell::swap_gate(int i, int j) const
{
   return this->swap_gate(0, i, 0, j);
}

UnitCellMPO
UnitCell::swap_gate_no_sign(int i, int j) const
{
   return this->swap_gate_no_sign(0, i, 0, j);
}

UnitCellMPO
UnitCell::swap_gate_no_sign(int Cell_i, int i, int Cell_j, int j) const
{
   CHECK(i >= 0 && i < this->size())("Site index is outside the unit cell!")(i)(this->size());
   CHECK(j >= 0 && j < this->size())("Site index is outside the unit cell!")(j)(this->size());
   // normal-order the sites
   if (Cell_i > Cell_j || (Cell_i == Cell_j && i > j))
   {
      std::swap(Cell_i, Cell_j);
      std::swap(i,j);
   }

   if (Cell_i == Cell_j && i == j)
   {
      return UnitCellMPO(Sites, identity_mpo(*Sites), LatticeCommute::Bosonic, Cell_i*this->size());
   }

   BasisList Basis_i = this->operator[](i).Basis1();
   BasisList Basis_j = this->operator[](j).Basis1();

   ProductBasis<BasisList, BasisList> Basis_ij(Basis_i, Basis_j);
   ProductBasis<BasisList, BasisList> Basis_ji(Basis_j, Basis_i);

   // The actual gate operator
   SimpleOperator Op = ::swap_gate(Basis_i, Basis_j, Basis_ji, Basis_ij);

   // decompose it back into sites
   OperatorComponent R1, R2;
   std::tie(R1, R2) = decompose_local_tensor_prod(Op, Basis_ji, Basis_ij);

   // now turn this into a FiniteMPO
   FiniteMPO Result(this->size() * (Cell_j-Cell_i+1));
   for (int n = 0; n <i; ++n)
   {
      Result[n] = OperatorComponent::make_identity(this->LocalBasis(n % this->size()));
   }
   Result[i] = R1;
   for (int n = i + 1; n < (Cell_j-Cell_i)*this->size() + j; ++n)
   {
      Result[n] = OperatorComponent::make_identity(this->LocalBasis(n % this->size()), R1.Basis2());
   }
   Result[(Cell_j-Cell_i)*this->size() + j] = R2;
   for (int n = (Cell_j-Cell_i)*this->size() + j + 1; n < (Cell_j-Cell_i+1)*this->size(); ++n)
   {
      Result[n] = OperatorComponent::make_identity(this->LocalBasis(n % this->size()));
   }

   Result.debug_check_structure();
   return UnitCellMPO(Sites, Result, LatticeCommute::Bosonic, Cell_i*this->size());
}

void
UnitCell::SetDefaultOperators()
{
   // do nothing if the UnitCell is null
   if (this->empty())
      return;

   Operators["I"] = UnitCellMPO(Sites, identity_mpo(*Sites), LatticeCommute::Bosonic, 0,
                                "Identity");

   // we don't need the R operator if we have string() in the parser.
   //   Operators["R"] = UnitCellMPO(Sites, string_mpo(Sites, ), LatticeCommute::Bosonic, 0);
}

#if 0
FiniteMPO
UnitCell::identity_mpo(QuantumNumbers::QuantumNumber const& q) const
{
   FiniteMPO Result(this->size());
   BasisList b = make_single_basis(q);
   for (int i = 0; i < this->size(); ++i)
   {
      Result[i] = OperatorComponent(this->LocalBasis(i), b, b);
      Result[i](0,0) = this->operator[](i)["I"];
   }
   return Result;
}

FiniteMPO
UnitCell::identity_mpo() const
{
   return this->identity_mpo(QuantumNumbers::QuantumNumber(this->GetSymmetryList()));
}

FiniteMPO
UnitCell::string_mpo(std::string const& OpName, QuantumNumbers::QuantumNumber const& Trans) const
{
   FiniteMPO Result(this->size());

   BasisList Vacuum = make_single_basis(Trans);

   // Assemble the JW-string
   for (int i = 0; i < this->size(); ++i)
   {
      if (!this->operator_exists(OpName))
      {
         WARNING("JW-string operator doesn't exist at a lattice site, using the identity")(i)(OpName);
      }
      SimpleOperator Op = (*Data)[i].operator_exists(OpName) ? (*Data)[i][OpName] : (*Data)[i].identity();
      Result[i] = OperatorComponent(Op.Basis1(), Op.Basis2(), Vacuum, Vacuum);
      Result[i](0,0) = Op;
   }
   return Result;
}

FiniteMPO
UnitCell::string_mpo(LatticeCommute Com, QuantumNumbers::QuantumNumber const& Trans) const
{
   return this->string_mpo(Com.SignOperator(), Trans);
}
#endif


std::complex<double>
UnitCell::arg(std::string const& a) const
{
   const_argument_iterator I = this->find_arg(a);
   if (I != this->end_arg())
      return I->second;
   return 0.0;
}

void
UnitCell::check_structure() const
{
   if (this->empty())
      return;

   SymmetryList SL = this->operator[](0).GetSymmetryList();
   for (auto const& x : *Sites)
   {
      CHECK_EQUAL(x.GetSymmetryList(), SL)("Sites in the unit cell do not have the same SymmetryList!");
      x.check_structure();
   }
}

void
UnitCell::debug_check_structure() const
{
#if !defined(NDEBUG)
   this->check_structure();
#endif
}

PStream::opstream& operator<<(PStream::opstream& out, UnitCell const& L)
{
   out << L.Sites;
   out << L.Operators;
   out << L.Arguments;
   out << L.Functions;
   return out;
}

PStream::ipstream& operator>>(PStream::ipstream& in, UnitCell& L)
{
   in >> L.Sites;
   in >> L.Operators;
   in >> L.Arguments;
   in >> L.Functions;
   return in;
}

UnitCell repeat(UnitCell const& x, int RepeatCount)
{
   return UnitCell(RepeatCount, x);
}

UnitCell join(UnitCell const& x, UnitCell const& y)
{
   return UnitCell(x, y);
}

UnitCell join(UnitCell const& x, UnitCell const& y, UnitCell const& z)
{
   return UnitCell(x,y,z);
}

UnitCell join(UnitCell const& x, UnitCell const& y, UnitCell const& z, UnitCell const& w)
{
   return UnitCell(x,y,z,w);
}

UnitCell join(UnitCell const& x, UnitCell const& y, UnitCell const& z, UnitCell const& w,
             UnitCell const& v)
{
   return join(UnitCell(x,y,z,w), v);
}

#if 0
FiniteMPO
UnitCell::Parse(std::string const& s)
{
}
#endif
