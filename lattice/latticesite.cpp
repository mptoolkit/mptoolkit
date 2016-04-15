// -*- C++ -*- $Id$

#include "latticesite.h"
#include "siteoperator-parser.h"
#include <boost/variant/get.hpp>
#include "pstream/variant.h"

std::complex<double>
LatticeSite::arg(std::string const& a) const
{
   const_argument_iterator I = this->find_arg(a);
   if (I != this->end_arg())
      return I->second;
   return 0.0;
}

LatticeSite::operator_type const&
LatticeSite::operator[](std::string const& s) const
{ 
   OperatorListType::const_iterator I = pImpl->Operators.find(s); 
   CHECK(I != pImpl->Operators.end()) << "The site does not contain any operator named " << s;
   return I->second; 
}

SiteOperator LatticeSite::identity() const
{
   OperatorListType::const_iterator I = pImpl->Operators.find("I"); 
   CHECK(I != pImpl->Operators.end()) << "The site does not contain the identity operator!";
   return I->second;
}

SymmetryList
LatticeSite::GetSymmetryList() const
{
   CHECK(!pImpl->Operators.empty());
   // if we're debugging, verify that the symmetry list is the same for all operators
#if !defined(NDEBUG)
   SymmetryList ThisSymmetryList = this->operator[]("I").GetSymmetryList();
   for (const_operator_iterator I = pImpl->Operators.begin(); I != pImpl->Operators.end(); ++I)
   {
      CHECK_EQUAL(I->second.GetSymmetryList(), ThisSymmetryList)(I->first);
   }
#endif
   return this->operator[]("I").GetSymmetryList();
}

LatticeSite::basis1_type const&
LatticeSite::Basis1() const
{
   CHECK(!pImpl->Operators.empty());
   // if we're debugging, verify that the basis is the same for all operators
#if !defined(NDEBUG)
   basis1_type ThisBasis1 = this->operator[]("I").Basis1();
   for (const_operator_iterator I = pImpl->Operators.begin(); I != pImpl->Operators.end(); ++I)
   {
      CHECK_EQUAL(I->second.Basis1(), ThisBasis1);
   }
#endif
   return this->operator[]("I").Basis1();
}

LatticeSite::basis2_type const&
LatticeSite::Basis2() const
{
   CHECK(!pImpl->Operators.empty());
   // if we're debugging, verify that the basis is the same for all operators
#if !defined(NDEBUG)
   basis1_type ThisBasis2 = this->operator[]("I").Basis2();
   for (const_operator_iterator I = pImpl->Operators.begin(); I != pImpl->Operators.end(); ++I)
   {
      CHECK_EQUAL(I->second.Basis2(), ThisBasis2);
   }
#endif
   return this->operator[]("I").Basis2();
}

void
LatticeSite::set_operator_descriptions(OperatorDescriptions const& Desc)
{
   // iterate through the descriptions
   for (OperatorDescriptions::const_iterator I = Desc.begin(); I != Desc.end(); ++I)
   {
      if (this->operator_exists(I->first))
      {
	 this->operator[](I->first).set_description(I->second);
      }
      else
      {
	 std::cerr << "warning: operator " << I->first 
		   << " has a description but is not defined in the lattice site.\n";
      }
   }

   // Now go through the operators and check for any that don't have a description
   for (const_operator_iterator I = this->begin_operator(); I != this->end_operator(); ++I)
   {
      if (I->second.description().empty())
      {
	 std::cerr << "warning: local operator " << I->first << " has no description.\n";
      }
   }
}


// functor to parse default arguments of a LatticeSite operator
struct ParseSiteExpression
{
   ParseSiteExpression(LatticeSite const& Site_) : Site(Site_) {}

   std::complex<double> operator()(Function::ArgumentList const& Args,
				   std::string const& Str) const
   {
      return ParseSiteNumber(Site, Str, Args);
   }

   LatticeSite const& Site;
};

boost::variant<SiteOperator, std::complex<double> >
LatticeSite::eval_function(Function::OperatorFunction const& Func, 
			   Function::ParameterList const& Params) const
{
   Function::ArgumentList Args = GetArguments(Func.Args, Params, ParseSiteExpression(*this));
   return ParseSiteElement(*this, Func.Def, Args);
}

boost::variant<SiteOperator, std::complex<double> >
LatticeSite::eval_function(std::string const& Func, 
			   Function::ParameterList const& Params) const
{
   const_function_iterator I = this->find_function(Func);
   CHECK(I != this->end_function())("Function is not defined!")(Func);
   return this->eval_function(I->second, Params);
}

PStream::opstream& operator<<(PStream::opstream& out, LatticeSite const& B)
{
   return out << B.pImpl;
}

PStream::ipstream& operator>>(PStream::ipstream& in, LatticeSite& B)
{
   return in >> B.pImpl;
}


PStream::opstream& operator<<(PStream::opstream& out, LatticeSite::ImplType const& Impl)
{
   out << Impl.Description;
   out << Impl.Operators;
   out << Impl.Arguments;
   out << Impl.Functions;
   return out;
}

PStream::ipstream& operator>>(PStream::ipstream& in, LatticeSite::ImplType& Impl)
{
   in >> Impl.Description;
   in >> Impl.Operators;
   in >> Impl.Arguments;
   in >> Impl.Functions;
   return in;
}

void
CoerceSymmetryListInPlace(LatticeSite& s, QuantumNumbers::SymmetryList const& sl)
{
   using ::CoerceSymmetryListInPlace;
   LatticeSite::ptr_type::lock_type Lock(s.pImpl.lock());
   for (auto& x : Lock->Operators)
   {
      CoerceSymmetryListInPlace(x.second, sl);
   }
}

LatticeSite
CoerceSymmetryList(LatticeSite const& s, QuantumNumbers::SymmetryList const& sl)
{
   LatticeSite Result(s);
   CoerceSymmetryListInPlace(Result, sl);
   return Result;
}

std::ostream& operator<<(std::ostream& out, LatticeSite const& s)
{
   out << "LatticeSite: " << s.Description() << '\n';
   out << "Symmetry: " << s.GetSymmetryList() << '\n';
   out << "Arguments:" << (s.arg_empty() ? " none\n" : "\n");
   for (LatticeSite::const_argument_iterator I = s.begin_arg(); I != s.end_arg(); ++I)
   {
      out << I->first << " = " << I->second << '\n';
   }
   out << "Operators:\n";
   for (LatticeSite::const_operator_iterator I = s.begin_operator(); I != s.end_operator(); ++I)
   {
      out << I->first << " transforms as " << I->second.TransformsAs() << '\n';
   }
   out << "Functions:" << (s.function_empty() ? " none\n" : "\n");
   for (LatticeSite::const_function_iterator I = s.begin_function(); I != s.end_function(); ++I)
   {
      out << I->first << I->second << '\n';
   }
   return out;
}

std::vector<BasisList>
Basis1FromSiteList(SiteListType const& s)
{
   std::vector<BasisList> Result;
   for (LatticeSite const& x : s)
   {
      Result.push_back(x.Basis1());
   }
   return Result;
}

void
LatticeSite::check_structure() const
{
   basis1_type ThisBasis1 = this->operator[]("I").Basis1();
   for (auto const& x : pImpl->Operators)
   {
      CHECK_EQUAL(x.second.Basis1(), ThisBasis1);
      x.second.check_structure();
   }
}

void
LatticeSite::debug_check_structure() const
{
#if !defined(NDEBUG)
   this->check_structure();
#endif
}
