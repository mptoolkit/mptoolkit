// -*- C++ -*- $Id$

#include "latticesite.h"
#include "siteoperator-parser.h"
#include <boost/variant/get.hpp>
#include "pstream/variant.h"

void
CoerceSymmetryList(LatticeSite& site, QuantumNumbers::SymmetryList const& sl)
{
   site.CoerceSymmetryList(sl);
}

std::complex<double>
LatticeSite::arg(std::string const& a) const
{
   const_argument_iterator I = this->find_arg(a);
   if (I != this->end_arg())
      return I->second;
   return 0.0;
}

void
LatticeSite::CoerceSymmetryList(QuantumNumbers::SymmetryList const& sl)
{
   using ::CoerceSymmetryList;
   ptr_type::lock_type Lock(pImpl.lock());
   for (operator_iterator I = Lock->Operators.begin(); I != Lock->Operators.end(); ++I)
   {
      CoerceSymmetryList(I->second, sl);
   }
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
