// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// lattice/siteoperator.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#if !defined(SITEOPERATOR_H_CHIUH38Y9GYEA89EY9)
#define SITEOPERATOR_H_CHIUH38Y9GYEA89EY9

#include "sitebasis.h"
#include "linearalgebra/eigen.h"
#include "tensor/tensor_exponential.h"

using Tensor::IrredTensor;

struct LatticeCommute
{
   enum Values { Fermionic = -1, None = 0, Bosonic = 1, Custom = 2 };
   
   LatticeCommute() : Value_(None) {}

   LatticeCommute(Values v) : Value_(v) {}

   LatticeCommute(std::string const& s) : Value_(Custom), SignOperator_(s) {}

   LatticeCommute& operator=(Values v)
   {
      Value_ = v;
      return *this;
   }

   std::string SignOperator() const 
   {
      if (Value_ == Fermionic)
	 return "P";
      else if (Value_ == Bosonic)
	 return "I";
      else if (Value_ == Custom)
	 return SignOperator_;
      PANIC("Undefined commutation operation");
      return "";
   }
   
   Values Value_;
   std::string SignOperator_;
};

inline
LatticeCommute operator*(LatticeCommute x, LatticeCommute y)
{
   if (x.Value_ == LatticeCommute::Custom || y.Value_ == LatticeCommute::Custom)
   {
      if (x.SignOperator() == y.SignOperator())
	 return LatticeCommute(LatticeCommute::Bosonic);
      else
	 PANIC("unsupported custom commutation string")(x.SignOperator())(y.SignOperator());
   }
      
   return LatticeCommute(LatticeCommute::Values(int(x.Value_) * (int(y.Value_))));
}

inline LatticeCommute operator*(LatticeCommute x, int N)
{
   // if x is None then return None.
   // Otherwise, if N is even then we must be bosonic.
   // Otherwise, if N is odd then return x.
   return (int(x.Value_) != 0 && N % 2 == 0) ? LatticeCommute::Bosonic : x;
}

inline LatticeCommute operator*(int N, LatticeCommute x)
{
   return x*N;
}

inline LatticeCommute& operator*=(LatticeCommute& x, LatticeCommute const& y)
{
   x = x*y;
   return x;
}

inline LatticeCommute& operator*=(LatticeCommute& x, int N)
{
   x = x*N;
   return x;
}

inline bool operator==(LatticeCommute const& x, LatticeCommute const& y)
{
   return x.Value_ == y.Value_ && (x.Value_ != LatticeCommute::Custom || x.SignOperator_ == y.SignOperator_);
}

inline bool operator!=(LatticeCommute const& x, LatticeCommute const& y)
{
   return !(x == y);
}

inline
std::ostream& operator<<(std::ostream& out, LatticeCommute const& x)
{
   if (x.Value_ == LatticeCommute::Fermionic)
      out << "Fermionic";
   else if (x.Value_ == LatticeCommute::Bosonic)
      out << "Bosonic";
   else if (x.Value_ == LatticeCommute::None)
      out << "None";
   else if (x.Value_ == LatticeCommute::Custom)
      out << "Custom(" << x.SignOperator_ << ")";
   return out;
}

inline
PStream::opstream& operator<<(PStream::opstream& out, LatticeCommute const& x)
{
   out << int(x.Value_);
   if (x.Value_ == LatticeCommute::Custom)
   {
      out << x.SignOperator_;
   }
   return out;
}

inline
PStream::ipstream& operator>>(PStream::ipstream& in, LatticeCommute& x)
{
   int Com;
   in >> Com;
   x.Value_ = LatticeCommute::Values(Com);
   if (x.Value_ == LatticeCommute::Custom)
   {
      in >> x.SignOperator_;
   }
   return in;
}

class SiteOperator : public IrredTensor<std::complex<double> >
{
   public:
      typedef IrredTensor<std::complex<double> > base_type;

      // To construct an operator that acts on a single site but is defined on
      // the whole lattice, it is necessary to know what the commutation relations are
      // for two operators acting on different sites.  In general this could be
      // arbitrary, however we distinguish two common cases: 
      //
      // Bosonic operators commute and the lattice representation of a site operator S 
      // acting on site n is
      // I_1 \otimes ... \otimes I_{n-1} \otimes S_n \times I_{n+1} \otimes ...
      // where I_j is the identity operator acting on site j.
      //
      // Fermionic operators anticommute, and have a phase factor (-1)^{\sum_{j<n} N_j}
      // for a lattice representation of
      // P_1 \otimes ... \otimes P_{n-1} \otimes S_n \times I_{n+1} \otimes ...
      // where P is the on-site parity operator = (-1)^N.

      typedef base_type::value_type value_type;

      typedef SiteBasis basis1_type;
      typedef SiteBasis basis2_type;

      using base_type::operator();

      SiteOperator() {}

      SiteOperator(SiteBasis const& B, QuantumNumber const& q, LatticeCommute Com = LatticeCommute::None,
		   std::string Description = "")
	 : base_type(B, q), Basis_(B), Com_(Com), Description_(Description) {}

      SiteOperator(SiteBasis const& B, base_type const& b, LatticeCommute Com = LatticeCommute::None,
		   std::string Description = "")
	 : base_type(b), Basis_(B), Com_(Com), Description_(Description)
	 { CHECK_EQUAL(b.Basis1(), b.Basis2()); CHECK_EQUAL(B, b.Basis1()); }

      SiteBasis const& Basis() const { return Basis_; }
      SiteBasis const& Basis1() const { return Basis_; }
      SiteBasis const& Basis2() const { return Basis_; }

      int Lookup(std::string const& Ident) const { return Basis_.Lookup(Ident); }

      value_type operator()(std::string const& si, std::string const& sj) const 
      { return this->operator()(Basis_.Lookup(si), Basis_.Lookup(sj)); }

      value_type& operator()(std::string const& si, std::string const& sj)
      { return this->operator()(Basis_.Lookup(si), Basis_.Lookup(sj)); }

      base_type& base() { return *this; }
      base_type const& base() const { return *this; }

      LatticeCommute Commute() const { return Com_; }

      void SetCommute(LatticeCommute x) { Com_ = x; }

      std::string const& description() const { return Description_; }
      std::string description_or_none() const 
      { return Description_.empty() ? "(no description)" : Description_; }
      void set_description(std::string const& s) { Description_ = s; }

      // Make an identity operator over the given basis
      static SiteOperator Identity(SiteBasis const& Basis);

      // Make an identity operator over the given basis - this is an error if Basis1 != Basis2
      static SiteOperator Identity(SiteBasis const& Basis1, SiteBasis const& Basis2);

   private:
      SiteBasis Basis_;
      LatticeCommute Com_;
      std::string Description_;

   friend PStream::opstream& operator<<(PStream::opstream& out, SiteOperator const& Op);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, SiteOperator& Op);
   friend void CoerceSymmetryListInPlace(SiteOperator& s, SymmetryList const& sl);
};

// Make an identity operator over the same basis as for operator x
SiteOperator MakeIdentityFrom(SiteOperator const& x);

std::ostream& operator<<(std::ostream& out, SiteOperator const& Op);

void CoerceSymmetryListInPlace(SiteOperator& s, SymmetryList const& sl);

inline
SiteOperator
CoerceSymmetryList(SiteOperator const& s, SymmetryList const& sl)
{
   SiteOperator Result(s);
   CoerceSymmetryListInPlace(Result, sl);
   return Result;
}

SiteOperator flip_conj(SiteOperator const& s, SiteBasis const& ReflectedBasis);

inline
SiteOperator flip_conj(SiteOperator const& s)
{
   return flip_conj(s, adjoint(s.Basis()));
}

// this should be a function in LinearAlgebra, not sure if it exists yet
// Power of a site operator, n must be positive.
SiteOperator
pow(SiteOperator const& Op, int n);

SiteOperator
dot(SiteOperator const& Op1, SiteOperator const& Op2);

// For completeness, the cross product
SiteOperator
cross(SiteOperator const& x, SiteOperator const& y);

// And the outer product
SiteOperator
outer(SiteOperator const& x, SiteOperator const& y);

namespace LinearAlgebra
{

template <>
struct interface<SiteOperator>
{
   typedef void type;
};

// unary operators

template <>
struct NormFrobSq<SiteOperator>
{
   typedef SiteOperator const& argument_type;
   typedef double result_type;

   result_type operator()(argument_type x) const 
   { 
      return norm_frob_sq(x.base());
   }
};

template <>
struct Herm<SiteOperator>
{
   typedef SiteOperator const& argument_type;
   typedef HermitianProxy<SiteOperator> result_type;

   result_type operator()(argument_type x) const 
   { return result_type(x); }
};

template <>
struct Conj<SiteOperator>
{
   typedef SiteOperator const& argument_type;
   typedef SiteOperator result_type;

   result_type operator()(argument_type x) const 
   { return SiteOperator(x.Basis(), conj(x.base()), x.Commute()); }
};

template <>
struct Adjoint<SiteOperator>
{
   typedef SiteOperator const& argument_type;
   typedef SiteOperator result_type;

   result_type operator()(argument_type x) const 
   { return SiteOperator(x.Basis(), adjoint(x.base()), x.Commute()); }
};

template <>
struct InvAdjoint<SiteOperator>
{
   typedef SiteOperator const& argument_type;
   typedef SiteOperator result_type;

   result_type operator()(argument_type x) const 
   { return SiteOperator(x.Basis(), inv_adjoint(x.base()), x.Commute()); }
};

template <>
struct Negate<SiteOperator>
{
   typedef SiteOperator const& argument_type;
   typedef SiteOperator result_type;

   result_type operator()(argument_type x) const 
   { return SiteOperator(x.Basis(), -x.base(), x.Commute()); }
};

template <>
struct Exp<SiteOperator>
{
   typedef SiteOperator const& argument_type;
   typedef SiteOperator result_type;

   result_type operator()(argument_type x) const 
   { 
      SiteOperator Result(x);
      Result.base() = Tensor::Exponentiate(x.base());
      return Result;
   }
};


// binary operators

template <>
struct Equal<SiteOperator, SiteOperator> : Equal<IrredTensor<std::complex<double> >, 
                                                 IrredTensor<std::complex<double> > > {};

template <>
struct Addition<SiteOperator, SiteOperator>
{
   typedef SiteOperator const& first_argument_type;
   typedef SiteOperator const& second_argument_type;
   typedef SiteOperator result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const 
   { 
      PRECONDITION_EQUAL(x.Commute(), y.Commute());
      return SiteOperator(x.Basis(), x.base() + y.base(), x.Commute()); 
   }
};

template <>
struct Subtraction<SiteOperator, SiteOperator>
{
   typedef SiteOperator const& first_argument_type;
   typedef SiteOperator const& second_argument_type;
   typedef SiteOperator result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const 
   { 
      PRECONDITION_EQUAL(x.Commute(), y.Commute());
      return SiteOperator(x.Basis(), x.base() - y.base(), x.Commute()); 
   }
};

template <>
struct Multiplication<SiteOperator, int>
{
   typedef SiteOperator const& first_argument_type;
   typedef int second_argument_type;
   typedef SiteOperator result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const 
   { 
      return SiteOperator(x.Basis(), x.base() * double(y), x.Commute()); 
   }
};

template <>
struct Multiplication<int, SiteOperator>
{
   typedef int first_argument_type;
   typedef SiteOperator const& second_argument_type;
   typedef SiteOperator result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const 
   { 
      return SiteOperator(y.Basis(), double(x) * y.base(), y.Commute()); 
   }
};

template <>
struct Multiplication<SiteOperator, double>
{
   typedef SiteOperator const& first_argument_type;
   typedef double second_argument_type;
   typedef SiteOperator result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const 
   { 
      return SiteOperator(x.Basis(), x.base() * y, x.Commute()); 
   }
};

template <>
struct Multiplication<double, SiteOperator>
{
   typedef double first_argument_type;
   typedef SiteOperator const& second_argument_type;
   typedef SiteOperator result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const 
   { 
      return SiteOperator(y.Basis(), x * y.base(), y.Commute()); 
   }
};

template <>
struct Multiplication<SiteOperator, std::complex<double> >
{
   typedef SiteOperator const& first_argument_type;
   typedef std::complex<double> second_argument_type;
   typedef SiteOperator result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const 
   { 
      return SiteOperator(x.Basis(), x.base() * y, x.Commute()); 
   }
};

template <>
struct Multiplication<std::complex<double>, SiteOperator>
{
   typedef std::complex<double>  first_argument_type;
   typedef SiteOperator const& second_argument_type;
   typedef SiteOperator result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const 
   { 
      return SiteOperator(y.Basis(), x * y.base(), y.Commute()); 
   }
};

template <>
struct InnerProd<SiteOperator, SiteOperator>
{
   typedef SiteOperator const& first_argument_type;
   typedef SiteOperator const& second_argument_type;
   typedef std::complex<double> result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const 
   { 
      return inner_prod(x.base(), y.base());
   }
};

template <>
struct ScalarProd<SiteOperator, HermitianProxy<SiteOperator> >
{
   typedef SiteOperator const& first_argument_type;
   typedef HermitianProxy<SiteOperator> const& second_argument_type;
   typedef SiteOperator result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const 
   { 
      return SiteOperator(x.Basis(), scalar_prod(x.base(), herm(y.base().base())), 
                          x.Commute()*y.base().Commute()); 
   }
};

template <>
struct ScalarProd<HermitianProxy<SiteOperator>, SiteOperator>
{
   typedef HermitianProxy<SiteOperator> const& first_argument_type;
   typedef SiteOperator const& second_argument_type;
   typedef SiteOperator result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const 
   { 
      return SiteOperator(y.Basis(), scalar_prod(herm(x.base().base()), y.base()), 
                          x.base().Commute()*y.Commute()); 
   }
};

template <typename Nest>
struct IrredProd<SiteOperator, SiteOperator, Nest>
{
   typedef SiteOperator const& first_argument_type;
   typedef SiteOperator const& second_argument_type;
   typedef SiteOperator result_type;

   IrredProd() {}
   IrredProd(QuantumNumbers::QuantumNumber const& Trans) : Trans_(Trans) {}

   result_type operator()(first_argument_type x, second_argument_type y) const 
   {
      return SiteOperator(x.Basis(), prod(x.base(), y.base(), Trans_), x.Commute()*y.Commute());
   }

   QuantumNumbers::QuantumNumber Trans_;
};

// forward ordinary multiplication to IrredProd.  Caution to the user if the
// result is ambiguous, as usual
template <>
struct Multiplication<SiteOperator, SiteOperator> : IrredProd<SiteOperator, SiteOperator> {};

} // namespace LinearAlgebra

SiteOperator 
tensor_prod(SiteOperator const& S1, SiteOperator const& S2,
	    SiteProductBasis const& SPBasis, QuantumNumber const& q);

SiteOperator 
tensor_prod(SiteOperator const& S1, SiteOperator const& S2, QuantumNumber const& q);

SiteOperator 
tensor_prod(SiteOperator const& S1, SiteOperator const& S2);

std::ostream& operator<<(std::ostream& out, SiteOperator const& Op);

void show_projections(std::ostream& out, SiteOperator const& Op);

std::vector<std::pair<SiteOperator, SiteOperator> >
decompose_tensor_prod(SiteOperator const& S, SiteProductBasis const& SPBasis);

#include "siteoperator.cc"

#endif
