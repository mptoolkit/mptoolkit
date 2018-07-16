// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// lattice/siteoperator.h
//
// Copyright (C) 2014-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#if !defined(SITEOPERATOR_H_CHIUH38Y9GYEA89EY9)
#define SITEOPERATOR_H_CHIUH38Y9GYEA89EY9

#include "sitebasis.h"
#include "tensor/tensor_exponential.h"

using Tensor::IrredTensor;

struct LatticeCommute
{
   enum Values { Fermionic = -1, None = 0, Bosonic = 1, Custom = 2 };

   LatticeCommute() : Value_(None) {}

   // These are redundant, but defined here for debugging
   LatticeCommute(LatticeCommute&&) noexcept = default;
   LatticeCommute(LatticeCommute const&) = default;
   LatticeCommute& operator=(LatticeCommute&&) noexcept = default;
   LatticeCommute& operator=(LatticeCommute const&) = default;
   
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

      SiteOperator() noexcept {}

      SiteOperator(SiteBasis const& B, QuantumNumber const& q, LatticeCommute Com = LatticeCommute::None,
                   std::string Description = "")
         : base_type(B, q), Basis_(B), Com_(Com), Description_(Description) {}

      SiteOperator(SiteBasis const& B, base_type b, LatticeCommute Com = LatticeCommute::None,
                   std::string Description = "")
         : base_type(std::move(b)), Basis_(B), Com_(Com), Description_(Description)
         { CHECK_EQUAL(b.Basis1(), b.Basis2()); CHECK_EQUAL(B, b.Basis1()); }

      SiteOperator(SiteOperator&& Other) noexcept
      : Basis_(std::move(Other.Basis_)),
	Com_(std::move(Other.Com_)),
	Description_(std::move(Other.Description_))
      {}

      SiteOperator(SiteOperator const& Other)
	 : base_type(copy(static_cast<base_type const&>(Other))),
	   Com_(Other.Com_),
	   Description_(Other.Description_)
      {}

      SiteOperator& operator=(SiteOperator&&) = default;
      SiteOperator& operator=(SiteOperator const&) = default;

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
CoerceSymmetryList(SiteOperator Result, SymmetryList const& sl)
{
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

inline
void
inplace_conj(SiteOperator& x)
{
   inplace_conj(x.base());
}

inline
SiteOperator conj(SiteOperator const& x)
{
   SiteOperator(x.Basis(), conj(copy(x.base())), x.Commute());
}


inline
SiteOperator adjoint(SiteOperator const& x)
{
   SiteOperator(x.Basis(), adjoint(copy(x.base())), x.Commute());
}

inline
SiteOperator inv_adjoint(SiteOperator const& x)
{
   SiteOperator(x.Basis(), inv_adjoint(copy(x.base())), x.Commute());
}

inline
SiteOperator operator-(SiteOperator const& x)
{
   return SiteOperator(x.Basis(), -copy(x.base()), x.Commute()); 
}

inline
SiteOperator exp(SiteOperator const& x)
{
   return SiteOperator(x.Basis(), Tensor::exp(x.base()), x.Commute()); 
}

inline
SiteOperator
operator+(SiteOperator const& x, SiteOperator const& y)
{
   return SiteOperator(x.Basis(), x.base() + y.base(), x.Commute());
}

inline
SiteOperator
operator-(SiteOperator const& x, SiteOperator const& y)
{
   return SiteOperator(x.Basis(), x.base() - y.base(), x.Commute());
}

inline
SiteOperator&
operator+=(SiteOperator& x, SiteOperator const& y)
{
   CHECK_EQUAL(x.Commute(), y.Commute());
   x.base() += y.base();
   return x;
}

inline
SiteOperator&
operator-=(SiteOperator& x, SiteOperator const& y)
{
   CHECK_EQUAL(x.Commute(), y.Commute());
   x.base() -= y.base();
   return x;
}

template <typename T>
inline
std::enable_if_t<blas::is_numeric_v<T>, SiteOperator&>
operator*=(SiteOperator& x, T a)
{
   x.base() += a;
}

template <typename T>
inline
std::enable_if_t<blas::is_numeric_v<T>, SiteOperator>
operator*(SiteOperator const& x, T a)
{
   return SiteOperator(x.Basis(), x.base() * a, x.Commute()); 
}

template <typename T>
inline
std::enable_if_t<blas::is_numeric_v<T>, SiteOperator>
operator*(T a, SiteOperator const& x)
{
   return SiteOperator(x.Basis(), a * x.base(), x.Commute()); 
}

inline
complex
inner_prod(SiteOperator const& x, SiteOperator const& y)
{
   return inner_prod(x.base(), y.base());
}

inline
SiteOperator
prod(SiteOperator const& x, SiteOperator const& y, QuantumNumber const& q)
{
   return SiteOperator(x.Basis(), prod(x.base(), y.base(), q), x.Commute()*y.Commute());
}

inline
SiteOperator
prod(SiteOperator const& x, SiteOperator const& y)
{
   return SiteOperator(x.Basis(), prod(x.base(), y.base()), x.Commute()*y.Commute());
}

inline
SiteOperator
operator*(SiteOperator const& x, SiteOperator const& y)
{
   return prod(x,y);
}


inline
SiteOperator
scalar_prod(SiteOperator const& x, Tensor::HermitianProxy<SiteOperator> const& y)
{
   return SiteOperator(x.Basis(), scalar_prod(x.base(), herm(y.base().base())),
		       x.Commute()*y.base().Commute());
}
   
inline
SiteOperator
scalar_prod(Tensor::HermitianProxy<SiteOperator> const& x, SiteOperator const& y)
{
   return SiteOperator(y.Basis(), scalar_prod(herm(x.base().base()), y.base()),
		       x.base().Commute()*y.Commute());
}

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
