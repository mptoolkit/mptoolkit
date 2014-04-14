// -*- C++ -*- $Id$

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

      SiteOperator(SiteBasis const& B, QuantumNumber const& q, LatticeCommute Com = LatticeCommute::None)
	 : base_type(B, q), Basis_(B), Com_(Com) {}

      SiteOperator(SiteBasis const& B, base_type const& b, LatticeCommute Com = LatticeCommute::None) 
	 : base_type(b), Basis_(B), Com_(Com)
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

      std::string const& Description() const { return Description_; }
      void SetDescription(std::string const& s) { Description_ = s; }

     static SiteOperator Identity(SiteBasis const& Basis);

   private:
      std::string Description_;
      SiteBasis Basis_;
      LatticeCommute Com_;

   friend PStream::opstream& operator<<(PStream::opstream& out, SiteOperator const& Op);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, SiteOperator& Op);
   friend void CoerceSymmetryList(SiteOperator& s, SymmetryList const& sl);
};

inline
SiteOperator SiteOperator::Identity(SiteBasis const& Basis)
{
   SiteOperator Result(Basis, QuantumNumber(Basis.GetSymmetryList()), LatticeCommute::Bosonic);
   for (std::size_t i = 0; i < Basis.size(); ++i)
   {
      Result(i,i) = 1.0;
   }
   return Result;
}

std::ostream& operator<<(std::ostream& out, SiteOperator const& Op);

inline
void CoerceSymmetryList(SiteOperator& s, SymmetryList const& sl)
{
   CoerceSymmetryList(s.Basis_, sl);
   CoerceSymmetryList(static_cast<IrredTensor<std::complex<double> >&>(s), sl);
}

SiteOperator flip_conj(SiteOperator const& s, SiteBasis const& ReflectedBasis);

inline
SiteOperator flip_conj(SiteOperator const& s)
{
   return flip_conj(s, adjoint(s.Basis()));
}

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
