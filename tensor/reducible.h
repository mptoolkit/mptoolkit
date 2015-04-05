/* -*- C++ -*- $Id$

  Reducible tensor.  This is represented as a sum of irriducible tensors,
  indexed by the quantum number.  Conceptually, it is similar to an
  MPS or a 3-index tensor, but the third index is constrained to be
  a simple list of distinct quantum numbers.

*/

#if !defined(REDUCIBLE_H_SDHCEWUYUIGHUIHYO89337489FD)
#define REDUCIBLE_H_SDHCEWUYUIGHUIHYO89337489FD
  
#include "tensor.h"
#include "tensorproduct.h"
#include <set>
#include "common/set_operations.h"
#include <boost/iterator/transform_iterator.hpp>

namespace Tensor
{

template <typename T>
struct TransformsAsCompare
{
   typedef T first_argument_type;
   typedef T second_argument_type;
   typedef bool result_type;

   result_type operator()(T const& x, T const& y) const
   {
      return x.TransformsAs() < y.TransformsAs();
   }
};

template <
    typename T, 
    typename Basis1T = BasisList, 
    typename Basis2T = Basis1T, 
    typename Structure = LinearAlgebra::SparseMatrix<T> >
class ReducibleTensor;

template <typename T, typename B1, typename B2, typename S>
PStream::opstream& operator<<(PStream::opstream& out, ReducibleTensor<T, B1, B2, S> const& x);

template <typename T, typename B1, typename B2, typename S>
PStream::ipstream& operator>>(PStream::ipstream& in, ReducibleTensor<T, B1, B2, S>& x);

template <typename T, typename B1, typename B2, typename S>
void
CoerceSymmetryList(ReducibleTensor<T, B1, B2, S>& t, SymmetryList const& sl);

template <typename T, typename Basis1T, typename Basis2T, typename Structure>
class ReducibleTensor
{
   public:
      typedef IrredTensor<T, Basis1T, Basis2T, Structure> IrredTensorType;
      typedef IrredTensorType value_type;

      typedef Basis1T basis1_type;
      typedef Basis2T basis2_type;

   private:
      typedef TransformsAsCompare<IrredTensorType> Compare;
      typedef std::map<QuantumNumber, IrredTensorType> data_type;

      struct GetSecond
      {
         typedef IrredTensorType& result_type;
         result_type
         operator()(std::pair<QuantumNumber const, IrredTensorType>& x) const
         {
            return x.second;
         }
      };
      struct GetSecondConst
      {
         typedef IrredTensorType const& result_type;
         result_type
         operator()(std::pair<QuantumNumber const, IrredTensorType> const& x) const
         {
            return x.second;
         }
      };

   public:
      typedef typename data_type::iterator map_iterator;
      typedef typename data_type::const_iterator const_map_iterator;
      typedef boost::transform_iterator<GetSecond, map_iterator> iterator; 
      typedef boost::transform_iterator<GetSecondConst, const_map_iterator> const_iterator; 

      ReducibleTensor();

      ReducibleTensor(basis1_type const& Basis);
      ReducibleTensor(basis1_type const& Basis1, basis2_type const& Basis2);

      template <typename U, typename US>
      ReducibleTensor(IrredTensor<U, basis1_type, basis2_type, US> const& x);

      template <typename U, typename US>
      ReducibleTensor(ReducibleTensor<U, basis1_type, basis2_type, US> const& x);

      bool empty() const { return data_.empty(); }

      bool is_null() const { return this->empty(); }

      // returns the number of irreducible components
      size_type size() const { return data_.size(); }

      SymmetryList GetSymmetryList() const { return Basis1_.GetSymmetryList(); }

      size_type size1() const { return Basis1_.size(); }
      size_type size2() const { return Basis2_.size(); }

      ReducibleTensor& operator+=(ReducibleTensor const& x);
      ReducibleTensor& operator-=(ReducibleTensor const& x);

      ReducibleTensor& operator+=(IrredTensorType const& x);
      ReducibleTensor& operator-=(IrredTensorType const& x);

      basis1_type const& Basis1() const { return Basis1_; }
      basis2_type const& Basis2() const { return Basis2_; }

      template <typename U>
      ReducibleTensor& operator*=(U const& x);

      // iterate over the irred tensors that comprise this object
      iterator begin() { return iterator(data_.begin(), GetSecond()); }
      iterator end() { return iterator(data_.end(), GetSecond()); }

      const_iterator begin() const { return const_iterator(data_.begin(), GetSecondConst()); }
      const_iterator end() const { return const_iterator(data_.end(), GetSecondConst()); }

      // returns an IrredTensor for the component that transforms as q
      value_type project(QuantumNumber const& q) const;
      value_type& project(QuantumNumber const& q);

      // projects onto the scalar component
      value_type scalar() const;
      value_type& scalar();

      // variant where we know in advance that the corresponding component
      // exists in the container
      value_type const& project_assert(QuantumNumber const& q) const;

      QuantumNumber const& qn1(size_type i) const { return Basis1_[i]; }
      QuantumNumber const& qn2(size_type j) const { return Basis2_[j]; }

      // return the quantum numbers for which this tensor has non-zero matrix elements
      std::set<QuantumNumber> components() const;

      void check_structure() const;

      void debug_check_structure() const;

      static ReducibleTensor<T, Basis1T, Basis2T, Structure> 
      make_identity(basis1_type const& b);

   private:
      basis1_type Basis1_;
      basis2_type Basis2_;
      data_type data_;

   template <typename U, typename B1, typename B2, typename S>
   friend class ReducibleTensor;

   friend PStream::opstream& operator<< <>(PStream::opstream& out, ReducibleTensor const& x);
   friend PStream::ipstream& operator>> <>(PStream::ipstream& in, ReducibleTensor& x);

   friend void CoerceSymmetryList<>(ReducibleTensor& t, SymmetryList const& sl);
};

// project an irred tensor onto some irreducible subspace.
template <typename T, typename B1, typename B2, typename S>
inline
IrredTensor<T, B1, B2, S>
project(ReducibleTensor<T, B1, B2, S> const& x, QuantumNumber const& q)
{
   return x.project(q);
}

// project an irred tensor onto some irreducible subspace, as an lvalue.
// If x has no component that transforms in this way, then one is created
// and initialized to zero.
template <typename T, typename B1, typename B2, typename S>
inline
IrredTensor<T, B1, B2, S>&
project(ReducibleTensor<T, B1, B2, S>& x, QuantumNumber const& q)
{
   return x.project(q);
}

template <typename T, typename B1, typename B2, typename S>
inline
IrredTensor<T, B1, B2, S>
scalar(ReducibleTensor<T, B1, B2, S> const& x)
{
   return x.scalar();
}

template <typename T, typename B1, typename B2, typename S>
inline
IrredTensor<T, B1, B2, S>&
scalar(ReducibleTensor<T, B1, B2, S>& x)
{
   return x.scalar();
}

// is_pure_scalar, returns true only if the operator contains only scalar (or zero) components
template <typename T, typename B1, typename B2, typename S>
inline
bool
is_pure_scalar(ReducibleTensor<T, B1, B2, S> const& x)
{
   std::set<QuantumNumber> comp = x.components();
   return comp.empty() || (comp.size() == 1 && comp.count(QuantumNumbers::QuantumNumber(x.GetSymmetryList())));
}

// text I/O

template <typename T, typename B1, typename B2, typename S>
std::ostream& operator<<(std::ostream& out, ReducibleTensor<T, B1, B2, S> const& Op);

template <typename T, typename B1, typename B2, typename S>
std::string show_projections(ReducibleTensor<T, B1, B2, S> const& Op);

// construction

template <typename T, typename B1, typename B2, typename S>
inline
ReducibleTensor<T, B1, B2, S>
ReducibleTensor<T, B1, B2, S>::make_identity(B1 const& b)
{
   return IrredTensor<T, B1, B2, S>::make_identity(b);
}

} // namespace Tensor

namespace LinearAlgebra
{

template <typename T, typename B1, typename B2, typename S>
struct interface<Tensor::ReducibleTensor<T, B1, B2, S> >
{
   typedef void type;
};


// zero

template <typename T, typename B1, typename B2, typename S>
struct Zero<Tensor::ReducibleTensor<T, B1, B2, S> >
{
   typedef Tensor::ReducibleTensor<T, B1, B2, S> result_type;
   result_type operator()() const { return result_type(); }
};

template <typename T, typename B1, typename B2, typename S>
struct IsZero<Tensor::ReducibleTensor<T, B1, B2, S> >
{
   typedef bool result_type;
   typedef Tensor::ReducibleTensor<T, B1, B2, S> const& argument_type;
   result_type operator()(argument_type x) const
   {
      return x.is_null();
   }
};

template <typename T, typename B1, typename B2, typename S>
struct ZeroAll<Tensor::ReducibleTensor<T, B1, B2, S>&>
{
   typedef void result_type;
   typedef Tensor::ReducibleTensor<T, B1, B2, S>& argument_type;
   void operator()(argument_type x) const
   {
      x = Tensor::ReducibleTensor<T, B1, B2, S>(x.Basis1(), x.Basis2());
   }
};

// comparison

template <typename T, typename B1, typename B2, typename S1, typename S2>
struct Equal<Tensor::ReducibleTensor<T, B1, B2, S1>,  Tensor::ReducibleTensor<T, B1, B2, S2> >
{
   typedef Tensor::ReducibleTensor<T, B1, B2, S1> const& first_argument_type;
   typedef Tensor::ReducibleTensor<T, B1, B2, S2> const& second_argument_type;
   typedef bool result_type;

   Equal(double tol = LinearAlgebra::default_tolerance()) : eq_(tol) {}

   bool operator()(first_argument_type x, second_argument_type y) const
   {
      if (x.is_null() && y.is_null()) return true;

      std::set<QuantumNumbers::QuantumNumber> q1 = x.components(), q2 = y.components();
      q1.insert(q2.begin(), q2.end());

      std::set<QuantumNumbers::QuantumNumber>::const_iterator IEnd = q1.end();
      for (std::set<QuantumNumbers::QuantumNumber>::const_iterator I = q1.begin(); 
           I != IEnd; ++I)
      {
         if (!eq_(project(x, *I), project(y, *I)))
            return false;
      }

      return true;
   }

   Equal<Tensor::IrredTensor<T, B1, B2, S1>,  Tensor::IrredTensor<T, B1, B2, S2> > eq_;
};

// arithmetic

template <typename T, typename B1, typename B2, typename S>
struct Negate<Tensor::ReducibleTensor<T, B1, B2, S> >
{
   typedef Tensor::ReducibleTensor<T, B1, B2, S> const& argument_type;
   typedef Tensor::ReducibleTensor<T, B1, B2, S> result_type;

   result_type operator()(Tensor::ReducibleTensor<T, B1, B2, S> const& x) const
   {
      result_type Result(x);
      for (typename result_type::iterator I = Result.begin(); I != Result.end(); ++I)
      {
         (*I) = -(*I);
      }
      return Result;
   }
};

template <typename T1, typename T2, typename B1, typename B2, typename S1, typename S2>
struct Addition<Tensor::ReducibleTensor<T1, B1, B2, S1>, 
                Tensor::ReducibleTensor<T2, B1, B2, S2> >
{
   typedef Tensor::ReducibleTensor<T1, B1, B2, S1> const& first_argument_type;
   typedef Tensor::ReducibleTensor<T2, B1, B2, S2> const& second_argument_type;
   typedef typename result_value<Addition<S1, S2> >::type ResultS;
   typedef Tensor::ReducibleTensor<typename value_type<ResultS>::type, 
                               B1, B2, ResultS> result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      if (x.is_null()) return y;
      if (y.is_null()) return x;
      PRECONDITION_EQUAL(x.Basis1(), y.Basis1());
      PRECONDITION_EQUAL(x.Basis2(), y.Basis2());
      result_type Result(x.Basis1(), x.Basis2());
      std::set<QuantumNumbers::QuantumNumber> q1 = x.components(), q2 = y.components();
      q1.insert(q2.begin(), q2.end());

      std::set<QuantumNumbers::QuantumNumber>::const_iterator IEnd = q1.end();
      for (std::set<QuantumNumbers::QuantumNumber>::const_iterator I = q1.begin(); 
           I != IEnd; ++I)
      {
         project(Result, *I) = project(x, *I) + project(y, *I);
      }
      
      return Result;
   }
};

template <typename T1, typename T2, typename B1, typename B2, typename S1, typename S2>
struct Subtraction<Tensor::ReducibleTensor<T1, B1, B2, S1>, 
                   Tensor::ReducibleTensor<T2, B1, B2, S2> >
{
   typedef Tensor::ReducibleTensor<T1, B1, B2, S1> const& first_argument_type;
   typedef Tensor::ReducibleTensor<T2, B1, B2, S2> const& second_argument_type;
   typedef typename result_value<Addition<S1, S2> >::type ResultS;
   typedef Tensor::ReducibleTensor<typename value_type<ResultS>::type, 
                               B1, B2, ResultS> result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      if (x.is_null()) return y;
      if (y.is_null()) return x;
      PRECONDITION_EQUAL(x.Basis1(), y.Basis1());
      PRECONDITION_EQUAL(x.Basis2(), y.Basis2());
      result_type Result(x.Basis1(), x.Basis2());
      std::set<QuantumNumbers::QuantumNumber> q1 = x.components(), q2 = y.components();
      q1.insert(q2.begin(), q2.end());

      std::set<QuantumNumbers::QuantumNumber>::const_iterator IEnd = q1.end();
      for (std::set<QuantumNumbers::QuantumNumber>::const_iterator I = q1.begin(); 
           I != IEnd; ++I)
      {
         project(Result, *I) = project(x, *I) - project(y, *I);
      }
      
      return Result;
   }
};

// multiply by scalar

template <typename T, typename B1, typename B2, typename S, typename U>
struct Multiplication<U, Tensor::ReducibleTensor<T, B1, B2, S> >
{
   typedef U first_argument_type;
   typedef Tensor::ReducibleTensor<T, B1, B2, S> const& second_argument_type;
   typedef typename result_value<Multiplication<U, S> >::type ResultS;
   typedef Tensor::ReducibleTensor<typename value_type<ResultS>::type, 
                               B1, B2, ResultS> result_type;

   result_type operator()(first_argument_type n, second_argument_type x) const
   {
      result_type Result(x.Basis1(), x.Basis2());
      for (typename result_type::const_iterator I = x.begin(); I != x.end(); ++I)
      {
         project(Result, I->TransformsAs()) = n * (*I);
      }
      return Result;
   }
};

template <typename T, typename B1, typename B2, typename S, typename U>
struct Multiplication<Tensor::ReducibleTensor<T, B1, B2, S>, U>
{
   typedef Tensor::ReducibleTensor<T, B1, B2, S> const& first_argument_type; 
  typedef U second_argument_type;
   typedef typename result_value<Multiplication<S, U> >::type ResultS;
   typedef Tensor::ReducibleTensor<typename value_type<ResultS>::type, 
                               B1, B2, ResultS> result_type;

   result_type operator()(first_argument_type x, second_argument_type n) const
   {
      result_type Result(x.Basis1(), x.Basis2());
      for (typename result_type::const_iterator I = x.begin(); I != x.end(); ++I)
      {
         project(Result, I->TransformsAs()) = (*I) * n;
      }
      return Result;
   }
};

template <typename T, typename B1, typename B2, typename S, typename Func>
struct Transform<Tensor::ReducibleTensor<T, B1, B2, S>, Func>
{
   typedef Tensor::ReducibleTensor<T, B1, B2, S> const& first_argument_type;
   typedef Func second_argument_type;
   typedef typename result_value<Transform<S, Func> >::type ResultS;
   typedef Tensor::ReducibleTensor<typename value_type<ResultS>::type, 
                                   B1, B2, ResultS> result_type;

   result_type operator()(first_argument_type x, second_argument_type f) const
   {
      result_type Result(x.Basis1(), x.Basis2());
      for (typename result_type::const_iterator I = x.begin(); I != x.end(); ++I)
      {
         project(Result, I->TransformsAs()) = transform(*I, f);
      }
      return Result;
   }
};


// complex conjugation

template <typename T, typename B1, typename B2, typename S>
struct Conj<Tensor::ReducibleTensor<T, B1, B2, S> >
{
   typedef Tensor::ReducibleTensor<T, B1, B2, S> const& argument_type;
   typedef Tensor::ReducibleTensor<T, B1, B2, S> result_type;
   
   result_type operator()(argument_type x) const
   {
      result_type Result(x);
      for (typename Tensor::ReducibleTensor<T, B1, B2, S>::iterator I = Result.begin(); 
           I != Result.end(); ++I)
      {
         (*I) = conj(*I);
      }
      return Result;
   }
};

// hermitian conjugation
  
template <typename T, typename B, typename B2, typename S>
inline
HermitianProxy<Tensor::ReducibleTensor<T, B, B2, S> >
herm(Tensor::ReducibleTensor<T, B, B2, S> const& x)
{
   return HermitianProxy<Tensor::ReducibleTensor<T, B, B2, S> >(x);
}

// flip conjugation

template <typename T, typename B1, typename B2, typename S>
struct FlipConjugate<Tensor::ReducibleTensor<T, B1, B2, S> > : TensorFlipConjugate<Tensor::ReducibleTensor<T, B1, B2, S> > {};

template <typename T, typename B1, typename B2, typename S, typename F>
struct TensorFlipConjugate<Tensor::ReducibleTensor<T, B1, B2, S>, F>
{
   TensorFlipConjugate(F f = F()) : f_(f) {}

   typedef Tensor::ReducibleTensor<typename result_value<F>::type::value_type, B2, B1> result_type;
   typedef Tensor::ReducibleTensor<T, B1, B2, S> const& argument_type;

   result_type operator()(Tensor::ReducibleTensor<T, B1, B2, S> const& x) const
   {
      result_type Result(adjoint(x.Basis1()), adjoint(x.Basis2()));
      for (typename Tensor::ReducibleTensor<T, B1, B2, S>::const_iterator I = x.begin(); I != x.end(); ++I)
      {
         Result.project(adjoint(I->TransformsAs())) = flip_conj(*I);
      }
      return Result;
   }
   F f_;
};

// trace

template <typename T, typename B1, typename B2, typename S>
struct Trace<Tensor::ReducibleTensor<T, B1, B2, S> >
{
   typedef Tensor::ReducibleTensor<T, B1, B2, S> const& argument_type;
   typedef typename Trace<Tensor::IrredTensor<T, B1, B2, S> >::result_type result_type;

   result_type operator()(argument_type x) const
   {
      return trace(project(x, QuantumNumber(x.GetSymmetryList())));
   }
};

// norm_frob

template <typename T, typename B1, typename B2, typename S>
struct NormFrobSq<Tensor::ReducibleTensor<T, B1, B2, S> >
{
   typedef Tensor::ReducibleTensor<T, B1, B2, S> const& argument_type;
   typedef typename result_value<NormFrobSq<T> >::type result_type;

   result_type operator()(argument_type x) const
   {
      result_type Result = 0;
      for (typename Tensor::ReducibleTensor<T, B1, B2, S>::const_iterator I = x.begin(); 
           I != x.end(); ++I)
      {
         Result += norm_frob_sq(*I);
      }
      return Result;
   }
};

// inner_prod

template <typename T1, typename B1, typename B2, typename S1, typename T2, typename S2>
struct InnerProd<Tensor::ReducibleTensor<T1, B1, B2, S1>, 
                 Tensor::ReducibleTensor<T2, B1, B2, S2> >
{
   typedef Tensor::ReducibleTensor<T1, B1, B2, S1> first_argument_type;
   typedef Tensor::ReducibleTensor<T2, B1, B2, S2> second_argument_type;
   typedef typename result_value<InnerProd<T1, T2> >::type result_type;

   result_type operator()(first_argument_type const& x, second_argument_type const& y) const
   {
      using QuantumNumbers::QuantumNumber;
      result_type Result = 0;
      std::set<QuantumNumber> q1 = x.components(), q2 = y.components();
      intersection_iterator<std::set<QuantumNumber> > IEnd = set_intersection_end(q1,q2);
      for (intersection_iterator<std::set<QuantumNumber> > I = set_intersection_begin(q1,q2);
	   I != IEnd; ++I)
      {
	 Result += inner_prod(x.project_assert(*I), y.project_assert(*I));
      }
      return Result;
   }
};

// parallel_prod - do we need this?

// scalar_prod

template <typename T1, typename B1, typename B2, typename S1, 
          typename T2, typename B3, typename S2, typename Functor>
struct ScalarProd<Tensor::ReducibleTensor<T1, B1, B2, S1>, 
		  HermitianProxy<Tensor::ReducibleTensor<T2, B3, B2, S2> >, 
                  Functor>
{
   typedef Tensor::ReducibleTensor<T1, B1, B2, S1> const& first_argument_type;
   typedef HermitianProxy<Tensor::ReducibleTensor<T2, B3, B2, S2> > const& second_argument_type;
   typedef Tensor::IrredTensor<T1, B1, B3> result_type;

   ScalarProd(Functor f = Functor()) : Impl_(f) {}

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      using QuantumNumbers::QuantumNumber;
      std::set<QuantumNumber> q1 = x.components(), q2 = y.base().components();
      intersection_iterator<std::set<QuantumNumber> > IEnd = set_intersection_end(q1,q2);
      intersection_iterator<std::set<QuantumNumber> > I = set_intersection_begin(q1,q2);
      if (I == IEnd)
	 return result_type();

      result_type Result = scalar_prod(x.project_assert(*I), herm(y.base().project_assert(*I)));
      //Impl_(x.project_assert(*I), herm(y.base().project_assert(*I)));
      ++I;
      while (I != IEnd)
      {
	 Result += scalar_prod(x.project_assert(*I), herm(y.base().project_assert(*I)));
         //Impl_(x.project_assert(*I), herm(y.base().project_assert(*I)));
	 ++I;
      }
      return Result;
   }

   ScalarProd<Tensor::IrredTensor<T1, B1, B2, S1>, 
	      HermitianProxy<Tensor::IrredTensor<T2, B3, B2, S2> >, 
	      Functor> Impl_;
};
      
template <typename T1, typename B1, typename B2, typename S1, 
          typename T2, typename B3, typename S2, typename Functor>
struct ScalarProd<HermitianProxy<Tensor::ReducibleTensor<T1, B1, B2, S1> >, 
                  Tensor::ReducibleTensor<T2, B1, B3, S2>, 
                  Functor>
{
   typedef HermitianProxy<Tensor::ReducibleTensor<T1, B1, B2, S1> > const& first_argument_type;
   typedef Tensor::ReducibleTensor<T2, B1, B3, S2> const& second_argument_type;
   typedef Tensor::IrredTensor<T1, B2, B3> result_type;

   ScalarProd(Functor f = Functor()) : Impl_(f) {}

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      using QuantumNumbers::QuantumNumber;
      std::set<QuantumNumber> q1 = x.base().components(), q2 = y.components();
      intersection_iterator<std::set<QuantumNumber> > IEnd = set_intersection_end(q1,q2);
      intersection_iterator<std::set<QuantumNumber> > I = set_intersection_begin(q1,q2);
      if (I == IEnd)
	 return result_type();

      result_type Result = scalar_prod(herm(x.base().project_assert(*I)), y.project_assert(*I));
      ++I;
      while (I != IEnd)
      {
	 Result += scalar_prod(herm(x.base().project_assert(*I)), y.project_assert(*I));
	 ++I;
      }
      return Result;
   }
   
   ScalarProd<HermitianProxy<Tensor::IrredTensor<T1, B1, B2, S1> >, 
	      Tensor::IrredTensor<T2, B1, B3, S2>, 
	      Functor> Impl_;
};

//
// this needs to use IrredProd directly.
//

template <typename S, typename T, typename Nest>
struct ReducibleProd {};


template <typename S, typename T, typename Nest>
inline
typename ReducibleProd<S, T, Nest>::result_type
prod(S const& x, T const& y, Nest const& n)
{
   return ReducibleProd<S, T, Nest>()(x,y,n);
}

template <typename T1, typename B1, typename B2, typename S1, 
          typename T2, typename B3, typename S2, typename Nest>
struct ReducibleProd<Tensor::ReducibleTensor<T1, B1, B2, S1>, 
                     Tensor::ReducibleTensor<T2, B2, B3, S2>, Nest >
{
   typedef Tensor::ReducibleTensor<T1, B1, B2, S1> const& first_argument_type;
   typedef Tensor::ReducibleTensor<T2, B2, B3, S2> const& second_argument_type;
   typedef Nest third_argument_type;
   typedef Tensor::ReducibleTensor<typename result_value<Nest>::type, B1, B3> result_type;
   typedef IrredProd<Tensor::IrredTensor<T1, B1, B2, S1>, Tensor::IrredTensor<T2,B2,B3,S2> > IrredProdType;

   result_type operator()(first_argument_type x, second_argument_type y, third_argument_type f) const
   {
      using QuantumNumbers::QuantumNumber;
      using QuantumNumbers::QuantumNumberList;
      typedef Tensor::ReducibleTensor<T1, B1, B2, S1> t1;
      typedef Tensor::ReducibleTensor<T2, B2, B3, S2> t2;
      result_type Result(x.Basis1(), y.Basis2());
      for (typename t1::const_iterator I1 = x.begin(); I1 != x.end(); ++I1)
      {
	 for (typename t2::const_iterator I2 = y.begin(); I2 != y.end(); ++I2)
	 {
	    QuantumNumberList QL = transform_targets(I1->TransformsAs(),I2->TransformsAs());
	    for (QuantumNumberList::const_iterator I = QL.begin(); I != QL.end(); ++I)
	    {
	       project(Result, *I) += IrredProdType(*I)(*I1, *I2, f);
	    }
	 }
      }
      return Result;
   }

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      return this->operator()(x,y, Nest());
   }

};

template <typename T1, typename B1, typename B2, typename S1, 
          typename T2, typename B3, typename S2, typename Nest>
struct ReducibleProd<Tensor::ReducibleTensor<T1, B1, B2, S1>, 
                     Tensor::IrredTensor<T2, B2, B3, S2>, Nest >
{
   typedef Tensor::ReducibleTensor<T1, B1, B2, S1> const& first_argument_type;
   typedef Tensor::IrredTensor<T2, B2, B3, S2> const& second_argument_type;
   typedef Nest third_argument_type;
   typedef Tensor::ReducibleTensor<typename result_value<Nest>::type, B1, B3> result_type;
   typedef IrredProd<Tensor::IrredTensor<T1, B1, B2, S1>, Tensor::IrredTensor<T2,B2,B3,S2> > IrredProdType;

   result_type operator()(first_argument_type x, second_argument_type y,
                          Nest F) const
   {
      using QuantumNumbers::QuantumNumber;
      using QuantumNumbers::QuantumNumberList;
      typedef Tensor::ReducibleTensor<T1, B1, B2, S1> t1;
      typedef Tensor::IrredTensor<T2, B2, B3, S2> t2;
      result_type Result(x.Basis1(), y.Basis2());
      for (typename t1::const_iterator I1 = x.begin(); I1 != x.end(); ++I1)
      {
         QuantumNumberList QL = transform_targets(I1->TransformsAs(), y.TransformsAs());
         for (QuantumNumberList::const_iterator I = QL.begin(); I != QL.end(); ++I)
	 {
            project(Result, *I) += IrredProdType(*I)(*I1, y, *I, F);
         }
      }
      return Result;
   }
};

template <typename T1, typename B1, typename B2, typename S1, 
          typename T2, typename B3, typename S2, typename Nest>
struct ReducibleProd<Tensor::IrredTensor<T1, B1, B2, S1>, 
		     Tensor::ReducibleTensor<T2, B2, B3, S2>, Nest >
{
   typedef Tensor::IrredTensor<T1, B1, B2, S1> const& first_argument_type;
   typedef Tensor::ReducibleTensor<T2, B2, B3, S2> const& second_argument_type;
   typedef Nest third_argument_type;
   typedef Tensor::ReducibleTensor<typename result_value<Nest>::type, B1, B3> result_type;
   typedef IrredProd<Tensor::IrredTensor<T1, B1, B2, S1>, Tensor::IrredTensor<T2,B2,B3,S2> > IrredProdType;

   result_type operator()(first_argument_type x, second_argument_type y,
                          Nest F) const
   {
      using QuantumNumbers::QuantumNumber;
      using QuantumNumbers::QuantumNumberList;
      typedef Tensor::IrredTensor<T1, B1, B2, S1> t1;
      typedef Tensor::ReducibleTensor<T2, B2, B3, S2> t2;
      result_type Result(x.Basis1(), y.Basis2());
      for (typename t2::const_iterator I2 = y.begin(); I2 != y.end(); ++I2)
      {
         QuantumNumberList QL = transform_targets(x.TransformsAs(), I2->TransformsAs());
         for (QuantumNumberList::const_iterator I = QL.begin(); I != QL.end(); ++I)
	 {
            project(Result, *I) += IrredProdType(*I)(x, *I2, *I, F);
         }
      }
      return Result;
   }
};

template <typename T1, typename B1, typename B2, typename S1, 
          typename T2, typename B3, typename S2>
struct Multiplication<Tensor::ReducibleTensor<T1, B1, B2, S1>, Tensor::ReducibleTensor<T2, B2, B3, S2> >
   : ReducibleProd<Tensor::ReducibleTensor<T1, B1, B2, S1>, 
                   Tensor::ReducibleTensor<T2, B2, B3, S2>, Multiplication<T1, T2> > {};

template <typename T1, typename B1, typename B2, typename S1, 
          typename T2, typename B3, typename S2>
struct Multiplication<Tensor::ReducibleTensor<T1, B1, B2, S1>, Tensor::IrredTensor<T2, B2, B3, S2> >
   : ReducibleProd<Tensor::ReducibleTensor<T1, B1, B2, S1>, 
                   Tensor::IrredTensor<T2, B2, B3, S2>,
                   Multiplication<T1, T2> > {};

template <typename T1, typename B1, typename B2, typename S1, 
          typename T2, typename B3, typename S2>
struct Multiplication<Tensor::IrredTensor<T1, B1, B2, S1>, Tensor::ReducibleTensor<T2, B2, B3, S2> >
   : ReducibleProd<Tensor::IrredTensor<T1, B1, B2, S1>, 
                   Tensor::ReducibleTensor<T2, B2, B3, S2>,
                   Multiplication<T1, T2> > {};

template <typename T, typename B1, typename B2, typename S, typename F>
struct TensorAdjoint<Tensor::ReducibleTensor<T, B1, B2, S>, F>
{
   TensorAdjoint(F f = F()) : f_(f) {}

   typedef Tensor::ReducibleTensor<typename result_value<F>::type, B2, B1> result_type;
   typedef Tensor::ReducibleTensor<T, B1, B2, S> const& argument_type;

   result_type operator()(Tensor::ReducibleTensor<T, B1, B2, S> const& x) const
   {
      result_type Result(x.Basis2(), x.Basis1());
      for (typename Tensor::ReducibleTensor<T, B1, B2, S>::const_iterator I = x.begin(); I != x.end(); ++I)
      {
         Result += adjoint(*I, f_);
      }
      return Result;
   }

   F f_;
};

template <typename T, typename B1, typename B2, typename S>
struct Adjoint<Tensor::ReducibleTensor<T, B1, B2, S> > : TensorAdjoint<Tensor::ReducibleTensor<T, B1, B2, S>,
                                                                       Adjoint<typename Tensor::IrredTensor<T, B1, B2, S>::value_type> > {};

template <typename T, typename B1, typename B2, typename S, typename F>
struct TensorInvAdjoint<Tensor::ReducibleTensor<T, B1, B2, S>, F>
{
   TensorInvAdjoint(F f = F()) : f_(f) {}

   typedef Tensor::ReducibleTensor<typename result_value<F>::type, B2, B1> result_type;
   typedef Tensor::ReducibleTensor<T, B1, B2, S> const& argument_type;

   result_type operator()(Tensor::ReducibleTensor<T, B1, B2, S> const& x) const
   {
      result_type Result(x.Basis2(), x.Basis1());
      for (typename Tensor::ReducibleTensor<T, B1, B2, S>::const_iterator I = x.begin(); I != x.end(); ++I)
      {
         Result += inv_adjoint(*I, f_);
      }
      return Result;
   }

   F f_;
};

template <typename T, typename B1, typename B2, typename S>
struct InvAdjoint<Tensor::ReducibleTensor<T, B1, B2, S> > : TensorInvAdjoint<Tensor::ReducibleTensor<T, B1, B2, S>,
                                                                       InvAdjoint<typename Tensor::IrredTensor<T, B1, B2, S>::value_type> > {};

} // namespace LinearAlgebra

namespace Tensor
{

template <typename T, typename B1, typename B2, typename S1, typename T2, typename B3, typename B4, typename S2>
struct TensorProd<ReducibleTensor<T, B1, B2, S1>, ReducibleTensor<T2, B3, B4, S2> >
{
   typedef ReducibleTensor<T, B1, B2, S1> first_argument_type;
   typedef ReducibleTensor<T2, B3, B4, S2> second_argument_type;
   typedef IrredTensor<T, B1, B2, S1> Irred1Type;
   typedef IrredTensor<T2, B3, B4, S2> Irred2Type;
   typedef typename TensorProd<Irred1Type, Irred2Type>::result_type IrredResultType;
   typedef ReducibleTensor<typename IrredResultType::value_type,
			   typename IrredResultType::basis1_type,
			   typename IrredResultType::basis2_type,
			   typename IrredResultType::MatrixType> result_type;

   result_type operator()(first_argument_type const& x, second_argument_type const& y) const
   {
      Tensor::ProductBasis<B1, B3> LB1(x.Basis1(), y.Basis1());
      Tensor::ProductBasis<B2, B4> LB2(x.Basis2(), y.Basis2());
      
      result_type Result(LB1.Basis(), LB2.Basis());
      for (typename first_argument_type::const_iterator xI = x.begin(); xI != x.end(); ++xI)
      {
	 for (typename second_argument_type::const_iterator yI = y.begin(); yI != y.end(); ++yI)
	 {
	    QuantumNumbers::QuantumNumberList QL = transform_targets(xI->TransformsAs(), yI->TransformsAs());
	    for (unsigned i = 0; i < QL.size(); ++i)
	    {
	       Result += tensor_prod(*xI, *yI, QL[i]);
	    }
	 }
      }
      return Result;
   }
};

template <typename T, typename B1, typename B2, typename S1, typename T2, typename B3, typename B4, typename S2>
typename TensorProd<ReducibleTensor<T, B1, B2, S1>, ReducibleTensor<T2, B2, B3, S2> >::result_type
tensor_prod(ReducibleTensor<T, B1, B2, S1> const& x, ReducibleTensor<T2, B3, B4, S2> const& y)
{
   return  TensorProd<ReducibleTensor<T, B1, B2, S1>, ReducibleTensor<T2, B3, B4, S2> >()(x,y);
}

} // namespace Tensor

#include "reducible.cc"

#endif
