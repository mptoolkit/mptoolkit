// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/reducible.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

/*
  Reducible tensor.  This is represented as a sum of irriducible tensors,
  indexed by the quantum number.  Conceptually, it is similar to an
  MPS or a 3-index tensor, but the third index is constrained to be
  a simple list of distinct quantum numbers.

*/

#if !defined(MPTOOLKIT_TENSOR_REDUCIBLE_H)
#define MPTOOLKIT_TENSOR_REDUCIBLE_H

#include "tensor.h"
#include "tensorproduct.h"
#include <set>
#include "common/set_operations.h"
#include <boost/iterator/transform_iterator.hpp>
#include "blas/number_traits.h"

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
    typename Structure = DefaultStructure>
class ReducibleTensor;

template <typename T, typename B1, typename B2, typename S>
PStream::opstream& operator<<(PStream::opstream& out, ReducibleTensor<T, B1, B2, S> const& x);

template <typename T, typename B1, typename B2, typename S>
PStream::ipstream& operator>>(PStream::ipstream& in, ReducibleTensor<T, B1, B2, S>& x);

template <typename T, typename B1, typename B2, typename S>
ReducibleTensor<T, B1, B2, S>
CoerceSymmetryList(ReducibleTensor<T, B1, B2, S> const& t, SymmetryList const& sl);

template <typename T, typename B1, typename B2, typename S>
void
CoerceSymmetryListInPlace(ReducibleTensor<T, B1, B2, S>& t, SymmetryList const& sl);

template <typename T, typename B1, typename B2, typename S>
inline
ReducibleTensor<T, B1, B2, S>
copy(ReducibleTensor<T, B1, B2, S> const& x);

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

      // need default ctor for persistent streaming
      ReducibleTensor() = default;

      ReducibleTensor(ReducibleTensor const&) = delete;

      ReducibleTensor(ReducibleTensor&& Other) noexcept = default;

      ReducibleTensor& operator=(ReducibleTensor&& Other) noexcept = default;

      ReducibleTensor(basis1_type const& Basis);
      ReducibleTensor(basis1_type const& Basis1, basis2_type const& Basis2);

      ReducibleTensor(basis1_type const& Basis1, basis2_type const& Basis2,
		      data_type&& Data)
	 : Basis1_(Basis1), Basis2_(Basis2), data_(std::move(Data)) {}


      template <typename U, typename US>
      ReducibleTensor(IrredTensor<U, basis1_type, basis2_type, US> x);

#if 0
      template <typename U, typename US>
      ReducibleTensor(ReducibleTensor<U, basis1_type, basis2_type, US> const& x);
#endif

      bool empty() const { return data_.empty(); }

      bool is_null() const { return this->empty(); }

      // returns the number of irreducible components
      int size() const { return data_.size(); }

      SymmetryList GetSymmetryList() const { return Basis1_.GetSymmetryList(); }

      int size1() const { return Basis1_.size(); }
      int size2() const { return Basis2_.size(); }

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

      void insert(value_type x);

      void add(value_type x);

      // projects onto the scalar component
      value_type scalar() const;
      value_type& scalar();

      void trim()
      {
         map_iterator I = data_.begin();
         while (I != data_.end())
         {
            if (norm_frob_sq(I->second) == 0)
            {
               map_iterator Temp = I;
               ++I;
               data_.erase(Temp);
            }
            else
               ++I;
         }
      }

      // front element of the container
      value_type const& front() const
      {
	 return data_.begin()->second;
      }

      // variant where we know in advance that the corresponding component
      // exists in the container
      value_type const& project_assert(QuantumNumber const& q) const;

      QuantumNumber const& qn1(int i) const { return Basis1_[i]; }
      QuantumNumber const& qn2(int j) const { return Basis2_[j]; }

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

      static_assert(std::is_nothrow_move_constructible<basis1_type>::value, "");
      static_assert(std::is_nothrow_move_constructible<basis2_type>::value, "");
      static_assert(std::is_nothrow_move_constructible<data_type>::value, "");

   template <typename U, typename B1, typename B2, typename S>
   friend class ReducibleTensor;

   friend PStream::opstream& operator<< <>(PStream::opstream& out, ReducibleTensor const& x);
   friend PStream::ipstream& operator>> <>(PStream::ipstream& in, ReducibleTensor& x);

   friend void CoerceSymmetryListInPlace<>(ReducibleTensor& t, SymmetryList const& sl);
   friend ReducibleTensor CoerceSymmetryList<>(ReducibleTensor const& t, SymmetryList const& sl);

      friend ReducibleTensor copy<T,Basis1T,Basis2T,Structure>(ReducibleTensor const& x);
};

} // namespace Tensor

template <typename T, typename B1, typename B2, typename S>
struct ScalarTypes<Tensor::ReducibleTensor<T, B1, B2, S>> : ScalarTypes<T> {};

namespace Tensor
{

// arithmetic

template <typename T, typename B1, typename B2, typename S>
ReducibleTensor<T, B1, B2, S> operator+(ReducibleTensor<T, B1, B2, S> const& x,
					ReducibleTensor<T, B1, B2, S> const& y)
{
   ReducibleTensor<T, B1, B2, S> Result = copy(x);
   Result += y;
   return Result;
}
template <typename T, typename B1, typename B2, typename S>
ReducibleTensor<T, B1, B2, S> operator-(ReducibleTensor<T, B1, B2, S> const& x,
					ReducibleTensor<T, B1, B2, S> const& y)
{
   ReducibleTensor<T, B1, B2, S> Result = copy(x);
   Result -= y;
   return Result;
}

template <typename T, typename B1, typename B2, typename S>
ReducibleTensor<T, B1, B2, S> operator-(ReducibleTensor<T, B1, B2, S> const& x)
{
   ReducibleTensor<T, B1, B2, S> Result(x.Basis1(), x.Basis2());
   for (auto const& c : x)
   {
      Result.insert(-c);
   }
   return Result;
}

template <typename T, typename B1, typename B2, typename S>
bool is_scalar(ReducibleTensor<T, B1, B2, S> const& x)
{
   return x.size() == 1 && is_scalar(x.front());
}

template <typename T, typename B1, typename B2, typename S>
inline
ReducibleTensor<T, B1, B2, S>
copy(ReducibleTensor<T, B1, B2, S> const& x)
{
   typename ReducibleTensor<T, B1, B2, S>::data_type Copy;
   for (auto const& c : x.data_)
   {
      Copy.insert(std::make_pair(c.first, copy(c.second)));
   }
   return ReducibleTensor<T, B1, B2, S>(x.Basis1(), x.Basis2(), std::move(Copy));
}

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

// text I/O

template <typename T, typename B1, typename B2, typename S>
std::ostream& operator<<(std::ostream& out, ReducibleTensor<T, B1, B2, S> const& Op);

// construction

template <typename T, typename B1, typename B2, typename S>
inline
ReducibleTensor<T, B1, B2, S>
ReducibleTensor<T, B1, B2, S>::make_identity(B1 const& b)
{
   return IrredTensor<T, B1, B2, S>::make_identity(b);
}

// complex conjugation

template <typename T, typename B1, typename B2, typename S>
void
inplace_conj(ReducibleTensor<T, B1, B2, S>& x)
{
   for (auto& c : x)
   {
      inplace_conj(c);
   }
}

template <typename T, typename B1, typename B2, typename S>
ReducibleTensor<T, B2, B1, S>
conj(ReducibleTensor<T, B1, B2, S> const& x)
{
   ReducibleTensor<T, B2, B1, S> Result(x.Basis1(), x.Basis2());
   for (auto const& c : x)
   {
      Result.insert(conj(c));
   }
   return Result;
}

template <typename T, typename B1, typename B2, typename S>
ConjugateProxy<ReducibleTensor<T, B1, B2, S>>
conj(Tensor::ReducibleTensor<T, B1, B2, S> const& x)
{
   return ConjugateProxy<ReducibleTensor<T, B1, B2, S>>(x);
}

// hermitian conjugation

template <typename T, typename B1, typename B2, typename S>
HermitianProxy<ReducibleTensor<T, B1, B2, S>>
herm(ReducibleTensor<T, B1, B2, S> const& x)
{
   return HermitianProxy<ReducibleTensor<T, B1, B2, S>>(x);
}

// flip conjugation
template <typename T, typename B1, typename B2, typename S>
ReducibleTensor<T, B2, B1, S>
flip_conj(ReducibleTensor<T, B1, B2, S> const& x)
{
   ReducibleTensor<T, B2, B1, S> Result(adjoint(x.Basis1()), adjoint(x.Basis2()));
   for (auto const& c : x)
   {
      Result.insert(flip_conj(c));
   }
   return Result;
}

template <typename T, typename B1, typename B2, typename S>
ReducibleTensor<T, B2, B1, S>
flip_conj(ReducibleTensor<T, B1, B2, S>&& x)
{
   ReducibleTensor<T, B2, B1, S> Result(adjoint(x.Basis1()), adjoint(x.Basis2()));
   for (auto&& c : x)
   {
      Result.insert(flip_conj(std::move(c.value)));
   }
   return Result;
}

// trace

template <typename T, typename B1, typename B2, typename S>
auto
trace(Tensor::ReducibleTensor<T, B1, B2, S> const& x)
{
   return trace(project(x,  QuantumNumber(x.GetSymmetryList())));
}

// norm_frob

template <typename T, typename B1, typename B2, typename S>
real_t<T>
norm_frob_sq(Tensor::ReducibleTensor<T, B1, B2, S> const& x)
{
   real_t<T> Result = 0;
   for (auto const& Component : x)
   {
      Result += norm_frob_sq(Component);
   }
   return Result;
}

template <typename T, typename B1, typename B2, typename S>
inline
real_t<T>
norm_frob(Tensor::ReducibleTensor<T, B1, B2, S> const& x)
{
   using std::sqrt;
   return sqrt(norm_frob_sq(x));
}

// inner_prod

template <typename T1, typename B1, typename B2, typename S1, typename T2, typename S2>
scalar_t<decltype(std::declval<T1>()*std::declval<T2>())>
inner_prod(ReducibleTensor<T1, B1, B2, S1> const& x, ReducibleTensor<T2, B1, B2, S2> const& y)
{
   scalar_t<decltype(std::declval<T1>()*std::declval<T2>())> Result = 0.0;
   std::set<QuantumNumber> q1 = x.components(), q2 = y.components();
   std::set<QuantumNumber> Intersect;
   std::set_intersection(q1.begin(), q1.end(), q2.begin(), q2.end(), std::inserter(Intersect, Intersect.begin()));
   for (auto const& q : Intersect)
   {
      Result += inner_prod(x.project_assert(q), y.project_assert(q));
   }
   return Result;
}

// scalar_prod

template <typename T1, typename B1, typename B2, typename S1, typename T2, typename S2>
IrredTensor<blas::remove_proxy_t<decltype(std::declval<T1>()*std::declval<T2>())>,B1, B1>
scalar_prod(ReducibleTensor<T1, B1, B2, S1> const& x, HermitianProxy<ReducibleTensor<T2, B1, B2, S2>> const& y)
{
   IrredTensor<blas::remove_proxy_t<decltype(std::declval<T1>()*std::declval<T2>())>,B1, B1>
      Result(x.Basis1(), y.base().Basis1());
   std::set<QuantumNumber> q1 = x.components(), q2 = y.base().components();
   intersection_iterator<std::set<QuantumNumber> > IEnd = set_intersection_end(q1,q2);
   intersection_iterator<std::set<QuantumNumber> > I = set_intersection_begin(q1,q2);

   while (I != IEnd)
   {
      Result += scalar_prod(x.project_assert(*I), herm(y.base().project_assert(*I)));
      ++I;
   }
   return Result;
}

template <typename T1, typename B1, typename B2, typename S1, typename T2, typename S2>
IrredTensor<blas::remove_proxy_t<decltype(std::declval<T1>()*std::declval<T2>())>,B1, B1>
scalar_prod(HermitianProxy<ReducibleTensor<T1, B1, B2, S1>> const& x, ReducibleTensor<T2, B1, B2, S2> const& y)
{
   IrredTensor<blas::remove_proxy_t<decltype(std::declval<T1>()*std::declval<T2>())>,B1, B1>
      Result(x.base().Basis2(), y.Basis2());
   std::set<QuantumNumber> q1 = x.base().components(), q2 = y.components();
   intersection_iterator<std::set<QuantumNumber> > IEnd = set_intersection_end(q1,q2);
   intersection_iterator<std::set<QuantumNumber> > I = set_intersection_begin(q1,q2);
   while (I != IEnd)
   {
      Result += scalar_prod(herm(x.base().project_assert(*I)), y.project_assert(*I));
      ++I;
   }
   return Result;
}

template <typename T, typename B, typename S>
ReducibleTensor<T, B, B, S>
prod(ReducibleTensor<T, B, B, S> const& x, ReducibleTensor<T, B, B, S> const& y)
{
   using QuantumNumbers::QuantumNumber;
   using QuantumNumbers::QuantumNumberList;
   ReducibleTensor<T, B, B, S> Result(x.Basis1(), y.Basis2());
   for (auto const& Ix : x)
   {
      for (auto const& Iy : y)
      {
	 for (auto const& q : Ix.TransformsAs() * Iy.TransformsAs())
	 {
	    Result.add(prod(Ix, Iy, q));
	 }
      }
   }
   return Result;
}

template <typename T, typename B, typename S>
inline
ReducibleTensor<T, B, B, S>
operator*(ReducibleTensor<T, B, B, S> const& x, ReducibleTensor<T, B, B, S> const& y)
{
   return prod(x,y);
}

// multiplication by scalar
template <typename T, typename B, typename S, typename U>
inline
std::enable_if_t<blas::is_numeric_v<U>, ReducibleTensor<T, B, B, S>>
operator*(ReducibleTensor<T, B, B, S> const& x, U a)
{
   ReducibleTensor<T, B, B, S> Result(copy(x));
   Result *= a;
   return std::move(Result);
}

template <typename T, typename B, typename S, typename U>
inline
std::enable_if_t<blas::is_numeric_v<U>, ReducibleTensor<T, B, B, S>>
operator*(ReducibleTensor<T, B, B, S>&& Result, U a)
{
   Result *= a;
   return std::move(Result);
}


template <typename T, typename B, typename S, typename U>
inline
std::enable_if_t<blas::is_numeric_v<U>, ReducibleTensor<T, B, B, S>>
operator*(U a, ReducibleTensor<T, B, B, S> const& x)
{
   ReducibleTensor<T, B, B, S> Result(copy(x));
   Result *= a;
   return std::move(Result);
}

template <typename T, typename B, typename S, typename U>
inline
std::enable_if_t<blas::is_numeric_v<U>, ReducibleTensor<T, B, B, S>>
operator*(U a, ReducibleTensor<T, B, B, S>&& Result)
{
   Result *= a;
   return std::move(Result);
}

template <typename T, typename B1, typename B2, typename S>
ReducibleTensor<T,B2,B1,S>
adjoint(ReducibleTensor<T, B1, B2, S> const& x)
{
   ReducibleTensor<T,B2,B1,S> Result(x.Basis2(), x.Basis1());
   for (auto const& xc : x)
   {
      Result.insert(adjoint(xc));
   }
}

template <typename T, typename B1, typename B2, typename S>
ReducibleTensor<T,B2,B1,S>
inv_adjoint(ReducibleTensor<T, B1, B2, S> const& x)
{
   ReducibleTensor<T,B2,B1,S> Result(x.Basis2(), x.Basis1());
   for (auto const& xc : x)
   {
      Result.insert(inv_adjoint(xc));
   }
}

#if 0
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
#endif

} // namespace Tensor

#include "reducible.icc"

#endif
