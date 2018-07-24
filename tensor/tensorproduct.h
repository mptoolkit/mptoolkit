// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/tensorproduct.h
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

#if !defined(TENSORPRODUCT_H_JHCIUWHFIUH98743Y9843YP9)
#define TENSORPRODUCT_H_JHCIUWHFIUH98743Y9843YP9

#include "tensor.h"
#include "pstream/pstream.h"
#include <tuple>
#include "blas/matrix.h"
#include "blas/functors.h"
#include "blas/number_traits.h"

namespace Tensor
{

template <typename BL, typename BR = BL>
class ProductBasis;

template <typename BL, typename BR>
PStream::opstream& operator<<(PStream::opstream& out, ProductBasis<BL, BR> const& B);

template <typename BL, typename BR>
PStream::ipstream& operator>>(PStream::ipstream& in, ProductBasis<BL, BR>& B);

template <typename BL, typename BR, typename Base>
class ProductBasisBase
{
   public:
      typedef BL                              left_basis_type;
      typedef BR                              right_basis_type;
      typedef Base                            basis_type;
      typedef typename basis_type::value_type value_type;

      // in the 'forward' direction, this maps a pair (s1, s2) onto a container of target subspaces
      typedef int target_type;
      // in the opposite direction, we map a subspace into a pair of original subspaces
      typedef std::pair<int, int> source_type;

   private:
      typedef std::list<target_type>      TargetListType;
      typedef std::vector<source_type>    SourceListType;

   public:
      typedef TargetListType::const_iterator const_iterator;

      ProductBasisBase(ProductBasisBase&& Other) = default;

      left_basis_type const& Left() const { return Left_; }
      right_basis_type const& Right() const { return Right_; }

      basis_type const& Basis() const { return Basis_; }

      int size() const { return Basis_.size(); }

      bool is_null() const { return Basis_.is_null(); }

      const_iterator begin(int s1, int s2) const { return TransformData_(s1,s2).begin(); }
      const_iterator end(int s1, int s2) const { return TransformData_(s1,s2).end(); }

      TargetListType const& operator()(int s1, int s2) const
      {
	 return TransformData_(s1,s2);
      }

      const_iterator begin(source_type const& s12) const
        { return this->begin(s12.first, s12.second); }

      const_iterator end(source_type const& s12) const
        { return this->end(s12.first, s12.second); }

      source_type rmap(int s) const { return ReverseMapping_[s]; }

      value_type const& operator[](int s) const { return Basis_[s]; }

   protected:
      ProductBasisBase() {}
      ProductBasisBase(left_basis_type const& B1, right_basis_type const& B2);
      ProductBasisBase(left_basis_type const& B1, right_basis_type const& B2,
                       QuantumNumber const& q);

      // Named constructor for making a basis for a lower-triangular operator that projects onto some irreducible
      // component in the final element.  This is the same as a normal product basis between Basis1_ and Basis2_,
      // except that for the last element in Basis1_ and Basis2_, we take only the single projection q.
      static
      ProductBasisBase
      MakeTriangularProjected(left_basis_type const& Basis1_,
                              right_basis_type const& Basis2_,
                              QuantumNumbers::QuantumNumber const& q);

      basis_type Basis_;
      left_basis_type Left_;
      right_basis_type Right_;

      blas::Matrix<TargetListType> TransformData_;
      SourceListType ReverseMapping_;
};

template <>
class ProductBasis<BasisList, BasisList>
   : public ProductBasisBase<BasisList, BasisList, BasisList>
{
   private:
      typedef ProductBasisBase<BasisList, BasisList, BasisList> base_type;

      ProductBasis(base_type&& x) : base_type(std::move(x)) {}

   public:
      //      ProductBasis() {}

      ProductBasis(left_basis_type const& B1, right_basis_type const& B2)
         : base_type(B1, B2) {}

      // A product basis that projects onto a single quantum number only
      ProductBasis(left_basis_type const& B1, right_basis_type const& B2,
                   QuantumNumber const& q)
         : base_type(B1, B2, q) {}

      // Named constructor for making a basis for a lower-triangular operator that projects onto some irreducible
      // component in the final element.  This is the same as a normal product basis between Basis1_ and Basis2_,
      // except that for the last element in Basis1_ and Basis2_, we take only the single projection q.
      static
      ProductBasis MakeTriangularProjected(left_basis_type Basis1_, right_basis_type Basis2_,
                                           QuantumNumbers::QuantumNumber const& q);

   friend PStream::opstream& operator<< <>(PStream::opstream& out, ProductBasis const& B);
   friend PStream::ipstream& operator>> <>(PStream::ipstream& in, ProductBasis& B);
};

template <>
class ProductBasis<VectorBasis, VectorBasis>
   : public ProductBasisBase<VectorBasis, VectorBasis, VectorBasis>
{
   private:
      typedef ProductBasisBase<VectorBasis, VectorBasis, VectorBasis> base_type;

   public:
      //      ProductBasis() {}
      ProductBasis(left_basis_type const& B1, right_basis_type const& B2);
      ProductBasis(left_basis_type const& B1, right_basis_type const& B2, QuantumNumber const& q);

      // forwarders to VectorBasis
      int total_dimension() const { return Basis_.total_dimension(); }
      int total_degree() const { return Basis_.total_degree(); }
      int dim(int s) const { return Basis_.dim(s); }

   friend PStream::opstream& operator<< <>(PStream::opstream& out, ProductBasis const& B);
   friend PStream::ipstream& operator>> <>(PStream::ipstream& in, ProductBasis& B);
};

template <>
class ProductBasis<BasisList, VectorBasis>
   : public ProductBasisBase<BasisList, VectorBasis, VectorBasis>
{
   private:
      typedef ProductBasisBase<BasisList, VectorBasis, VectorBasis> base_type;

   public:
      //      ProductBasis() {}
      ProductBasis(left_basis_type const& B1, right_basis_type const& B2);
      ProductBasis(left_basis_type const& B1, right_basis_type const& B2, QuantumNumber const& q);

      // forwarders to VectorBasis
      int total_dimension() const { return Basis_.total_dimension(); }
      int total_degree() const { return Basis_.total_degree(); }
      int dim(int s) const { return Basis_.dim(s); }

   friend PStream::opstream& operator<< <>(PStream::opstream& out, ProductBasis const& B);
   friend PStream::ipstream& operator>> <>(PStream::ipstream& in, ProductBasis& B);
};

template <>
class ProductBasis<VectorBasis, BasisList>
   : public ProductBasisBase<VectorBasis, BasisList, VectorBasis>
{
   private:
      typedef ProductBasisBase<VectorBasis, BasisList, VectorBasis> base_type;

   public:
      //      ProductBasis() {}
      ProductBasis(left_basis_type const& B1, right_basis_type const& B2);
      ProductBasis(left_basis_type const& B1, right_basis_type const& B2, QuantumNumber const& q);

      // forwarders to VectorBasis
      int total_dimension() const { return Basis_.total_dimension(); }
      int total_degree() const { return Basis_.total_degree(); }
      int dim(int s) const { return Basis_.dim(s); }

   friend PStream::opstream& operator<< <>(PStream::opstream& out, ProductBasis const& B);
   friend PStream::ipstream& operator>> <>(PStream::ipstream& in, ProductBasis& B);
};

template <typename B1, typename B2>
inline
ProductBasis<B1, B2>
make_product_basis(B1 const& b1, B2 const& b2)
{
   return ProductBasis<B1, B2>(b1, b2);
}

template <typename B1, typename B2>
inline
ProductBasis<B1, B2>
make_product_basis(B1 const& b1, B2 const& b2, QuantumNumber const& q)
{
   return ProductBasis<B1, B2>(b1, b2, q);
}

template <typename B1a, typename B1b, typename B2a, typename B2b>
inline
double
tensor_coefficient(ProductBasis<B1a, B1b> const& P1, ProductBasis<B2a, B2b> const& P2,
                   QuantumNumber const& k1, QuantumNumber const& k2, QuantumNumber const& k,
                   int s1prime, int s2prime, int sprime, int s1, int s2, int s)
{
   return  tensor_coefficient(P2.Left()[s1],      P2.Right()[s2],      P2[s],
                              k1,                 k2,                  k,
                              P1.Left()[s1prime], P1.Right()[s2prime], P1[sprime]);
}

template <typename T1, typename B1, typename B2, typename S1,
          typename T2, typename B3, typename B4, typename S2, typename ProductFunctor>
	  IrredTensor<blas::remove_proxy_t<decltype(std::declval<ProductFunctor>()(std::declval<T1>(), std::declval<T2>()))>,
	   typename ProductBasis<B1, B3>::basis_type,
	   typename ProductBasis<B2, B4>::basis_type>
tensor_prod(IrredTensor<T1, B1, B2, S1> const& ML, IrredTensor<T2, B3, B4, S2> const& MR,
            ProductBasis<B1, B3> const& PBasis1,
            ProductBasis<B2, B4> const& PBasis2, QuantumNumber q,
            ProductFunctor ProdFunctor = blas::Multiplication())
{
   DEBUG_PRECONDITION_EQUAL(ML.Basis1(), PBasis1.Left());
   DEBUG_PRECONDITION_EQUAL(MR.Basis1(), PBasis1.Right());
   DEBUG_PRECONDITION_EQUAL(ML.Basis2(), PBasis2.Left());
   DEBUG_PRECONDITION_EQUAL(MR.Basis2(), PBasis2.Right());

   using blas::norm_frob;

   if (q.is_null())
   {
      QuantumNumbers::QuantumNumberList QL =
         transform_targets(ML.TransformsAs(), MR.TransformsAs());
      CHECK(QL.size() == 1)("Transform product is not specified and not unique")
         (ML.TransformsAs())(MR.TransformsAs())(QL);
      q = QL[0];
   }

   using result_value = blas::remove_proxy_t<decltype(std::declval<ProductFunctor>()(std::declval<T1>(), std::declval<T2>()))>;

   using ResultType = IrredTensor<result_value,
      typename ProductBasis<B1, B3>::basis_type,
      typename ProductBasis<B2, B4>::basis_type>;

   ResultType Result(PBasis1.Basis(), PBasis2.Basis(), q);

   for (auto const& rML : ML)
   {
      for (auto const& rMR : MR)
      {
	 for (auto const& cML : rML)
	 {
	    for (auto const& cMR : rMR)
	    {
	       for (auto const& b1 : PBasis1(rML.row(), rMR.row()))
	       {
		  for (auto const& b2 : PBasis2(cML.col(), cMR.col()))
		  {
		     if (is_transform_target(q, PBasis2[b2], PBasis1[b1]))
		     {
                        auto Coeff = tensor_coefficient(PBasis1, PBasis2,
							ML.TransformsAs(), MR.TransformsAs(), q,
							rML.row(), rMR.row(), b1,
							cML.col(), cMR.col(), b2);
			
                        if (norm_frob(Coeff) > 1E-14)
			   Result.insert(b1, b2, Coeff * ProdFunctor(cML.value, cMR.value));
		     }
                  }
               }
            }
         }
      }
   }
   return Result;
}

template <typename T1, typename B1, typename B2, typename S1,
          typename T2, typename B3, typename B4, typename S2>
inline
auto
tensor_prod(IrredTensor<T1, B1, B2, S1> const& ML, IrredTensor<T2, B3, B4, S2> const& MR,
            ProductBasis<B1,B3> const& PBasis1, ProductBasis<B2,B4> const& PBasis2,
            QuantumNumber const& q = QuantumNumber())
{
   return tensor_prod(ML, MR, PBasis1, PBasis2,  q, blas::Multiplication());
}

template <typename T1, typename B1, typename B2, typename S1,
          typename T2, typename B3, typename B4, typename S2>
inline
auto
tensor_prod(IrredTensor<T1, B1, B1, S1> const& ML, IrredTensor<T2, B2, B2, S2> const& MR,
            ProductBasis<B1, B2> const& PBasis,
            QuantumNumber const& q = QuantumNumber())
{
   return tensor_prod(ML, MR, PBasis, PBasis,  q, blas::Prod());
}

template <typename T1, typename B1, typename B2, typename S1,
          typename T2, typename B3, typename B4, typename S2>
inline
auto
tensor_prod(IrredTensor<T1, B1, B2, S1> const& ML,
            IrredTensor<T2, B3, B4, S2> const& MR,
            QuantumNumber const& q)
{
   return tensor_prod(ML, MR,
                      ProductBasis<B1, B3>(ML.Basis1(), MR.Basis1()),
                      ProductBasis<B2, B4>(ML.Basis2(), MR.Basis2()),
                      q);
}

template <typename T1, typename B1, typename B2, typename S1,
          typename T2, typename B3, typename B4, typename S2>
inline
auto
tensor_prod(IrredTensor<T1, B1, B2, S1> const& ML, IrredTensor<T2, B3, B4, S2> const& MR)
{
   return tensor_prod(ML, MR,
                      ProductBasis<B1, B3>(ML.Basis1(), MR.Basis1()),
                      ProductBasis<B2, B4>(ML.Basis2(), MR.Basis2()),
                      QuantumNumber());
}

//
// partial transpose
//
// Given a tensor x^[\gamma]_{ab}, derived from a tensor product:
// a = (i j), b = (k l).
// Give also the cross product basis: alpha = (i \bar{k}), beta = (j \bar{l}).
// Then the partial transpose of x is a tensor
// y^[\gamma]_{\bar{alpha} beta}
// which we can then decompose into a product
// y^[\gamma]_{\bar{alpha} beta} = E^{alpha}_{ik} . F^{beta}_{jl} [no coupling coefficent here!]
// In particular, E and F are given by the singular decomposition of y.
//
// Currently, this only works properly when T is a scalar type;
// for a non-scalar T we would need to figure out how to
// do the nested partial_transpose.
//

template <typename T, typename B1, typename B2, typename S,
          typename B3, typename B4, typename B5, typename B6>
IrredTensor<T,
            typename ProductBasis<B3, B4>::basis_type,
            typename ProductBasis<B5, B6>::basis_type>
partial_transpose(IrredTensor<T, B1, B2, S> const& x,
                  ProductBasis<B3, B5> const& a,
                  ProductBasis<B4, B6> const& b,
                  ProductBasis<B3, B4> const& alpha,
                  ProductBasis<B5, B6> const& beta);

//
// swap_product_basis
// Constructs a unitary transformation that maps between a product
// basis ProductBasis<B1,B2> and ProductBasis<B2,B1>
//
// ** 2016-07-08: untested **

template <typename B1, typename B2>
IrredTensor<double,
            typename ProductBasis<B1,B2>::basis_type,
            typename ProductBasis<B2,B1>::basis_type>
swap_product_basis(ProductBasis<B1,B2> const& a,
                   ProductBasis<B2,B1> const& b);

template <typename B1, typename B2>
IrredTensor<double,
            typename ProductBasis<B1,B2>::basis_type,
            typename ProductBasis<B2,B1>::basis_type>
swap_product_basis(ProductBasis<B1,B2> const& a,
                   ProductBasis<B2,B1> const& b)
{
   IrredTensor<double,
      typename ProductBasis<B1,B2>::basis_type,
      typename ProductBasis<B2,B1>::basis_type>
      Result(a, b, QuantumNumber(a.GetSymmetryList()));

   CHECK_EQUAL(a.Left(), b.Right());
   CHECK_EQUAL(a.Right(), b.Left());

   for (int i = 0; i < a.size(); ++i)
   {
      std::pair<int,int> p = a.rmap(i);
      // Find a mapping in b(p.second, p.first) that has the same quantum number as i
      typename ProductBasis<B2,B1>::const_iterator J = b.begin(p.second, p.first);
      typename ProductBasis<B2,B1>::const_iterator Jend = b.end(p.second, p.first);
      while (J != Jend && b[*J] != a[i])
         ++J;
      if (J != Jend)
      {
         Result(i, *J) = conj_phase(a.Left()[p.first], a[i], a.Right()[p.second]);
      }
      else
      {
         WARNING("swap_product_basis")("component is missing in the right-hand basis");
      }
   }
   return Result;
}

//
// decompose_tensor_prod
//
// Decompose a tensor product matrix element <(j_1' j_2') j' || M^k || (j_1 j_2) j>
// into < j_1' || M^{k_1} || j_1 >  < j_2' || M^{k_2} || j_2 >

struct PartialProdIndex
{
   typedef QuantumNumbers::QuantumNumber QuantumNumber;

   QuantumNumber qLeft;
   int Left1, Left2;

   QuantumNumber qRight;
   int Right1, Right2;

   PartialProdIndex(QuantumNumber const& qLeft_, int Left1_, int Left2_,
                    QuantumNumber const& qRight_, int Right1_, int Right2_)
      : qLeft(qLeft_), Left1(Left1_), Left2(Left2_),
        qRight(qRight_), Right1(Right1_), Right2(Right2_)
   {}
};

inline
std::ostream& operator<<(std::ostream& out, PartialProdIndex const& p)
{
   out << p.qLeft << ": (" << p.Left1 << "," << p.Left2 << ")  "
       << p.qRight << ": (" << p.Right1 << "," << p.Right2 << ")  ";
   return out;
}

inline
bool operator<(PartialProdIndex const& x, PartialProdIndex const& y)
{
   return (x.qLeft < y.qLeft ||
           (x.qLeft == y.qLeft && (x.Left1 < y.Left1 ||
           (x.Left1 == y.Left1 && (x.Left2 < y.Left2 ||
           (x.Left2 == y.Left2 && (x.qRight < y.qRight ||
           (x.qRight == y.qRight && (x.Right1 < y.Right1 ||
           (x.Right1 == y.Right1 && x.Right2 < y.Right2))))))))));
}

template <typename T, typename B1, typename B2, typename S,
          typename B3, typename B4, typename B5, typename B6>
std::map<PartialProdIndex, T>
decompose_tensor_prod(IrredTensor<T, B1, B2, S> const& x,
                      ProductBasis<B3, B5> const& Basis1,
                      ProductBasis<B4, B6> const& Basis2)
{
   DEBUG_PRECONDITION_EQUAL(x.Basis1(), Basis1.Basis());
   DEBUG_PRECONDITION_EQUAL(x.Basis2(), Basis2.Basis());

   using QuantumNumbers::QuantumNumberList;

   typedef typename  ProductBasis<B3, B4>::const_iterator basis1_iter;
   typedef typename  ProductBasis<B5, B6>::const_iterator basis2_iter;

   std::map<PartialProdIndex, T> Result;

   for (auto const& rx : x)
   {
      for (auto const& cx : rx)
      {
         QuantumNumber q1 = Basis1[rx.row()];
         QuantumNumber q2 = Basis2[cx.col()];

         int Left1, Left2, Right1, Right2;
         std::tie(Left1,Right1) = Basis1.rmap(rx.row());
         std::tie(Left2,Right2) = Basis2.rmap(cx.col());

         QuantumNumber qLeft1 = Basis1.Left()[Left1];
         QuantumNumber qRight1 = Basis1.Right()[Right1];
         QuantumNumber qLeft2 = Basis2.Left()[Left2];
         QuantumNumber qRight2 = Basis2.Right()[Right2];

         QuantumNumberList qLeft = inverse_transform_targets(qLeft2, qLeft1);
         QuantumNumberList qRight = inverse_transform_targets(qRight2, qRight1);

	 for (auto const& qL : qLeft)
	 {
	    for (auto const& qR : qRight)
	    {
               if (!is_transform_target(qL, qR, x.TransformsAs()))
                  continue;

               double Coeff = inverse_tensor_coefficient(qLeft2, qRight2, q2,
                                                         qL, qR, x.TransformsAs(),
                                                         qLeft1, qRight1, q1);

               Result[PartialProdIndex(qL, Left1, Left2, qR, Right1, Right2)] += Coeff * cx.value;
            }
         }
      }
   }

   return Result;
}

} // namespace Tensor

#include "tensorproduct.icc"

#endif
