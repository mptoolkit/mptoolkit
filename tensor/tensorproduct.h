// -*- C++ -*- $Id$

#if !defined(TENSORPRODUCT_H_JHCIUWHFIUH98743Y9843YP9)
#define TENSORPRODUCT_H_JHCIUWHFIUH98743Y9843YP9

#include "tensor.h"
#include "linearalgebra/matrix.h"
#include "pstream/pstream.h"
#include <boost/tuple/tuple.hpp>

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
      typedef std::list<target_type>                TargetListType;
      typedef LinearAlgebra::Vector<source_type>    SourceListType;

   public:
      typedef TargetListType::const_iterator const_iterator;

      left_basis_type const& Left() const { return Left_; }
      right_basis_type const& Right() const { return Right_; }

      basis_type const& Basis() const { return Basis_; }

      size_type size() const { return Basis_.size(); }

      bool is_null() const { return Basis_.is_null(); }

      const_iterator begin(int s1, int s2) const { return TransformData_(s1,s2).begin(); }
      const_iterator end(int s1, int s2) const { return TransformData_(s1,s2).end(); }

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

      LinearAlgebra::Matrix<TargetListType> TransformData_;
      SourceListType ReverseMapping_;
};

template <>
class ProductBasis<BasisList, BasisList>
   : public ProductBasisBase<BasisList, BasisList, BasisList>
{
   private:
      typedef ProductBasisBase<BasisList, BasisList, BasisList> base_type;

      ProductBasis(base_type const& x) : base_type(x) {}

   public:
      ProductBasis() {}

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
      ProductBasis() {}
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
      ProductBasis() {}
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
      ProductBasis() {}
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
IrredTensor<typename LinearAlgebra::result_value<ProductFunctor>::type, 
            typename ProductBasis<B1, B3>::basis_type,
            typename ProductBasis<B2, B4>::basis_type>
tensor_prod(IrredTensor<T1, B1, B2, S1> const& ML, IrredTensor<T2, B3, B4, S2> const& MR, 
	    ProductBasis<B1, B3> const& PBasis1, 
	    ProductBasis<B2, B4> const& PBasis2, QuantumNumber q,
	    ProductFunctor TensorProd)
{
   DEBUG_PRECONDITION_EQUAL(ML.Basis1(), PBasis1.Left());
   DEBUG_PRECONDITION_EQUAL(MR.Basis1(), PBasis1.Right());
   DEBUG_PRECONDITION_EQUAL(ML.Basis2(), PBasis2.Left());
   DEBUG_PRECONDITION_EQUAL(MR.Basis2(), PBasis2.Right());

   if (q.is_null())
   {
      QuantumNumbers::QuantumNumberList QL = 
	 transform_targets(ML.TransformsAs(), MR.TransformsAs());
      CHECK(QL.size() == 1)("Transform product is not specified and not unique")
	 (ML.TransformsAs())(MR.TransformsAs())(QL);
      q = QL[0];
   }

   typedef IrredTensor<typename  LinearAlgebra::result_value<ProductFunctor>::type,
      typename ProductBasis<B1, B3>::basis_type,
      typename ProductBasis<B2, B4>::basis_type> ResultType;

   ResultType Result(PBasis1.Basis(), PBasis2.Basis(), q);

   typedef typename const_iterator<IrredTensor<T1, B1, B2, S1> >::type Liter;
   typedef typename const_iterator<IrredTensor<T2, B3, B4, S2> >::type Riter;

   typedef typename const_inner_iterator<IrredTensor<T1, B1, B2, S1> >::type Linner;
   typedef typename const_inner_iterator<IrredTensor<T2, B3, B4, S2> >::type Rinner;

   for (Liter iL = iterate(ML); iL; ++iL)
   {
      for (Riter iR = iterate(MR); iR; ++iR)
      {

         for (Linner jL = iterate(iL); jL; ++jL)
         {
            for (Rinner jR = iterate(iR); jR; ++jR)
            {

               typename ProductBasis<B1,B3>::const_iterator TiIter, 
                  TiEnd = PBasis1.end(jL.index1(), jR.index1());
               for (TiIter = PBasis1.begin(jL.index1(), jR.index1()); TiIter != TiEnd; ++TiIter)
               {
		  typename ProductBasis<B2,B4>::const_iterator TjIter, 
		     TjEnd = PBasis2.end(jL.index2(), jR.index2());
		  for (TjIter = PBasis2.begin(jL.index2(), jR.index2()); TjIter != TjEnd; 
		       ++TjIter)
		  {
                     if (is_transform_target(q, PBasis2[*TjIter], PBasis1[*TiIter]))
                     {
                        double Coeff = tensor_coefficient(PBasis1, PBasis2, 
                                                          ML.TransformsAs(), MR.TransformsAs(), q,
                                                          jL.index1(), jR.index1(), *TiIter, 
                                                          jL.index2(), jR.index2(), *TjIter);
                        
                        if (LinearAlgebra::norm_2(Coeff) > 1E-14)
                           set_new_element(Result.data(), 
                                           *TiIter, *TjIter, 
                                           Coeff * TensorProd(*jL, *jR));
                        else if (Coeff != 0)
                        {
                           DEBUG_WARNING("Coeff is small")(Coeff);
                        }
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
IrredTensor<typename LinearAlgebra::result_value<LinearAlgebra::DirectProduct<T1, T2> >::type, 
            typename ProductBasis<B1, B3>::basis_type,
            typename ProductBasis<B2, B4>::basis_type>
tensor_prod(IrredTensor<T1, B1, B2, S1> const& ML, IrredTensor<T2, B3, B4, S2> const& MR, 
	    ProductBasis<B1,B3> const& PBasis1, ProductBasis<B2,B4> const& PBasis2, 
	    QuantumNumber const& q = QuantumNumber())
{
   return tensor_prod(ML, MR, PBasis1, PBasis2,  q, LinearAlgebra::DirectProduct<T1, T2>());
}

template <typename T1, typename B1, typename B2, typename S1, 
          typename T2, typename B3, typename B4, typename S2>
inline
IrredTensor<typename LinearAlgebra::result_value<LinearAlgebra::DirectProduct<T1, T2> >::type, 
            typename ProductBasis<B1, B2>::basis_type,
            typename ProductBasis<B1, B2>::basis_type>
tensor_prod(IrredTensor<T1, B1, B1, S1> const& ML, IrredTensor<T2, B2, B2, S2> const& MR, 
	    ProductBasis<B1, B2> const& PBasis, 
	    QuantumNumber const& q = QuantumNumber())
{
   return tensor_prod(ML, MR, PBasis, PBasis,  q, LinearAlgebra::DirectProduct<T1, T2>());
}

template <typename T1, typename B1, typename B2, typename S1, 
          typename T2, typename B3, typename B4, typename S2>
inline
IrredTensor<typename LinearAlgebra::result_value<LinearAlgebra::DirectProduct<T1, T2> >::type, 
            typename ProductBasis<B1, B3>::basis_type,
            typename ProductBasis<B2, B4>::basis_type>
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
IrredTensor<typename LinearAlgebra::result_value<LinearAlgebra::DirectProduct<T1, T2> >::type, 
            typename ProductBasis<B1, B3>::basis_type,
            typename ProductBasis<B2, B4>::basis_type>
tensor_prod(IrredTensor<T1, B1, B2, S1> const& ML, IrredTensor<T2, B3, B4, S2> const& MR)
{
   return tensor_prod(ML, MR, 
                      ProductBasis<B1, B3>(ML.Basis1(), MR.Basis1()), 
		      ProductBasis<B2, B4>(ML.Basis2(), MR.Basis2()),
		      QuantumNumber());
}

//
// TensorProd
//
// A binary functor to calculate the tensor product of two operators, given
// the target quantum number and the product basis types as fixed.
//
// TL and TR are assumed to be IrredTensors's.
//

template <typename TL, typename TR>
struct TensorProd
{
};

template <typename TL, typename B1, typename B2, typename SL, 
          typename TR, typename C1, typename C2, typename SR>
struct TensorProd<IrredTensor<TL, B1, B2, SL>, IrredTensor<TR, C1, C2, SR> >
{
   TensorProd(ProductBasis<B1, C1> const& b1, ProductBasis<B2, C2> const& b2, QuantumNumber q) 
      : Basis1_(b1), Basis2_(b2), q_(q) {}

   typedef IrredTensor<TL, B1, B2, SL> const& first_argument_type;
   typedef IrredTensor<TR, C1, C2, SR> const& second_argument_type;

   typedef IrredTensor<typename LinearAlgebra::result_value<
      LinearAlgebra::DirectProduct<TL, TR> >::type,
                       typename ProductBasis<B1, C1>::basis_type,
                       typename ProductBasis<B2, C2>::basis_type> result_type;

   result_type operator()(first_argument_type A, second_argument_type B) const
   {
      return tensor_prod(A, B, Basis1_, Basis2_, q_);
   }

   ProductBasis<B1, C1> Basis1_;
   ProductBasis<B2, C2> Basis2_;
   QuantumNumber q_;
};

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

   typedef typename const_iterator<IrredTensor<T, B1, B2, S> >::type iter;
   typedef typename const_inner_iterator<IrredTensor<T, B1, B2, S> >::type inner;

   typedef typename  ProductBasis<B3, B4>::const_iterator basis1_iter;
   typedef typename  ProductBasis<B5, B6>::const_iterator basis2_iter;

   std::map<PartialProdIndex, T> Result;

   for (iter I = iterate(x); I; ++I)
   {
      for (inner J = iterate(I); J; ++J)
      {
         QuantumNumber q1 = Basis1[J.index1()];
         QuantumNumber q2 = Basis2[J.index2()];

         int Left1, Left2, Right1, Right2;
         boost::tie(Left1,Right1) = Basis1.rmap(J.index1());
         boost::tie(Left2,Right2) = Basis2.rmap(J.index2());

         QuantumNumber qLeft1 = Basis1.Left()[Left1];
         QuantumNumber qRight1 = Basis1.Right()[Right1];
         QuantumNumber qLeft2 = Basis2.Left()[Left2];
         QuantumNumber qRight2 = Basis2.Right()[Right2];

	 QuantumNumberList qLeft = inverse_transform_targets(qLeft2, qLeft1);
	 QuantumNumberList qRight = inverse_transform_targets(qRight2, qRight1);

	 for (QuantumNumberList::const_iterator qL = qLeft.begin(); qL != qLeft.end(); ++qL)
	 {
	    for (QuantumNumberList::const_iterator qR = qRight.begin(); qR != qRight.end(); ++qR)
	    {
	       if (!is_transform_target(*qL, *qR, x.TransformsAs()))
		  continue;

	       double Coeff = inverse_tensor_coefficient(qLeft2, qRight2, q2,
							 *qL, *qR, x.TransformsAs(),
							 qLeft1, qRight1, q1);

	       Result[PartialProdIndex(*qL, Left1, Left2, *qR, Right1, Right2)] += Coeff * (*J);
	    }
	 }
      }
   }

   return Result;
}

} // namespace Tensor

#include "tensorproduct.cc"

#endif
