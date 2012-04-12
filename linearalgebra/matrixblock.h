// -*- C++ -*- $Id$
//
// Class MatrixBlockRef
//
// This is similar to a MatrixMemProxy, but the header information (size, stride etc)
// is not stored locally, but is instead passed in via a pointer to a HeaderType
// object.  The HeaderType can be anything, but must have free functions
// size(), size1(), size2(), stride1(), stride2() defined on it.
//

#if !defined(MATRIXBLOCK_H_DHJFKHURWEIRY4572893475489YRUI34897)
#define MATRIXBLOCK_H_DHJFKHURWEIRY4572893475489YRUI34897

#include "matrixfwd.h"
#include "matrixinterface.h"
#include "matrixoperations.h"
#include "datablockreference.h"
#include "matrixptriterator.h"
#include "crtp_matrix.h"

namespace LinearAlgebra
{

template <typename Scalar, typename Orientation = RowMajor, typename HeaderType,
          typename Derived = void>
class MatrixBlockRef;

template <typename Scalar, typename Orient, typename HeaderType, typename Derived>
struct MatrixBlockRefDerivedType
{
   typedef Derived type;
};

template <typename Scalar, typename Orient, typename HeaderType>
struct MatrixBlockRefDerivedType<Scalar, Orient, HeaderType, void>
{
   typedef MatrixBlockRef<Scalar, Orient, HeaderType> type;
};

template <typename Scalar, typename Orientation, typename HeaderType, typename Derived>
class MatrixBlockRef 
   : public MatrixBase<typename MatrixBlockRefDerivedType<Scalar, Orientation, HeaderType, Derived>::type>
{
   private:
      typedef MatrixBase
              <
                 typename MatrixBlockRefDerivedType
                 <
                    Scalar
                  , Orientation
                  , HeaderType
                  , Derived
                 >::type
              > 
          base_type;

   public:
      using base_type::operator();

      typedef Scalar value_type;
      typedef Scalar& reference;
      typedef Scalar* pointer;

      typedef Scalar const& const_reference;
      typedef Scalar const* const_pointer;

      typedef HeaderType header_type;

      // No default constructor

      MatrixBlockRef(Scalar* B, MatrixOffsetDimensions const* Dims)
         : Block(B), Dims_(Dims) {}

      MatrixBlockRef(MatrixBlockRef const& V) 
	: Block(V.Block_), Dims_(V.Dims_) {}

      MatrixBlockRef& operator=(MatrixBlockRef const& V)
      {
	 this->assign(V);
	 return *this;
      }

      template <typename U>
      typename boost::enable_if<is_matrix<U>, MatrixBlockRef&>::type
      operator=(U const& x);

      template <typename U>
      typename boost::enable_if<is_matrix<U>, MatrixBlockRef&>::type
      operator=(NoAliasProxy<U> const& x);

      size_type size() const { return size(*Dims_); }
      size_type size1() const { return size1(*Dims_); }
      size_type size2() const { return size2(*Dime_); }

      // members defined in MatrixStride
      Scalar* data() { return Block.get() + Dims_.offset; }
      Scalar const* data() const { return Block_.get() + Dims_.offset; }
      difference_type stride1() const { return stride1(*Dims_); }
      difference_type stride2() const { return stride2(*Dims_); }

      Scalar const& operator()(size_type i, size_type j) const 
      { DEBUG_RANGE_CHECK_OPEN(i, 0U, this->size1()); 
        DEBUG_RANGE_CHECK_OPEN(j, 0U, this->size2());
	return *(this->data() + this->stride1() * i + this->stride2() * j); }

      Scalar& operator()(size_type i, size_type j)
      { DEBUG_RANGE_CHECK_OPEN(i, 0U, this->size1()); 
        DEBUG_RANGE_CHECK_OPEN(j, 0U, this->size2());
        this->cow(); return *(this->data() + this->stride1() * i + this->stride2() * j); }

      difference_type leading_dimension() const 
         { return std::max(this->stride1(), this->stride2()); }

   private:
      Scalar* Block_;
      header_type const* Dims_;

   template <typename U, typename Orient, typename D> friend class MatrixBlockRef;
   //   template <typename U, typename Orient, > friend class MatrixConstRef;
};

// interface

template <typename Scalar, typename Orientation, typename HeaderType>
struct interface<MatrixBlockRef<Scalar, Orientation, HeaderType> >
{
   typedef CONTIGUOUS_MATRIX(Scalar, Orientation, void) type;
   typedef Scalar value_type;
};

// iterators

template <typename Scalar, typename Orient, typename Header, typename Derived>
struct Iterate<MatrixBlockRef<Scalar, Orient, Header, Derived>&>
{
   typedef MatrixPtrIterator<Scalar, Orient> result_type;
   typedef MatrixBlockRef<Scalar, Orient, Header, Derived>& argument_type;
   result_type operator()(argument_type x) const
   {
      return result_type(x.data(), x.size1(), x.stride1(), x.size2(), x.stride2());
   }
};

template <typename Scalar, typename Orient, typename Header, typename Derived>
struct Iterate<MatrixBlockRef<Scalar, Orient, Header, Derived> >
{
   typedef MatrixPtrIterator<Scalar const, Orient> result_type;
   typedef MatrixBlockRef<Scalar, Orient, Header, Derived> argument_type;
   result_type operator()(argument_type const& x) const
   {
      return result_type(x.data(), x.size1(), x.stride1(), x.size2(), x.stride2());
   }
};

} // namespace LinearAlgebra

#include "matrixblock.cc"

#endif
