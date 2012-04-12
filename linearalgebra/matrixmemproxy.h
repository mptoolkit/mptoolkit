/* -*- C++ -*- $Id$

  matrixmemproxy.h

*/

#if !defined(MATRIXMEMPROXY_H_CJKY89YU89P5Y89)
#define MATRIXMEMPROXY_H_CJKY89YU89P5Y89

#include "matrixptriterator.h"
#include "matrixiterators.h"
#include "matrixoperationsbase.h"

namespace LinearAlgebra
{

template <typename T, typename Orientation,
	  typename Size1 = tagVariable, 
	  typename Stride1 = tagVariable, 
	  typename Size2 = tagVariable, 
	  typename Stride2 = tagVariable>
class MatrixMemProxy;


template <typename T, typename Orientation>
class MatrixMemProxy<T, Orientation,
		     tagVariable, tagVariable, tagVariable, tagVariable>
{
   public:
      typedef boost::is_const<T>                    const_proxy;
      typedef boost::mpl::not_<boost::is_const<T> > proxy;

      typedef typename boost::remove_const<T>::type value_type;
      typedef T* pointer;
      typedef T const* const_pointer;
      typedef T& reference;
      typedef T const& const_reference;

      MatrixMemProxy(T* Data, size_type Size1, difference_type Stride1,
		     size_type Size2, difference_type Stride2)
	 : Data_(Data), Size1_(Size1), Size2_(Size2), Stride1_(Stride1), Stride2_(Stride2) 
     	 {}
      //	 { TRACE(this)(Data_)(Size1_)(Stride1_)(Size2_)(Stride2_); }

      MatrixMemProxy& operator=(MatrixMemProxy const& Other)
      {
	 assign(*this, Other);
      }

      template <typename U>
      typename boost::enable_if<is_matrix<U>, MatrixMemProxy&>::type
      operator=(U const& x)
      {
	 // TODO: better temp type
	 typename make_value<U>::type Temp(x);

	 assign(*this, Temp);
	 return *this;
      }

      template <typename U>
      typename boost::enable_if<is_matrix<U>, MatrixMemProxy&>::type
      operator=(NoAliasProxy<U> const& x)
      {
	 assign(*this, x.value());
	 return *this;
      }

      reference operator()(size_type i, size_type j)
	 { return Data_[difference_type(i) * Stride1_ + difference_type(j) * Stride2_]; }

      size_type size1() const { return Size1_; }
      size_type size2() const { return Size2_; }

      pointer data() { return Data_; }
      const_pointer data() const { return Data_; }
      difference_type stride1() const { return Stride1_; }
      difference_type stride2() const { return Stride2_; }

   private:
      T* Data_;
      size_type Size1_, Size2_;
      difference_type Stride1_, Stride2_;
};

// interface

template <typename Scalar, typename Orientation>
struct interface<MatrixMemProxy<Scalar, Orientation> >
{
   typedef typename boost::remove_const<Scalar>::type value_type;
   typedef STRIDE_MATRIX(value_type, Orientation, void) type;
};

// iterators

template <typename T, typename Orient,
	  typename Size1, typename Stride1, typename Size2, typename Stride2>
struct Iterate<MatrixMemProxy<T, Orient, Size1, Stride1, Size2, Stride2>&>
{
   typedef MatrixPtrIterator<T, Orient> result_type;
   typedef MatrixMemProxy<T, Orient, Size1, Stride1, Size2, Stride2>& argument_type;
   result_type operator()(argument_type x) const
   {
      //TRACE(x.data());
      return result_type(x.data(), x.size1(), x.stride1(), x.size2(), x.stride2());
   }
};

template <typename T, typename Orient,
	  typename Size1, typename Stride1, typename Size2, typename Stride2>
struct Iterate<MatrixMemProxy<T, Orient, Size1, Stride1, Size2, Stride2> >
{
   typedef MatrixPtrIterator<T const, Orient> result_type;
   typedef MatrixMemProxy<T, Orient, Size1, Stride1, Size2, Stride2> argument_type;
   result_type operator()(argument_type const& x) const
   {
      //      TRACE(x.data());
      return result_type(x.data(), x.size1(), x.stride1(), x.size2(), x.stride2());
   }
};

// transpose

template <typename T, typename Tv, typename Orient, typename Ti>
struct TransposeInterface<MatrixMemProxy<T, Orient>, STRIDE_MATRIX(Tv, Orient, Ti)>
{
   typedef MatrixMemProxy<Tv const, typename SwapOrientation<Orient>::type> result_type;
   typedef MatrixMemProxy<T, Orient> argument_type;
   result_type operator()(argument_type const& m) const 
      { return result_type(data(m), size2(m), stride2(m), size1(m), stride1(m)); }
};

// defaults for some operations

// swap_sort_order

template <typename T, typename Tv, typename Orient, typename Ti>
struct SwapSortOrderInterface<T, STRIDE_MATRIX(Tv, Orient, Ti)>
{
   typedef MatrixMemProxy<Tv const, typename SwapOrientation<Orient>::type> result_type;
   typedef T const& argument_type;
   result_type operator()(argument_type m) const 
      { return result_type(data(m), size1(m), stride1(m), size2(m), stride2(m)); }
};

template <typename T, typename Tv, typename Orient, typename Ti>
struct SwapSortOrderInterface<T&, STRIDE_MATRIX(Tv, Orient, Ti)>
{
   typedef MatrixMemProxy<Tv, typename SwapOrientation<Orient>::type> result_type;
   typedef T& argument_type;
   result_type operator()(argument_type m) const 
      { return result_type(data(m), size1(m), stride1(m), size2(m), stride2(m)); }
};

// hack Real for STRIDE_MATRIX

template <typename T, typename Tv, typename To, typename Ti>
struct RealInterface<T, STRIDE_MATRIX(std::complex<Tv>, To, Ti)>
{
   typedef MatrixMemProxy<Tv const, To> result_type;
   typedef T const& argument_type;
   result_type operator()(T const& x) const 
      { return result_type(reinterpret_cast<Tv const*>(data(x)), 
			   size1(x), 2 * stride1(x),
			   size2(x), 2 * stride2(x)); }
};

template <typename T, typename Tv, typename To, typename Ti>
struct RealInterface<T&, STRIDE_MATRIX(std::complex<Tv>, To, Ti)>
{
   typedef MatrixMemProxy<Tv, To> result_type;
   typedef T& argument_type;
   result_type operator()(T& x) const 
      { return result_type(reinterpret_cast<Tv*>(data(x)), 
			   size1(x), 2 * stride1(x),
			   size2(x), 2 * stride2(x)); }
};

// hack Imag for STRIDE_MATRIX

template <typename T, typename Tv, typename To, typename Ti>
struct ImagInterface<T, STRIDE_MATRIX(std::complex<Tv>, To, Ti)>
{
   typedef MatrixMemProxy<Tv const, To> result_type;
   typedef T const& argument_type;
   result_type operator()(T const& x) const 
      { return result_type(reinterpret_cast<Tv const*>(data(x))+1, 
			   size1(x), 2 * stride1(x),
			   size2(x), 2 * stride2(x)); }
};

template <typename T, typename Tv, typename To, typename Ti>
struct ImagInterface<T&, STRIDE_MATRIX(std::complex<Tv>, To, Ti)>
{
   typedef MatrixMemProxy<Tv, To> result_type;
   typedef T& argument_type;
   result_type operator()(T& x) const 
      { return result_type(reinterpret_cast<Tv*>(data(x))+1, 
			   size1(x), 2 * stride1(x),
			   size2(x), 2 * stride2(x)); }
};




// an experiment to handle fixed vs variable sizes in a sane way

template <typename Sz1>
class MatrixStrideDims1
{
   public:
      MatrixStrideDims1(size_type Size1)
	 : Size1_(Size1) {}

      size_type size1() const { return Size1_; }

   private:
      size_type Size1_;
};

template <int Size1_>
class MatrixStrideDims1<boost::mpl::int_<Size1_> >
{
   public:
      MatrixStrideDims1(size_type Size1)
	 { DEBUG_PRECONDITION_EQUAL(Size1, Size1_); }

      size_type size1() const { return Size1_; }
      static size_type const static_size1 = Size1_;
};



template <typename Sz1, typename Sr1>
class MatrixStrideDims2 : public MatrixStrideDims1<Sz1>
{
   public:
      MatrixStrideDims2(size_type Size1, difference_type Stride1)
	 : MatrixStrideDims1<Sz1>(Size1),
	 Stride1_(Stride1) {}

      difference_type stride1() const { return Stride1_; }

   private:
      difference_type Stride1_;
};

template <typename Sz1, int Stride1_>
class MatrixStrideDims2<Sz1, boost::mpl::int_<Stride1_> > : public MatrixStrideDims1<Sz1>
{
   public:
      MatrixStrideDims2(size_type Size1, difference_type Stride1)
	 : MatrixStrideDims1<Sz1>(Size1)
	 { DEBUG_PRECONDITION_EQUAL(Stride1, Stride1_); }

      difference_type stride1() const { return Stride1_; }
      static difference_type const static_stride1 = Stride1_;
};



template <typename Sz1, typename Sr1, typename Sz2>
class MatrixStrideDims3 : public MatrixStrideDims2<Sz1, Sr1>
{
   public:
      MatrixStrideDims3(size_type Size1, difference_type Stride1,
		       size_type Size2)
	 : MatrixStrideDims2<Sz1, Sr1>(Size1, Stride1),
	 Size2_(Size2) {}

      size_type size2() const { return Size2_; }

   private:
      size_type Size2_;
};

template <typename Sz1, typename Sr1, int Size2_>
class MatrixStrideDims3<Sz1, Sr1, boost::mpl::int_<Size2_> > : public MatrixStrideDims2<Sz1, Sr1>
{
   public:
      MatrixStrideDims3(size_type Size1, difference_type Stride1,
		       size_type Size2)
	 : MatrixStrideDims2<Sz1, Sr1>(Size1, Stride1)
	 { DEBUG_PRECONDITION_EQUAL(Size2, Size2_); }

      size_type size2() const { return Size2_; }
      static size_type const static_size2 = Size2_;
};



template <typename Sz1, typename Sr1, typename Sz2, typename Sr2>
class MatrixStrideDims : public MatrixStrideDims3<Sz1, Sr1, Sz2>
{
   public:
      MatrixStrideDims(size_type Size1, difference_type Stride1,
		       size_type Size2, difference_type Stride2)
	 : MatrixStrideDims3<Sz1, Sr1, Sz2>(Size1, Stride1, Size2),
	 Stride2_(Stride2) {}

      difference_type stride2() const { return Stride2_; }

   private:
      difference_type Stride2_;
};

template <typename Sz1, typename Sr1, typename Sz2, int Stride2_>
class MatrixStrideDims<Sz1, Sr1, Sz2, boost::mpl::int_<Stride2_> > 
   : public MatrixStrideDims3<Sz1, Sr1, Sz2>
{
   public:
      MatrixStrideDims(size_type Size1, difference_type Stride1,
		       size_type Size2, difference_type Stride2)
	 : MatrixStrideDims3<Sz1, Sr1, Sz2>(Size1, Stride1, Size2)
	 { DEBUG_PRECONDITION_EQUAL(Stride2, Stride2_); }

      difference_type stride2() const { return Stride2_; }
      static difference_type const static_stride2 = Stride2_;
};

#if 0

// linear stride case
template <typename T, typename Size1, int Stride1_, int Size2_, int Stride2_>
struct MatrixStrideRowIterator<T, Size1, 
			       boost::mpl::int_<Stride1_>, 
			       boost::mpl::int_<Size2_>, 
			       boost::mpl::int_<Stride2_>, 
			       typename boost::enable_if<Stride1_ == Size2_ * Stride2_>::type>
{
   //   typedef 
};


// contiguous stride case
template <typename T, typename Size1, int Stride1_, int Size2_>
struct MatrixStrideRowIterator<T, Size1, 
			       boost::mpl::int_<Stride1_>, 
			       boost::mpl::int_<Size2_>, 
			       boost::mpl::int_<1>, 
			       typename boost::enable_if<Stride1_ == Size2_>::type>
{
   //   typedef 
};

#endif

} // namespace LinearAlgebra

#endif
