/* -*- C++ -*- $Id$

  matrixbinaryiterator.h

  Created 2005-01-10, Ian McCulloch
*/

#if !defined(MATRIXBINARYITERATOR_H_CHHUIYTH879YO)
#define MATRIXBINARYITERATOR_H_CHHUIYTH879YO

#include "vectortransformiterator.h"  // for operator_arrow_proxy

namespace LinearAlgebra
{


template <typename Iter1, typename Iter2, typename Func,
	  typename Cat1 = typename Iter1::category,
	  typename Cat2 = typename Iter2::category>
class MatrixBinaryInnerIterator;

template <typename Iter1, typename Iter2, typename Func,
	  typename Cat1 = typename Iter1::category,
	  typename Cat2 = typename Iter2::category>
class MatrixBinaryOuterIterator;

template <typename Iter1, typename Iter2, typename Func>
class MatrixBinaryInnerIterator<Iter1, Iter2, Func, vector_iterator_dense, vector_iterator_dense>
{
   public:
      typedef Func functor_type;
      typedef typename make_value<typename functor_type::result_type>::type value_type;
      typedef typename functor_type::result_type reference;
      typedef operator_arrow_proxy<value_type> pointer;
      typedef vector_iterator_dense category;

      typedef Iter1 iterator1_type;
      typedef Iter2 iterator2_type;

      MatrixBinaryInnerIterator() {}

      explicit MatrixBinaryInnerIterator(iterator1_type const& i1, iterator2_type const& i2,
				    functor_type f = functor_type())
	 : i1_(i1), i2_(i2), f_(f) {}

      MatrixBinaryInnerIterator& operator++() 
	 { ++i1_; ++i2_; return *this; }

      MatrixBinaryInnerIterator& operator++(int) 
	 { return MatrixBinaryInnerIterator(i1_++, i2_++, f_); }

      MatrixBinaryInnerIterator& operator--() 
	 { --i1_; --i2_; return *this; }

      MatrixBinaryInnerIterator& operator--(int) 
	 { return MatrixBinaryInnerIterator(i1_--, i2_--, f_); }

      MatrixBinaryInnerIterator& operator+=(difference_type n)
	 { i1_ += n; i2_ += n; return *this; }

      MatrixBinaryInnerIterator& operator-=(difference_type n)
	 { i1_ -= n; i2_ -= n; return *this; }

      size_type index1() const 
      { DEBUG_CHECK_EQUAL(i1_.index1(), i2_.index1()); return i1_.index1(); }

      size_type index2() const 
      { DEBUG_CHECK_EQUAL(i1_.index2(), i2_.index2()); return i1_.index2(); }

      reference operator*() const 
      { 
         DEBUG_CHECK_EQUAL(i1_.index1(), i2_.index1());
         DEBUG_CHECK_EQUAL(i1_.index2(), i2_.index2());
         return f_(*i1_, *i2_); 
      }

      pointer operator->() const 
      { 
         DEBUG_CHECK_EQUAL(i1_.index1(), i2_.index1()); 
         DEBUG_CHECK_EQUAL(i1_.index2(), i2_.index2());
         return pointer(&f_(*i1_, *i2_)); 
      }

      operator bool() const { return i1_; }

      iterator1_type const& iterator1() const { return i1_; }
      iterator2_type const& iterator2() const { return i2_; }

      functor_type const& functor() const { return f_; }

      // derived concepts
      size_type size() const { return i1_.size(); }
      reference operator[](difference_type n) const { return f_(i1_[n], i2_[n]); }

   private:
      iterator1_type i1_;
      iterator2_type i2_;
      functor_type f_;
};
   
template <typename Iter1, typename Iter2, typename Func>
class MatrixBinaryOuterIterator<Iter1, Iter2, Func, vector_iterator_dense, vector_iterator_dense>
{
   public:
      typedef Func functor_type;
      //      typedef typename make_value<typename functor_type::result_type>::type value_type;
      //      typedef typename functor_type::result_type reference;
      //      typedef operator_arrow_proxy<value_type> pointer;
      typedef vector_iterator_dense category;

      typedef typename Iter1::reference reference1;
      typedef typename Iter2::reference reference2;

      typedef typename basic_type<reference1>::type arg1;
      typedef typename basic_type<reference2>::type arg2;

      //      typedef typename reference_to_arg<reference1>::type arg1;
      //      typedef typename reference_to_arg<reference2>::type arg2;

      typedef typename BinaryTransform<arg1, arg2, Func>::result_type reference;

      //   typedef typename boost::mpl::print<arg1>::type d1;
      //   typedef typename boost::mpl::print<arg2>::type d2;

      typedef typename make_value<reference>::type value_type;
      typedef operator_arrow_proxy<value_type> pointer;

      typedef Iter1 iterator1_type;
      typedef Iter2 iterator2_type;

      typedef typename iterator<iterator1_type>::type inner_iterator1;
      typedef typename iterator<iterator2_type>::type inner_iterator2;

      typedef MatrixBinaryInnerIterator<inner_iterator1, inner_iterator2, Func> iterator;

      MatrixBinaryOuterIterator() {}

      explicit MatrixBinaryOuterIterator(iterator1_type const& i1, iterator2_type const& i2,
				    functor_type f = functor_type())
	 : i1_(i1), i2_(i2), f_(f) { }

      MatrixBinaryOuterIterator& operator++() 
	 { ++i1_; ++i2_; return *this; }

      MatrixBinaryOuterIterator& operator++(int) 
	 { return MatrixBinaryOuterIterator(i1_++, i2_++, f_); }

      MatrixBinaryOuterIterator& operator--() 
	 { --i1_; --i2_; return *this; }

      MatrixBinaryOuterIterator& operator--(int) 
	 { return MatrixBinaryOuterIterator(i1_--, i2_--, f_); }

      MatrixBinaryOuterIterator& operator+=(difference_type n)
	 { i1_ += n; i2_ += n; return *this; }

      MatrixBinaryOuterIterator& operator-=(difference_type n)
	 { i1_ -= n; i2_ -= n; return *this; }

      size_type index() const 
      { DEBUG_PRECONDITION(bool(i1_));
        DEBUG_CHECK_EQUAL(i1_.index(), i2_.index()); 
        return i1_.index(); }

      reference operator*() const { return transform(*i1_, *i2_, f_); }

      pointer operator->() const { return pointer(&transform(*i1_, *i2_, f_)); }

      operator bool() const { return i1_; }

      iterator1_type const& iterator1() const { return i1_; }
      iterator2_type const& iterator2() const { return i2_; }

      functor_type const& functor() const { return f_; }

      iterator iterate() const 
      { using LinearAlgebra::iterate; return iterator(iterate(i1_), iterate(i2_), f_); }

      // derived concepts
      size_type size() const { return i1_.size(); }
      reference operator[](difference_type n) const { return transform(i1_[n], i2_[n], f_); }

   private:
      iterator1_type i1_;
      iterator2_type i2_;
      functor_type f_;
};

// iterators

template <typename Iter1, typename Iter2, typename Func>
struct Iterate<MatrixBinaryOuterIterator<Iter1, Iter2, Func> >
{
   typedef typename MatrixBinaryOuterIterator<Iter1, Iter2, Func>::iterator result_type;
   typedef MatrixBinaryOuterIterator<Iter1, Iter2, Func> const& argument_type;
   result_type operator()(argument_type x) const 
   { return x.iterate(); }
};

} // namespace LinearAlgebra

#endif
