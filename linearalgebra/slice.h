/* -*- C++ -*- $Id$

  classes Range and Slice to perform basic indirection

  Created March 2004 Ian McCulloch

  Range and Slice are themselves Vector types, with a value_type of LinearAlgebra::size_type.

  The constructor Range(size_type first, size_type last) constructs a Range object
  which acts like a sequence vector; the elements of the vector being
  [first, first+1, ..., last).  The value last is not itself a member of the container.

  The constructor Slice(size_type Start, size_type Size, difference_type Stride)
  constructs a Slice object which acts like a vector with elements
  [Start, Start + Stride, Start + Stride*2, ..., Start + Stride*Size).  The
  end element (Start + Stride*Size) is not a member of the container.
*/

#if !defined(SLICE_H_IUIU43895U4398UP98UPWP)
#define SLICE_H_IUIU43895U4398UP98UPWP

#include "vectoroperationsbase.h"

namespace LinearAlgebra
{

//
// Range
//
// A Range represents some segment of another vector.
// It is itself a dense vector, as well as a unary function.
//

class Range;

class RangeIterator
{
   public:
      typedef size_type value_type;
      typedef size_type reference;
      typedef size_type const* pointer;
      typedef vector_iterator_dense category;

      RangeIterator() {}

      RangeIterator(size_type First, size_type Last)
         : First_(First), Last_(Last), n_(First) {}

      RangeIterator(size_type First, size_type Last, size_type n)
         : First_(First), Last_(Last), n_(n) {}

      RangeIterator& operator++() { ++n_; return *this; }
      RangeIterator operator++(int) { return RangeIterator(First_, Last_, n_++); }
      RangeIterator& operator+=(difference_type i) { n_ += i; return *this; }

      size_type operator[](difference_type i) const { return difference_type(n_)+i; }

      size_type operator*() const { return n_; }

      // TODO: operator->

      size_type index() const { return n_ - First_; }

      operator bool() const { return n_ != Last_; }

      size_type first() const { return First_; }
      size_type last() const { return Last_; }
  
   private:
      size_type First_, Last_, n_;
};

class Range
{
   public:
      typedef size_type value_type;
      typedef RangeIterator const_iterator;
      typedef size_type const_reference;

      // a range acts as a proxy
      typedef boost::mpl::true_ const_proxy;

      // force the assumption that no aliasing of a Range can occur
      typedef boost::mpl::true_ noalias;

      Range() : First_(0), Last_(0) {}

      Range(size_type First, size_type Last)
	: First_(First), Last_(Last) { DEBUG_PRECONDITION(First <= Last); }

      Range& operator=(Range const& r) { First_ = r.First_; Last_ = r.Last_; return *this; }

      size_type first() const { return First_; }
      size_type last() const { return Last_; }

      size_type size() const { return Last_-First_; }

      size_type operator()(size_type n) const { return n + First_; }
      size_type operator[](size_type n) const { return n + First_; }

      const_iterator iterate() const { return RangeIterator(First_, Last_); }

      difference_type stride() const { return 1; }

      bool operator==(Range const& r) const { return First_ == r.First_ && Last_ == r.Last_; }
      bool operator!=(Range const& r) const { return First_ != r.First_ || Last_ != r.Last_; }

   private:
      size_type First_, Last_;
};

inline
Range range(size_type first, size_type last)
{
   return Range(first, last);
}

// Blitz++ style arithmetic on ranges

inline
Range
operator+(Range const& r, difference_type x)
{
   return Range(difference_type(r.first())+x, difference_type(r.last())+x);
}

inline
Range
operator-(Range const& r, difference_type x)
{
   return Range(difference_type(r.first())-x, difference_type(r.last())-x);
}

// FIXME: these should be implemented as specializations of Multiplication

inline
Range
operator*(Range const& r, size_type n)
{
   return Range(r.first()*n, r.last()*n);
}

inline
Range
operator*(size_type n, Range const& r)
{
   return Range(n*r.first(), n*r.last());
}

template <>
struct Stride<Range, DENSE_VECTOR(size_type, void)>
{
   typedef difference_type result_type;
   typedef Range const& argument_type;
   result_type operator()(argument_type x) const
   {
      return 1;
   }
};

// interface

template <>
struct interface<Range>
{
   typedef DENSE_VECTOR(size_type, void) type;
   typedef size_type value_type;
};

template <>
struct StreamInsert<Range, DENSE_VECTOR(size_type, Range)>
{
   typedef std::ostream& result_type;
   typedef std::ostream& first_argument_type;
   typedef Range second_argument_type;
   result_type operator()(std::ostream& out, Range const& r) const
   {
      return out << '[' << r.first() << ", " << r.last() << '[';
   }
};


//
// Slice
//
// A Slice represents a section view of a container, slice(i) = start + stride * i.
// It is itself a vector, implementing the VectorConstStride interface.
//

class SliceIterator
{
   public:
      typedef size_type value_type;
      typedef size_type reference;
      typedef size_type const* pointer;
      typedef vector_iterator_dense category;

      SliceIterator() {}
      SliceIterator(size_type Start, size_type Size, difference_type Stride, size_type n = 0) 
	 : Start_(Start), Size_(Size), Stride_(Stride), n_(n) {}

      SliceIterator& operator++() { ++n_; return *this; }
      SliceIterator operator++(int) { return SliceIterator(Start_, Size_, Stride_, n_++); }

      SliceIterator& operator--() { --n_; return *this; }
      SliceIterator operator--(int) { return SliceIterator(Start_, Size_, Stride_, n_--); }

      SliceIterator& operator+=(difference_type i) { n_ += i; return *this; }
      SliceIterator& operator-=(difference_type i) { n_ -= i; return *this; }

      size_type operator*() const 
      { return difference_type(Start_) + difference_type(n_) * Stride_; }

      size_type operator[](difference_type n) const
      {
         return difference_type(Start_) + (difference_type(n_) + n) * Stride_;
      }

      // TODO: operator->

      size_type index() const { return n_; }

      operator bool() const { return n_ != Size_; }

      size_type size() const { return Size_; }
      size_type last() const { return difference_type(Start_) + difference_type(Size_) * Stride_; }

      difference_type stride() const { return Stride_; }
  
   private:
      size_type Start_;
      size_type Size_;
      difference_type Stride_;
      size_type n_;
};

class Slice
{
   public:
      typedef size_type value_type;
      typedef SliceIterator const_iterator;
      typedef size_type const_reference;

      // a slice acts as a proxy
      typedef boost::mpl::true_ const_proxy;

      // force the assumption that no aliasing of a Slice can occur
      typedef boost::mpl::true_ noalias;

      Slice() : Start_(0), Size_(0), Stride_(0) {}
      Slice(size_type Start, size_type Size, difference_type Stride)
        : Start_(Start), Size_(Size), Stride_(Stride) {}

      size_type start() const { return Start_; }
      size_type size() const { return Size_; }
      difference_type stride() const { return Stride_; }

      size_type operator()(size_type n) const 
      { return difference_type(Start_) + difference_type(n)*Stride_; }

      size_type operator[](size_type n) const 
      { return difference_type(Start_) + difference_type(n)*Stride_; }

      bool operator==(Slice const& s) const 
      { return Start_ == s.Start_ && Size_ == s.Size_ && Stride_ == s.Stride_; }

      bool operator!=(Slice const& s) const 
      { return Start_ != s.Start_ || Size_ != s.Size_ || Stride_ != s.Stride_; }

      const_iterator iterate() const 
      { return const_iterator(Start_, Size_, Stride_); }

   private:
      size_type Start_, Size_;
      difference_type Stride_;
};

template <>
struct Stride<Slice, DENSE_VECTOR(size_type, void)>
{
   typedef difference_type result_type;
   typedef Slice const& argument_type;
   result_type operator()(argument_type x) const
   {
      return x.stride();
   }
};

inline
Slice slice(size_type start, size_type size, difference_type stride)
{
   return Slice(start, size, stride);
}

inline
Slice range(size_type first, size_type last, size_type increment)
{
   return Slice(first, (last-first)/increment, increment);
}

// Blitz++ style arithmetic on slices

inline
Slice
operator+(Slice const& s, difference_type x)
{
   return Slice(difference_type(s.start())+x, s.size(), s.stride());
}

// FIXME: these should be implemented as specializations of Multiplication

inline
Slice
operator-(Slice const& s, difference_type x)
{
   return Slice(difference_type(s.start())-x, s.size(), s.stride());
}

inline
Slice
operator*(Slice const& s, size_type n)
{
   return Slice(s.start()*n, s.size(), s.stride()*n);
}

inline
Slice
operator*(size_type n, Slice const& s)
{
   return Slice(n*s.start(), s.size(), n*s.stride());
}

// interface

template <>
struct interface<Slice>
{
   typedef DENSE_VECTOR(size_type, void) type;
   typedef size_type value_type;
};

template <>
struct StreamInsert<Slice, DENSE_VECTOR(size_type, Slice)>
{
   typedef std::ostream& result_type;
   typedef std::ostream& first_argument_type;
   typedef Slice second_argument_type;
   result_type operator()(std::ostream& out, Slice const& s) const
   {
      return out << "slice(start = " << s.start() 
		 << ", size = " << s.size() 
		 << ", stride = " << s.stride() << ')';
   }
};

// functional composition of slices.
// Returns a slice S such that S(i) = s1(s2(i))
// This is taking the subset s2 of s1.
inline
Slice compose(Slice const& s1, Slice const& s2)
{
   PRECONDITION(s2(s2.size()) < s1.size() + s2.stride())(s1)(s2);
   return Slice(s1(s2.start()), s2.size(), s1.stride() * s2.stride()); 
}

// functional composition of a slice and a range.
// Returns a slice S such that S(i) = s(r(i))
inline
Slice compose(Slice const& s, Range const& r)
{
   PRECONDITION(r.last() <= s(s.size()))(s)(r);
   return Slice(s(r.first()), r.size(), s.stride()); 
}

// functional composition of a range and a slice.
// Returns a slice S such that S(i) = r(s(i))
inline
Slice compose(Range const& r, Slice const& s) 
{
   PRECONDITION(s(s.size()) < r.last() + s.stride())(r)(s);
   return Slice(r.first() + s.start(), s.size(), s.stride()); 
}

// functional composition of ranges.
// Returns a range R such that R(i) = r1(r2(i))
inline
Range compose(Range const& r1, Range const& r2)
{
   PRECONDITION(r1.first() + r2.last() <= r1.last())(r1)(r2);
   return Range(r1.first() + r2.first(), r1.first() + r2.last());
}

} // namespace LinearAlgebra

#endif
