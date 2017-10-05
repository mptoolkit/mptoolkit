// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/vectortypes.h
//
// Copyright (C) 2000-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
  expression templates fast vector library

  Created 2000-08 Ian McCulloch

  This library uses a non-virtual heirachy, based around the curiously recursive template
  pattern.
  That is, the derived class itself is a template parameter of the base class.
  The classes that have an explicit template parameter Derived should never be constructed
  directly.  These classes should have no public constructors if possible.
  The heirachy rooted at VectorExpression<Scalar const, Derived>.  This is
  a generic R-value, functions that accept a minimal interface non-mutable vector
  should use a parameter of this form.

  L-value objects have the root VectorExpression<Scalar, Derived>, which inherits from
  VectorExpression<Scalar const, Derived>.  This is a generic L-value object, with assignment,
  multiply by scalar etc etc.

  VectorExpression does not declare any constructors, dtor, or data members.  This is left
  totally up to
  the derived classes.  VectorExpression L-value assumes reference semantic copy and
  copy-assignment.
  For example,
  Vector<double> MyVector(20, 0);
  MyVector.sub_vector(5, 15) = Vector<double>(10, 1);

  Here, sub_vector() returns a VectorSlice object which is derived from
  VectorExpression<double>,
  which encapsulates a reference to a sub-vector of MyVector.
  The assignment operator implements reference counted assignment, and assigns 1.0 to the
  relevant portions
  of MyVector.

  Vector<Scalar const> objects themselves cannot appear on the left hand side of an
  expression, thus
  it is purely up to the derived classes how to implement assignment (if it is implemented at
  all).

  Vector<Scalar> inherits from VectorExpression<Scalar const> and provides deep copy
  assignment semantics.

  Derived classes:
  +  VectorSlice
        Implements a slice of fixed stride, can appear as an L-value.
  +  VectorRef
        Implements a reference to a vector, can appear as an L-value.
  +  Vector
        Basic vector class, deep copy assignment semantics.

  Examples:

  Vector<double> MyVector(10, 0);     // declare a vector of length 10, initialized to zero
  VectorRef<double> MyRef(MyVector);  // MyRef acts as a reference to MyVector

  MyVector[0] = 1;                    // changes the first element of MyVector and MyRef
  MyRef[1] = 2;                       // changes the second element of MyVector and MyRef

  MyRef = Vector<double>(10, 3);      // sets every element of MyVector and MyRef to 3
  // MyRef =  Vector<double>(11, 3); **ILLEGAL** trying to assign a length 11 vector to a
  //                                 length 10 vector

  MyVector = Vector<double>(11, 4);   // OK, MyVector is now length 11, MyRef is *unchanged*
  MyRef[3] = 5;                       // changes the 4th element of MyRef, MyVector is unchanged

  The reference counted semantics of VectorRef and Vector might seem unusual, but they're not
  really.
  Essentially, VectorRef is equivalent to a reference counted iterator to begin() of an
  std::vector.
  Iterator invalidation corresponds to the reference becoming detached from the original
  object,
  but here, because of the reference counting, the original data does not become invalid.

  TODO: add a lot more classes, VectorAddExpression, VectorScaleExpression,
  VectorGeneralSlice,
  and some sequence vectors - ConstSequence, ArithmeticSequence, GeometricSequence etc etc.

  Finish off the VectorExpression<Scalar> members, should have fill(), clear() etc etc.

  Finish off the Vector and VectorRef members - should have sub_vector(), slice() etc etc.

  Document some more about the iterator requirements - we are a lot more restrictive than STL.

        Document the VectorTraits class.

        copy and assignment of the base classes is a problem - really we want the minimum to
        be defined
        and the rest to be private and not defined.  What can we get away with?  We need to
        have public copy for VectorExpression<Scalar> at minimum.  Not really needed for
        VectorExpression<Scalar const>
        we can pass it by const reference.  Or do we need it for something?
        If we do make the copy and copy-assignment private, the base classes need to
        explicitly define
        them so we don't try to call the base class versions.

        Vector and VectorRef should interit from VectorSlice.

        Vector doesn't need to inherit from VectorExpression<Scalar const>, it is ok for it to
        inherit from VectorExpression<Scalar>.
*/

#if !defined(VECTORTYPES_H_DHJFKHURWEIRY4572893475489YRUI34897)
#define VECTORTYPES_H_DHJFKHURWEIRY4572893475489YRUI34897

#include "config.h"
#include "common/trace.h"
#include "vectorinterface.h"
#include "iteratortypes.h"
#include "vectorderived.h"
#include "datablockreference.h"

#include <boost/utility/enable_if.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/type_traits.hpp>

#if defined(USE_PSTREAM)
#include "pstream/pstream.h"
#endif

namespace LinearAlgebra
{

//
// Vector
//
// vector with copy-on-write semantics.  This is the 'default' vector class.
//

//template <typename Scalar>
//class Vector;

//
// VectorConstRef
//

template <class Scalar>
struct VectorTraits<VectorConstRef<Scalar, ConcreteClass> >
{
   typedef Scalar const*                                const_iterator;
   typedef VectorConstRef<Scalar>                       container_const_reference;
   typedef Scalar const&                                const_result_type;
   typedef VectorConstIndex<Scalar>                     const_index_type;
   typedef VectorConstSlice<Scalar>                     const_slice_type;
   typedef VectorConstRange<Scalar>                     const_range_type;
};

template <typename Scalar>
class VectorConstRef<Scalar, ConcreteClass>
  : public VectorSlice<Scalar, VectorRef<Scalar, ConcreteClass>, VectorRef<Scalar, ConcreteClass> >
{
   public:
      typedef VectorSlice<Scalar, VectorRef<Scalar, ConcreteClass>,
                          VectorRef<Scalar, ConcreteClass> > base_class;
      typedef Scalar*       iterator;
      typedef Scalar const* const_iterator;

      using base_class::operator[];

  //      VectorConstRef(DataBlock<Scalar const> const& Data_, size_type Start, size_type Size)
  //    : base_class(), Data(Data_) {}

      VectorConstRef(VectorConstRef const& V)
        : base_class(), Block(V.Block) {}

      template <typename S>
      VectorConstRef(VectorConstRef<S> const& V)
        : base_class(), Block(V.Block) {}

      template <typename S>
      VectorConstRef(VectorRef<S> const& V)
        : base_class(), Block(V.Block) {}

      size_type size() const { return Block.size(); }

      const_iterator begin() const { return Block.get(); }
      const_iterator end() const { return Block.get() + this->size(); }

      Scalar const& operator[](size_type n) const
         { DEBUG_RANGE_CHECK_OPEN(n,0,this->size()); return *(Block.get() + n); }

      Scalar const* data() const { return Block.get(); }

      // members defined in VectorConstStride
      difference_type stride() const { return 1; }
      size_type start() const { return 0; }

   private:
      VectorConstRef(); // not implemented

      BlockReference<Scalar const> Block;

   template <typename U, typename V> friend class VectorConstRef;
};

//
// VectorRef
//

template <class Scalar>
struct VectorTraits<VectorRef<Scalar, ConcreteClass> >
  : public VectorTraits<VectorConstRef<Scalar, ConcreteClass> >
{
   typedef VectorConstRef<Scalar, VectorRef<Scalar, ConcreteClass> > rvalue_type;
   typedef Scalar*                                                   iterator;
   typedef VectorRef<Scalar>                                         container_reference;
   typedef Scalar&                                                   result_type;
   typedef VectorIndex<Scalar>                                       index_type;
   typedef VectorSlice<Scalar>                                       slice_type;
   typedef VectorRange<Scalar>                                       range_type;
};

template <typename Scalar>
class VectorRef<Scalar, ConcreteClass>
  : public VectorSlice<Scalar, VectorRef<Scalar, ConcreteClass>, VectorRef<Scalar, ConcreteClass> >
{
   public:
      typedef VectorSlice<Scalar, VectorRef<Scalar, ConcreteClass>,
                          VectorRef<Scalar, ConcreteClass> >          base_class;
      typedef Scalar*       iterator;
      typedef Scalar const* const_iterator;

      using base_class::operator[];

      VectorRef(VectorRef const& V)
        : base_class(), Block(V.Block) {}

      template <class S2, class D2>
      VectorRef& operator=(VectorConstExpression<S2, D2> const& V)
      {
         this->assign(V);   // safe because assign() calls begin(), which calls cow()
         return *this;
      }

      VectorRef& operator=(VectorRef const& V)
      {
         this->assign(V);   // safe because assign() calls begin(), which calls cow()
         return *this;
      }

      size_type size() const { return Block.size(); }

      iterator begin() { this->cow(); return Block.get(); }
      iterator end()   { this->cow(); return Block.get() + this->size(); }

      const_iterator begin() const { return Block.get(); }
      const_iterator end() const { return Block.get() + this->size(); }

      Scalar const& operator[](size_type n) const
         { DEBUG_RANGE_CHECK_OPEN(n,0,this->size()); return *(Block.get() + n); }

      Scalar& operator[](size_type n)
         { DEBUG_RANGE_CHECK_OPEN(n,0,this->size());
         this->cow(); return *(Block.get() + n); }

      Scalar* data() { this->cow(); return Block.get(); }
      Scalar const* data() const { return Block.get(); }

      // members defined in VectorStride
      difference_type stride() const { return 1; }
      size_type start() const { return 0; }

   protected:
      BlockReference<Scalar> const& block() const { return Block; }
      BlockReference<Scalar>& block() { return Block; }  // danger: this does not call cow()

      VectorRef() {}

      explicit VectorRef(size_type Size_)
        : Block(Size_) {}

      VectorRef(size_type Size_, Scalar const& Fill)
        : Block(Size_, Fill) {}

     explicit VectorRef(BlockReference<Scalar> const& Block_) : Block(Block_) {}

     void resize(size_type NewSize, bool PreserveValues);

     void swap(VectorRef& Other) { std::swap(Block, Other.Block); }

   protected:
      void cow() { Block.cow(); }

      BlockReference<Scalar> Block;

   template <typename U, typename V> friend class VectorRef;
   template <typename U, typename V> friend class VectorConstRef;
};

//
// Vector
//
// vector with copy-on-write semantics.  This is the 'default' vector class.
//

template <typename Scalar>
class Vector;

template <class Scalar>
struct VectorTraits<Vector<Scalar> >
{
   typedef Scalar const*                                             const_iterator;
   typedef VectorConstRef<Scalar>                                    container_const_reference;
   typedef Scalar const&                                             const_result_type;
   typedef VectorIndex<Scalar, VectorConstRef<Scalar> >              const_index_type;
   typedef VectorSlice<Scalar, VectorConstRef<Scalar> >              const_slice_type;
   typedef VectorRef<Scalar>                                         const_range_type;

   typedef VectorConstRef<Scalar, VectorRef<Scalar> >                rvalue_type;
   typedef Scalar*                                                   iterator;
   typedef VectorRef<Scalar>                                         container_reference;
   typedef Scalar&                                                   result_type;
   typedef VectorIndex<Scalar, VectorRef<Scalar> >                   index_type;
   typedef VectorSlice<Scalar, VectorRef<Scalar> >                   slice_type;
   typedef VectorRef<Scalar, ConcreteClass>                          range_type;
};

template <class Scalar>
class Vector : public VectorRef<Scalar>
{
   public:
      typedef Scalar*       iterator;
      typedef Scalar const* const_iterator;

      Vector() {}

      explicit Vector(size_type Size_);

      Vector(size_type Size_, Scalar const& fill_);

      template <typename Iter>
      Vector(Iter first, Iter last,
             typename boost::enable_if<boost::mpl::not_<boost::is_arithmetic<Iter> > >::type* dummy = 0);

      Vector(Vector<Scalar> const& V)
        : VectorRef<Scalar>(V.block().copy()) {}

      Vector(VectorRef<Scalar> const& V)
        : VectorRef<Scalar>(V.block().copy()) {}

   //      Vector(VectorConstRef<Scalar> const& V)
   //   : VectorRef<Scalar>(V.block().copy()) {}

      template <class S2, class D2>
      Vector(VectorConstExpression<S2, D2> const& V);

      Vector<Scalar>& operator=(Vector<Scalar> const& V);

      template <class S2, class D2>
      Vector<Scalar>& operator=(VectorConstExpression<S2 const, D2> const& V);

      void resize(size_type NewSize, bool PreserveValues = false)
        { this->VectorRef<Scalar, ConcreteClass>::resize(NewSize, PreserveValues); }

      void swap(Vector& Other) { this->cow(); Other.cow(); this->VectorRef<Scalar>::swap(Other); }
};

#if defined(USE_PSTREAM)
template <int Format, typename T>
PStream::opstreambuf<Format>& operator<<(PStream::opstreambuf<Format>& out, Vector<T> const& Vec);

template <int Format, typename T>
PStream::ipstreambuf<Format>& operator>>(PStream::ipstreambuf<Format>& in, Vector<T>& Vec);
#endif

} // namespace LinearAlgebra

#include "vectortypes.cc"

#endif
