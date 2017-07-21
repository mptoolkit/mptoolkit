// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/operatoratsite.h
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

#if !defined(OPERATORATSITE_H_HSIUYHUY4878Y78Y7)
#define OPERATORATSITE_H_HSIUYHUY4878Y78Y7

#include "mpoperatorlist.h"
#include <boost/mpl/if.hpp>
#include <boost/type_traits.hpp>

template <typename OListType, typename T1 = void, typename T2 = void, typename T3 = void>
class OperatorAtSite;

// all others are now depreciated, should use this version
template <typename OListType>
class OperatorAtSite<OListType, void>
{
   public:
      typedef OListType OperatorListType;
      typedef typename OperatorListType::OperatorType OperatorType;
      typedef OperatorType value_type;

      typedef typename boost::mpl::if_<boost::is_const<OperatorListType>,
                                       OperatorType,
                                       OperatorType&>::type reference;
      typedef OperatorType const_reference;

      OperatorAtSite(OperatorListType& OList, std::string const& OName)
         : OList_(OList), OName_(OName) {}

      template <typename T>
      std::string FullName(T const& x) const
      {
         return OName_ + '(' + boost::lexical_cast<std::string>(x) + ')';
      }

      template <typename T1, typename T2>
      std::string FullName(T1 const& x1, T2 const& x2) const
      {
         return OName_ + '(' + boost::lexical_cast<std::string>(x1) + ','
            + boost::lexical_cast<std::string>(x2) + ')';
      }

      template <typename T1, typename T2, typename T3>
      std::string FullName(T1 const& x1, T2 const& x2, T3 const& x3) const
      {
         return OName_ + '(' + boost::lexical_cast<std::string>(x1) + ','
            + boost::lexical_cast<std::string>(x2) + ','
            + boost::lexical_cast<std::string>(x3) + ')';
      }

      // 1-argument versions
      template <typename T>
      bool HasOperator(T const& x) const
      {
         return OList_.HasOperator(this->FullName(x));
      }

      template <typename T>
      const_reference operator()(T const& x) const
      {
         return OList_[this->FullName(x)];
      }

      template <typename T>
      reference operator()(T const& x)
      {
         return OList_[this->FullName(x)];
      }

      template <typename T>
      const_reference operator[](T const& x) const
      {
         return OList_[this->FullName(x)];
      }

      template <typename T>
      reference operator[](T const& x)
      {
         return OList_[this->FullName(x)];
      }

      // 2-argument versions
      template <typename T1, typename T2>
      bool HasOperator(T1 const& x1, T2 const& x2) const
      {
         return OList_.HasOperator(this->FullName(x1, x2));
      }

      template <typename T1, typename T2>
      const_reference operator()(T1 const& x1, T2 const& x2) const
      {
         return OList_[this->FullName(x1, x2)];
      }

      template <typename T1, typename T2>
      reference operator()(T1 const& x1, T2 const& x2)
      {
         return OList_[this->FullName(x1, x2)];
      }

      // 3-argument versions
      template <typename T1, typename T2, typename T3>
      bool HasOperator(T1 const& x1, T2 const& x2, T3 const& x3) const
      {
         return OList_.HasOperator(this->FullName(x1, x2, x3));
      }

      template <typename T1, typename T2, typename T3>
      const_reference operator()(T1 const& x1, T2 const& x2, T3 const& x3) const
      {
         return OList_[this->FullName(x1, x2, x3)];
      }

      template <typename T1, typename T2, typename T3>
      reference operator()(T1 const& x1, T2 const& x2, T3 const& x3)
      {
         return OList_[this->FullName(x1, x2, x3)];
      }

   private:
      OperatorListType& OList_;
      std::string OName_;
};


template <typename OListType, typename T>
class OperatorAtSite<OListType, T, void>
{
   public:
      typedef OListType OperatorListType;
      typedef typename OperatorListType::OperatorType OperatorType;
      typedef OperatorType value_type;

      typedef typename boost::mpl::if_<boost::is_const<OperatorListType>,
                                       OperatorType,
                                       OperatorType&>::type reference;
      typedef OperatorType const_reference;

      OperatorAtSite(OperatorListType& OList, std::string const& OName)
         : OList_(OList), OName_(OName) {}

      std::string FullName(T const& x) const
      {
         return OName_ + '(' + boost::lexical_cast<std::string>(x) + ')';
      }

      bool HasOperator(T const& x) const
      {
         return OList_.HasOperator(this->FullName(x));
      }

      const_reference operator()(T const& x) const
      {
         return OList_[this->FullName(x)];
      }

      reference operator()(T const& x)
      {
         return OList_[this->FullName(x)];
      }

      const_reference operator[](T const& x) const
      {
         return OList_[this->FullName(x)];
      }

      reference operator[](T const& x)
      {
         return OList_[this->FullName(x)];
      }

   private:
      OperatorListType& OList_;
      std::string OName_;
};

template <typename OListType, typename T1, typename T2>
class OperatorAtSite<OListType, T1, T2>
{
   public:
      typedef OListType OperatorListType;
      typedef typename OperatorListType::OperatorType OperatorType;
      typedef OperatorType value_type;

      typedef typename boost::mpl::if_<boost::is_const<OperatorListType>,
                                       OperatorType,
                                       OperatorType&>::type reference;
      typedef OperatorType const_reference;

      OperatorAtSite(OperatorListType& OList, std::string const& OName)
         : OList_(OList), OName_(OName) {}

      std::string FullName(T1 const& x1, T2 const& x2) const
      {
         return OName_ + '(' + boost::lexical_cast<std::string>(x1) + ','
            + boost::lexical_cast<std::string>(x2) + ')';
      }

      bool HasOperator(T1 const& x1, T2 const& x2) const
      {
         return OList_.HasOperator(this->FullName(x1, x2));
      }

      const_reference operator()(T1 const& x1, T2 const& x2) const
      {
         return OList_[this->FullName(x1, x2)];
      }

      reference operator()(T1 const& x1, T2 const& x2)
      {
         return OList_[this->FullName(x1, x2)];
      }

   private:
      OperatorListType& OList_;
      std::string OName_;
};

template <typename OListType, typename T1, typename T2, typename T3>
class OperatorAtSite //<OListType, T1, T2, T3>
{
   public:
      typedef OListType OperatorListType;
      typedef typename OperatorListType::OperatorType OperatorType;
      typedef OperatorType value_type;

      typedef typename boost::mpl::if_<boost::is_const<OperatorListType>,
                                       OperatorType,
                                       OperatorType&>::type reference;
      typedef OperatorType const_reference;

      OperatorAtSite(OperatorListType& OList, std::string const& OName)
         : OList_(OList), OName_(OName) {}

   std::string FullName(T1 const& x1, T2 const& x2, T3 const& x3) const
      {
         return OName_ + '(' + boost::lexical_cast<std::string>(x1) + ','
            + boost::lexical_cast<std::string>(x2) + ','
            + boost::lexical_cast<std::string>(x3) + ')';
      }

   bool HasOperator(T1 const& x1, T2 const& x2, T3 const& x3) const
      {
         return OList_.HasOperator(this->FullName(x1, x2, x3));
      }

   const_reference operator()(T1 const& x1, T2 const& x2, T3 const& x3) const
      {
         return OList_[this->FullName(x1, x2, x3)];
      }

   reference operator()(T1 const& x1, T2 const& x2, T3 const& x3)
      {
         return OList_[this->FullName(x1, x2, x3)];
      }

   private:
      OperatorListType& OList_;
      std::string OName_;
};

#endif
