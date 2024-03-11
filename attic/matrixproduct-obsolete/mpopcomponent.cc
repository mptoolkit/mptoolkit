// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/matrixproduct-obsolete/mpopcomponent.cc
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
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

template <typename T>
std::ostream& operator<<(std::ostream& out, MPOperatorComponent<T> const& Op)
{
  out << "Site basis:\n" << Op.SiteBasis() << "Basis1:\n" << Op.Basis1()
      << "Basis2:\n" << Op.Basis2() << '\n';
  for (typename MPOperatorComponent<T>::const_iterator I = Op.begin(); I != Op.end(); ++I)
  {
     out << "Component transforms as " << I->first << '\n' << I->second << '\n';
  }
  return out;
}

template <typename T>
PStream::opstream& operator<<(PStream::opstream& out, MPOperatorComponent<T> const& Op)
{
   return out << Op.SBasis_ << Op.Basis1_ << Op.Basis2_ << Op.Data_;
}

template <typename T>
PStream::ipstream& operator>>(PStream::ipstream& in, MPOperatorComponent<T>& Op)
{
   return in >> Op.SBasis_ >> Op.Basis1_ >> Op.Basis2_ >> Op.Data_;
}

template <typename T>
MPOperatorComponent<T>& MPOperatorComponent<T>::operator*=(double x)
{
   for (typename MPOperatorComponent<T>::iterator I = this->begin(); I != this->end(); ++I)
      I->second *= x;

   return *this;
}

template <typename T>
MPOperatorComponent<T>& MPOperatorComponent<T>::operator*=(std::complex<double> x)
{
   for (typename MPOperatorComponent<T>::iterator I = this->begin(); I != this->end(); ++I)
      I->second *= x;

   return *this;
}

template <typename T>
void MPOperatorComponent<T>::check_structure() const
{
   typedef typename MPOperatorComponent<T>::const_iterator const_iterator;
   typedef typename MPOperatorComponent<T>::mapped_type mapped_type;
   for (typename MPOperatorComponent<T>::const_iterator I = this->begin(); I != this->end(); ++I)
   {
      CHECK_EQUAL(I->first, I->second.TransformsAs());
      CHECK_EQUAL(I->second.Basis1(), this->SiteBasis());
      CHECK_EQUAL(I->second.Basis2(), this->SiteBasis());
      I->second.check_structure();

      for (typename LinearAlgebra::const_iterator<mapped_type>::type X = iterate(I->second); X; ++X)
      {
         for (typename LinearAlgebra::const_inner_iterator<mapped_type>::type Y = iterate(X); Y; ++Y)
         {
            CHECK_EQUAL(I->first, Y->TransformsAs());
            CHECK_EQUAL(Y->Basis1(), this->Basis1());
            CHECK_EQUAL(Y->Basis2(), this->Basis2());

            Y->check_structure();
         }
      }
   }
}

template <typename T>
void MPOperatorComponent<T>::set_operator(int i, int j, SimpleOperator const& x)
{
   DEBUG_CHECK_EQUAL(SBasis_, x.Basis1());
   DEBUG_CHECK(is_transform_target(Basis2_[j], x.TransformsAs(), Basis1_[i]))
      (Basis2_[j])(x.TransformsAs())(Basis1_[i]);
   for (typename LinearAlgebra::iterator<SimpleOperator const>::type a = iterate(x); a; ++a)
   {
      for (typename LinearAlgebra::inner_iterator<SimpleOperator const>::type b = iterate(a); b; ++b)
      {
         if (this->operator[](x.TransformsAs())(b.index1(), b.index2()).is_null())
         {
            Data_[x.TransformsAs()](b.index1(), b.index2())
               = SimpleOperator(Basis1_, Basis2_, x.TransformsAs());
         }
         Data_[x.TransformsAs()](b.index1(), b.index2())(i,j) = *b;
      }
   }
}

namespace LinearAlgebra
{

template <typename T>
struct ScalarProdABh
{
   typedef T result_type;
   typedef T const& first_argument_type;
   typedef T const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      return scalar_prod(x, herm(y));
   }
};

template <typename T>
struct ScalarProdAhB
{
   typedef T result_type;
   typedef T const& first_argument_type;
   typedef T const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      return scalar_prod(herm(x), y);
   }
};

template <typename T>
typename ScalarProd<MPOperatorComponent<T>,
                    HermitianProxy<MPOperatorComponent<T> > >::result_type
ScalarProd<MPOperatorComponent<T>, HermitianProxy<MPOperatorComponent<T> > >::
operator()(MPOperatorComponent<T> const& A,
           HermitianProxy<MPOperatorComponent<T> > const& B) const
{
   typedef typename MPOperatorComponent<T>::OperatorType OperatorType;

   DEBUG_PRECONDITION_EQUAL(A.SiteBasis(), B.base().SiteBasis());
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), B.base().Basis2());

   //   A.check_structure();
   //B.base().check_structure();

   QuantumNumber Ident(A.GetSymmetryList());
   OperatorType Result(A.Basis1(), B.base().Basis1(), Ident);

   for (typename MPOperatorComponent<T>::const_iterator I = A.begin(); I != A.end(); ++I)
   {
      typename MPOperatorComponent<T>::const_iterator J = B.base().find(I->first);
      if (J != B.base().end())
         Result += parallel_prod(I->second, J->second, ScalarProdABh<OperatorType>());
   }
   return Result;
}

template <typename T>
typename ScalarProd<HermitianProxy<MPOperatorComponent<T> >,
                    MPOperatorComponent<T> >::result_type
ScalarProd<HermitianProxy<MPOperatorComponent<T> >, MPOperatorComponent<T> >::
operator()(HermitianProxy<MPOperatorComponent<T> > const& A,
           MPOperatorComponent<T> const& B) const
{
   typedef typename MPOperatorComponent<T>::OperatorType OperatorType;

   DEBUG_PRECONDITION_EQUAL(A.base().SiteBasis(), B.SiteBasis());
   DEBUG_PRECONDITION_EQUAL(A.base().Basis1(), B.Basis1());

   QuantumNumber Ident(B.GetSymmetryList());
   OperatorType Result(A.base().Basis2(), B.Basis2(), Ident);

   for (typename MPOperatorComponent<T>::const_iterator I = B.begin(); I != B.end(); ++I)
   {
      typename MPOperatorComponent<T>::const_iterator J = A.base().find(I->first);
      if (J != A.base().end())
         Result += parallel_prod(J->second, I->second, ScalarProdAhB<OperatorType>());
   }
   return Result;
}

} // namespace LinearAlgebra
