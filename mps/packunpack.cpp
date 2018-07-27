// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mps/packunpack.cpp
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

#include "packunpack.h"
#include <cstring>
#include "blas/vector_view.h"

// helper function.  Returns the new iterator pointing to one-past the end

template <typename T>
T*
pack_to(double Factor, blas::Matrix<T> const& m, T* Iter)
{
   // pack in column-major format
   for (int i = 0; i < m.cols(); ++i)
   {
      blas::make_vector_view<T>(m.rows(), 1, Iter+i*m.leading_dimension())
	 = Factor * m.column(i);
   }
   return Iter+m.rows()*m.cols();
}

template <typename T>
T const*
unpack_from(double Factor, blas::Matrix<T>& m, T const* Iter)
{
   // pack in column-major format
   for (int i = 0; i < m.cols(); ++i)
   {
      m.column(i) = Factor *
         blas::make_const_vector_view<T>(m.rows(), 1, Iter+i*m.leading_dimension());
   }
   return Iter+m.rows()*m.cols();
}

// determine the array size of the packed form of a MatrixOperator
int
pack_size(MatrixOperator const& m)
{
   int Result = 0;
   for (int i = 0; i < m.Basis1().size(); ++i)
   {
      for (int j = 0; j < m.Basis2().size(); ++j)
      {
         if (is_transform_target(m.Basis2()[j], m.TransformsAs(), m.Basis1()[i]))
         {
            Result += m.Basis1().dim(i) * m.Basis2().dim(j);
         }
      }
   }
   return Result;
}

complex*
pack_to(MatrixOperator const& m, complex* Iter)
{
   for (int i = 0; i < m.Basis1().size(); ++i)
   {
      for (int j = 0; j < m.Basis2().size(); ++j)
      {
         if (is_transform_target(m.Basis2()[j], m.TransformsAs(), m.Basis1()[i]))
         {
            int Size = m.Basis1().dim(i) * m.Basis2().dim(j);
            auto J = m.row(i).find(j);
            if (J != m.row(i).end())
            {
               // normalization correction is sqrt(degree(B1[I->r]))
               pack_to(std::sqrt(double(degree(m.Basis1()[i]))), get_wait(J.value()), Iter);
            }
            else
            {
               std::memset(Iter, 0, Size * sizeof(complex));
            }
            Iter += Size;
         }
      }
   }
   return Iter;
}

blas::Vector<complex>
pack(MatrixOperator const& m)
{
   blas::Vector<complex> Result(pack_size(m));
   pack_to(m, Result.storage());
   return Result;
}

// unpack a linear vector into a MatrixOperator.  The operator
// must have the appropriate basis and quantum numbers set
complex const*
unpack_from(MatrixOperator& m, complex const* Iter)
{
   for (int i = 0; i < m.Basis1().size(); ++i)
   {
      for (int j = 0; j < m.Basis2().size(); ++j)
      {
         if (is_transform_target(m.Basis2()[j], m.TransformsAs(), m.Basis1()[i]))
         {
            double Factor = 1.0 / std::sqrt(double(degree(m.Basis1()[i])));
            blas::Matrix<complex> Temp(m.Basis1().dim(i), m.Basis2().dim(j));
            Iter = unpack_from(Factor, Temp, Iter);
            blas::Matrix<complex, MatrixOperator::tag_type>
               DeviceTemp(Temp.rows(), Temp.cols());
            set_wait(DeviceTemp, Temp);
            m.set(i,j,DeviceTemp);
         }
      }
   }
   return Iter;
}

void unpack(MatrixOperator& m, blas::Vector<complex> const& v)
{
   CHECK_EQUAL(v.size(), pack_size(m));
   unpack_from(m, v.storage());
}

int
pack_size(StateComponent const& m)
{
   int Result = 0;
   for (auto const& mm : m)
   {
      Result += pack_size(mm);
   }
   return Result;
}

blas::Vector<complex>
pack(StateComponent const& m)
{
   blas::Vector<complex> Result(pack_size(m));
   complex* Iter = Result.storage();
   for (auto const& mm : m)
   {
      Iter = pack_to(mm, Iter);
   }
   CHECK_EQUAL(Iter-Result.storage(), Result.size());
   return Result;
}

void unpack(StateComponent& m, blas::Vector<complex> const& v)
{
   complex const* Iter = v.storage();
   for (auto& mm : m)
   {
      Iter = unpack_from(mm, Iter);
   }
   CHECK_EQUAL(Iter-v.storage(), v.size());
}
