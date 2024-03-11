// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// tensor/tensor_eigen.cpp
//
// Copyright (C) 2004-2024 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#include "tensor_eigen.h"
#include "regularize.h"
#include "linearalgebra/eigen.h"
#include "linearalgebra/matrix_utility.h"

namespace Tensor
{

IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, VectorBasis>
DiagonalizeHermitianInPlace(IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, VectorBasis>& x)
{
   typedef IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, VectorBasis>
      TensorType;

   DEBUG_CHECK(is_scalar(x.TransformsAs()));
   DEBUG_CHECK_EQUAL(x.Basis1(), x.Basis2());
   DEBUG_CHECK(is_regular_basis(x.Basis1()))(x.Basis1());

   TensorType Result(x);
   for (iterator<TensorType>::type I = iterate(Result); I; ++I)
   {
      for (inner_iterator<TensorType>::type J = iterate(I); J; ++J)
      {
         DEBUG_CHECK_EQUAL(J.index1(), J.index2());
         LinearAlgebra::Vector<double> EVal = LinearAlgebra::DiagonalizeHermitian(*J);
         x(J.index1(), J.index2()) = LinearAlgebra::diagonal_matrix(EVal);
      }
   }

   return Result;
}


IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>
DiagonalizeHermitianInPlace(IrredTensor<LinearAlgebra::Matrix<std::complex<double>>,
                     VectorBasis, VectorBasis>& x)
{
   typedef IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>
      TensorType;

   DEBUG_CHECK(is_scalar(x.TransformsAs()));
   DEBUG_CHECK_EQUAL(x.Basis1(), x.Basis2());

   // handle the case of a non-regular basis
   if (!is_regular_basis(x.Basis1()))
   {
      DEBUG_WARNING("DiagonalizeHermitianInPlace with an irregular basis");
      Regularizer R(x.Basis1());
      auto M = RegularizeBasis12(R, x, R);
      auto Result = UnregularizeBasis12(R, DiagonalizeHermitianInPlace(M), R);
      M = UnregularizeBasis12(R, M, R);
      return Result;
   }

   TensorType Result(x);
   for (iterator<TensorType>::type I = iterate(Result); I; ++I)
   {
      for (inner_iterator<TensorType>::type J = iterate(I); J; ++J)
      {
         DEBUG_CHECK_EQUAL(J.index1(), J.index2());
         LinearAlgebra::Vector<double> EVal = LinearAlgebra::DiagonalizeHermitian(*J);
         x(J.index1(), J.index2()) = LinearAlgebra::diagonal_matrix(EVal);
      }
   }

   return Result;
}

std::tuple<IrredTensor<LinearAlgebra::DiagonalMatrix<double>, VectorBasis, VectorBasis, Tensor::DiagonalStructure>,
           IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, VectorBasis>>
DiagonalizeHermitian(IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, VectorBasis> x)
{
   typedef IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, VectorBasis>
      TensorType;

   DEBUG_CHECK(is_scalar(x.TransformsAs()));
   DEBUG_CHECK_EQUAL(x.Basis1(), x.Basis2());
   DEBUG_CHECK(is_regular_basis(x.Basis1()))(x.Basis1());

   IrredTensor<LinearAlgebra::DiagonalMatrix<double>, VectorBasis, VectorBasis, Tensor::DiagonalStructure> EValues(x.Basis1(), x.Basis2());

   for (iterator<TensorType>::type I = iterate(x); I; ++I)
   {
      for (inner_iterator<TensorType>::type J = iterate(I); J; ++J)
      {
         DEBUG_CHECK_EQUAL(J.index1(), J.index2());
         LinearAlgebra::Vector<double> EVal = LinearAlgebra::DiagonalizeHermitian(*J);
         EValues(J.index1(), J.index2()) = LinearAlgebra::diagonal_matrix(EVal);
      }
   }

   return std::make_tuple(EValues, x);
}

std::tuple<IrredTensor<LinearAlgebra::DiagonalMatrix<double>, VectorBasis, VectorBasis, Tensor::DiagonalStructure>,
           IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>>
DiagonalizeHermitian(IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis> x)
{
   DEBUG_CHECK(is_scalar(x.TransformsAs()));
   DEBUG_CHECK_EQUAL(x.Basis1(), x.Basis2());
   //DEBUG_CHECK(is_regular_basis(x.Basis1()))(x.Basis1());

   // handle the case of a non-regular basis
   if (!is_regular_basis(x.Basis1()))
   {
      Regularizer R(x.Basis1());
      x = RegularizeBasis12(R, x, R);
      auto Result = DiagonalizeHermitian(x);
      std::get<1>(Result) = UnregularizeBasis12(R, std::get<1>(Result), R);
      std::get<0>(Result) = UnregularizeBasis12(R, std::get<0>(Result), R);
      return Result;
   }

   IrredTensor<LinearAlgebra::DiagonalMatrix<double>, VectorBasis, VectorBasis, Tensor::DiagonalStructure> EValues(x.Basis1(), x.Basis2());

   for (auto I = iterate(x); I; ++I)
   {
      for (auto J = iterate(I); J; ++J)
      {
         DEBUG_CHECK_EQUAL(J.index1(), J.index2());
         LinearAlgebra::Vector<double> EVal = LinearAlgebra::DiagonalizeHermitian(*J);
         EValues(J.index1(), J.index2()).diagonal() = EVal;
      }
   }

   return std::make_tuple(EValues, x);
}

LinearAlgebra::Vector<double>
EigenvaluesHermitian(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                     VectorBasis, VectorBasis> const& x)
{
   LinearAlgebra::Vector<double> Result;
   for (std::size_t i = 0; i < x.Basis1().size(); ++i)
   {
      if (iterate_at(x.data(), i, i))
      {
         Result = direct_sum(Result, LinearAlgebra::EigenvaluesHermitian(x(i,i)));
      }
   }
   return Result;
}

void
InvertHPD(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>& x)
{
   DEBUG_CHECK(is_scalar(x.TransformsAs()));
   DEBUG_CHECK_EQUAL(x.Basis1(), x.Basis2());
   DEBUG_CHECK(is_regular_basis(x.Basis1()))(x.Basis1());

   // we deliberately index the matrix here: if any components are missing
   // then we have a precondition failure (the matrix would be singular).
   for (std::size_t i = 0; i < x.Basis1().size(); ++i)
   {
      InvertHPD(x(i,i));
   }
}

void
InvertGeneral(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>& x)
{
   DEBUG_CHECK(is_scalar(x.TransformsAs()));
   DEBUG_CHECK_EQUAL(x.Basis1(), x.Basis2());
   DEBUG_CHECK(is_regular_basis(x.Basis1()))(x.Basis1());

   // we deliberately index the matrix here: if any components are missing
   // then we have a precondition failure (the matrix would be singular).
   for (std::size_t i = 0; i < x.Basis1().size(); ++i)
   {
      InvertGeneral(x(i,i));
   }
}

void
InvertIrregularHPD(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>& x)
{
   DEBUG_CHECK(is_scalar(x.TransformsAs()));
   DEBUG_CHECK_EQUAL(x.Basis1(), x.Basis2());
   Regularizer R(x.Basis1());
   x = RegularizeBasis12(R, x, R);
   InvertHPD(x);
   x = UnregularizeBasis12(R, x, R);
}

void
SingularValueDecompositionRegular(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                  VectorBasis, VectorBasis> const& m,
                                  IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                  VectorBasis, VectorBasis>& U,
                                  IrredTensor<LinearAlgebra::DiagonalMatrix<double>,
                                  VectorBasis, VectorBasis, DiagonalStructure>& D,
                                     IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                                 VectorBasis, VectorBasis>& Vh)
{
   typedef IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
      VectorBasis, VectorBasis> IrredT;
   typedef IrredTensor<LinearAlgebra::DiagonalMatrix<double>,
      VectorBasis, VectorBasis, DiagonalStructure> DiagonalT;
   DEBUG_CHECK(is_scalar(m.TransformsAs()));
   DEBUG_CHECK(is_regular_basis(m.Basis1()))(m.Basis1());
   DEBUG_CHECK(is_regular_basis(m.Basis2()))(m.Basis2());
   // make the basis for the D matrix
   std::vector<int> BasisMap(m.Basis1().size(), 0);
   VectorBasis DBasis(m.GetSymmetryList());
   for (unsigned i = 0; i < m.Basis1().size(); ++i)
   {
      for (unsigned j = 0; j < m.Basis2().size(); ++j)
      {
         if (m.Basis1()[i] != m.Basis2()[j])
            continue;

         BasisMap[i] = DBasis.size();
         DBasis.push_back(m.Basis1()[i], std::min(m.Basis1().dim(i), m.Basis2().dim(j)));
      }
   }

   U = IrredT(m.Basis1(), DBasis);
   D = DiagonalT(DBasis, DBasis);
   Vh = IrredT(DBasis, m.Basis2());

   for (unsigned i = 0; i < m.Basis1().size(); ++i)
   {
      for (unsigned j = 0; j < m.Basis2().size(); ++j)
      {
         if (m.Basis1()[i] != m.Basis2()[j])
            continue;

         if (!iterate_at(m.data(), i, j))
            continue;
         set_element(U.data(), i,BasisMap[i], LinearAlgebra::Matrix<std::complex<double> >());
         set_element(Vh.data(), BasisMap[i], j, LinearAlgebra::Matrix<std::complex<double> >());
         LinearAlgebra::Vector<double> Dvec;
         LinearAlgebra::SingularValueDecomposition(m(i,j),
                                                   U(i,BasisMap[i]),
                                                   Dvec,
                                                   Vh(BasisMap[i], j));
         set_element(D, BasisMap[i], BasisMap[i], LinearAlgebra::DiagonalMatrix<double>(Dvec));
      }
   }
}

void
SingularValueDecompositionRegular(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                  VectorBasis, VectorBasis> const& m,
                                  IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                  VectorBasis, VectorBasis>& U,
                                  IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                              VectorBasis, VectorBasis>& D,
                                     IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                                 VectorBasis, VectorBasis>& Vh)
{
   IrredTensor<LinearAlgebra::DiagonalMatrix<double>, VectorBasis, VectorBasis, DiagonalStructure> Temp;
   SingularValueDecompositionRegular(m, U, Temp, Vh);
   D = Temp;
}

void
SingularValueDecomposition(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                       VectorBasis, VectorBasis> const& m,
                           IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                       VectorBasis, VectorBasis>& U,
                           IrredTensor<LinearAlgebra::DiagonalMatrix<double>,
                           VectorBasis, VectorBasis, DiagonalStructure>& D,
                           IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                           VectorBasis, VectorBasis>& Vh)
{
   typedef IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                       VectorBasis, VectorBasis> MatrixOperator;

   // If the basis is already regular, don't bother with the extra work
   if (is_regular_basis(m.Basis1()) && is_regular_basis(m.Basis2()))
   {
      SingularValueDecompositionRegular(m, U, D, Vh);
      return;
   }
   // else

   Regularizer R1(m.Basis1());
   Regularizer R2(m.Basis2());

   SingularValueDecompositionRegular(RegularizeBasis12(R1, m, R2), U, D, Vh);
   U = UnregularizeBasis1(R1, U);
   Vh = UnregularizeBasis2(Vh, R2);
}

void
SingularValueDecomposition(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                       VectorBasis, VectorBasis> const& m,
                           IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                       VectorBasis, VectorBasis>& U,
                           IrredTensor<LinearAlgebra::DiagonalMatrix<double>,
                                       VectorBasis, VectorBasis>& D,
                           IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                           VectorBasis, VectorBasis>& Vh)
{
   IrredTensor<LinearAlgebra::DiagonalMatrix<double>, VectorBasis, VectorBasis, DiagonalStructure> Temp;
   SingularValueDecomposition(m, U, Temp, Vh);
   D = Temp;
}

void
SingularValueDecomposition(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                       VectorBasis, VectorBasis> const& m,
                           IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                       VectorBasis, VectorBasis>& U,
                           IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                       VectorBasis, VectorBasis>& D,
                           IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                           VectorBasis, VectorBasis>& Vh)
{
   IrredTensor<LinearAlgebra::DiagonalMatrix<double>, VectorBasis, VectorBasis, DiagonalStructure> Temp;
   SingularValueDecomposition(m, U, Temp, Vh);
   D = Temp;
}

LinearAlgebra::Vector<double>
SingularValuesRegular(IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis> const& m, SingularValueNormalization n)
{
   // First pass - calculate the total number of singular vectors
   std::vector<int> Offset(m.Basis1().size()+1, 0); // offset into the singular value array, as a function of 'i'
   int Dim = 0; // number of singular values
   for (unsigned i = 0; i < m.Basis1().size(); ++i)
   {
      Offset[i] = Dim;
      for (unsigned j = 0; j < m.Basis2().size(); ++j)
      {
         if (m.Basis1()[i] != m.Basis2()[j])
            continue;

         // if the matrix is zero then so are the singular values, so skip it.
         if (!iterate_at(m.data(), i, j))
            continue;

         Dim += std::min(m.Basis1().dim(i), m.Basis2().dim(j));
      }
   }
   Offset.back() = Dim;

   LinearAlgebra::Vector<double> Result(Dim);

   // Second pass, calculate the singular values
   for (unsigned i = 0; i < m.Basis1().size(); ++i)
   {
      //
      for (unsigned j = 0; j < m.Basis2().size(); ++j)
      {
         if (m.Basis1()[i] != m.Basis2()[j])
            continue;

         // if the matrix is zero then so are the singular values, so skip it.
         if (!iterate_at(m.data(), i, j))
            continue;

         auto Proxy = Result[LinearAlgebra::range(Offset[i], Offset[i+1])];
         LinearAlgebra::SingularValues(m(i,j), Proxy);
         if (n == SingularValueNormalization::QDim)
         {
            Result[LinearAlgebra::range(Offset[i], Offset[i+1])] *= degree(m.Basis1()[i]);
         }
         else if (n == SingularValueNormalization::SqrtQDim)
         {
            Result[LinearAlgebra::range(Offset[i], Offset[i+1])] *= std::sqrt(degree(m.Basis1()[i]));
         }
      }
   }
   std::sort(data(Result), data(Result)+Dim, std::greater<>());
   return Result;
}

LinearAlgebra::Vector<double>
SingularValues(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
               VectorBasis, VectorBasis> const& m, SingularValueNormalization n)
{
   // If the basis is already regular, don't bother with the extra work
   if (is_regular_basis(m.Basis1()) && is_regular_basis(m.Basis2()))
   {
      return SingularValuesRegular(m, n);
   }
   // else
   return SingularValuesRegular(RegularizeBasis12(Regularizer(m.Basis1()), m, Regularizer(m.Basis2())), n);
}

void
SingularValueDecomposition(IrredTensor<std::complex<double>, BasisList, BasisList> const& m,
                           IrredTensor<std::complex<double>, BasisList, BasisList>& U,
                           IrredTensor<std::complex<double>, BasisList, BasisList>& D,
                           IrredTensor<std::complex<double>, BasisList, BasisList>& Vh)
{
   typedef IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                       VectorBasis, VectorBasis> MatrixOperator;

   typedef IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, BasisList> MixedOperator;

   MixedOperator U1 = Regularize(m.Basis1());
   MixedOperator U2 = Regularize(m.Basis2());

   MatrixOperator UMatrix, DMatrix, VhMatrix;
   SingularValueDecompositionRegular(triple_prod(U1, m, herm(U2)), UMatrix, DMatrix, VhMatrix);

   MixedOperator Splitter = SplitBasis(DMatrix.Basis1());

   U = map_1x1_operator(triple_prod(herm(U1), UMatrix, Splitter));
   D = map_1x1_operator(triple_prod(herm(Splitter), DMatrix, Splitter));
   Vh = map_1x1_operator(triple_prod(herm(Splitter), VhMatrix, U2));
   return;
}

void
SingularValueDecompositionFullRegular(IrredTensor<LinearAlgebra::Matrix<std::complex<double>>,
                                      VectorBasis, VectorBasis> const& m,
                                      IrredTensor<LinearAlgebra::Matrix<std::complex<double>>,
                                      VectorBasis, VectorBasis>& U,
                                      IrredTensor<LinearAlgebra::DiagonalMatrix<double>,
                                      VectorBasis, VectorBasis, DiagonalStructure>& D,
   IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
               VectorBasis, VectorBasis>& Vh)
{
   typedef IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
      VectorBasis, VectorBasis> IrredT;
   typedef IrredTensor<LinearAlgebra::DiagonalMatrix<double>,
      VectorBasis, VectorBasis, DiagonalStructure> DiagonalT;
   DEBUG_CHECK(is_scalar(m.TransformsAs()));
   DEBUG_CHECK(is_regular_basis(m.Basis1()))(m.Basis1());
   DEBUG_CHECK(is_regular_basis(m.Basis2()))(m.Basis2());
   // make the basis for the D matrix
   // For states that exist in both bases, BasisMap tells us where the state from Basis1
   // occurs in the DBasis
   std::vector<int> BasisMap(m.Basis1().size(), 0);
   VectorBasis DBasis(m.GetSymmetryList());

   // Find the superset of the basis 1 and basis 2.  We need to handle separately the
   // case where a quantum number exists in only one basis
   for (unsigned i = 0; i < m.Basis1().size(); ++i)
   {
      for (unsigned j = 0; j < m.Basis2().size(); ++j)
      {
         if (m.Basis1()[i] != m.Basis2()[j])
            continue;

         BasisMap[i] = DBasis.size();
         DBasis.push_back(m.Basis1()[i], std::max(m.Basis1().dim(i), m.Basis2().dim(j)));
      }
   }

   int const Mark = DBasis.size(); // all states in DBasis from here only exist in one basis
   // Find states that exist only in Basis1
   for (unsigned i = 0; i < m.Basis1().size(); ++i)
   {
      unsigned j = 0;
      while (j < m.Basis2().size() && m.Basis1()[i] != m.Basis2()[j])
         ++j;
      if (j == m.Basis2().size())
      {
         // state only exists in basis 1
         DBasis.push_back(m.Basis1()[i], m.Basis1().dim(i));
      }
   }
   // Find states that only exist in Basis2
   for (unsigned j = 0; j < m.Basis2().size(); ++j)
   {
      unsigned i = 0;
      while (i < m.Basis1().size() && m.Basis1()[i] != m.Basis2()[j])
         ++i;
      if (i == m.Basis1().size())
      {
         // state only exists in basis 2
         DBasis.push_back(m.Basis2()[j], m.Basis2().dim(j));
      }
   }

   U = IrredT(m.Basis1(), DBasis);
   D = DiagonalT(DBasis, DBasis);
   Vh = IrredT(DBasis, m.Basis2());

   for (unsigned i = 0; i < m.Basis1().size(); ++i)
   {
      for (unsigned j = 0; j < m.Basis2().size(); ++j)
      {
         if (m.Basis1()[i] != m.Basis2()[j])
            continue;

         if (!iterate_at(m.data(), i, j))
         {
            PANIC("SVD - not implemented SVD of a component that doesn't exist - please fix!");
         }
         set_element(U.data(), i,BasisMap[i], LinearAlgebra::Matrix<std::complex<double> >());
         set_element(Vh.data(), BasisMap[i], j, LinearAlgebra::Matrix<std::complex<double> >());
         LinearAlgebra::DiagonalMatrix<double> Dvec;
         LinearAlgebra::SingularValueDecompositionFull(m(i,j),
                                                       U(i,BasisMap[i]),
                                                       Dvec,
                                                       Vh(BasisMap[i], j));
         set_element(D, BasisMap[i], BasisMap[i], Dvec);
      }
   }

   // Now for the states that only exist in Basis1, set the matrix elements
   int n = Mark;
   for (unsigned i = 0; i < m.Basis1().size(); ++i)
   {
      unsigned j = 0;
      while (j < m.Basis2().size() && m.Basis1()[i] != m.Basis2()[j])
         ++j;
      if (j == m.Basis2().size())
      {
         // state only exists in basis 1
         U(i,n) = LinearAlgebra::DiagonalMatrix<double>(m.Basis1().dim(i), m.Basis1().dim(i), 1.0);
         D(n,n) = LinearAlgebra::DiagonalMatrix<double>(m.Basis1().dim(i), m.Basis1().dim(i), 0.0);
         ++n;
      }
   }
   // Now for states that only exist in Basis2
   for (unsigned j = 0; j < m.Basis2().size(); ++j)
   {
      unsigned i = 0;
      while (i < m.Basis1().size() && m.Basis1()[i] != m.Basis2()[j])
         ++i;
      if (i == m.Basis1().size())
      {
         // state only exists in basis 2
         Vh(n,j) = LinearAlgebra::DiagonalMatrix<double>(m.Basis2().dim(j), m.Basis2().dim(j), 1.0);
         D(n,n) = LinearAlgebra::DiagonalMatrix<double>(m.Basis2().dim(j), m.Basis2().dim(j), 0.0);
         ++n;
      }
   }
   U.debug_check_structure();
   D.debug_check_structure();
   Vh.debug_check_structure();
}

void
SingularValueDecompositionKeepBasis1Regular(IrredTensor<LinearAlgebra::Matrix<std::complex<double>>,
                                      VectorBasis, VectorBasis> const& m,
                                      IrredTensor<LinearAlgebra::Matrix<std::complex<double>>,
                                      VectorBasis, VectorBasis>& U,
                                      IrredTensor<LinearAlgebra::DiagonalMatrix<double>,
                                      VectorBasis, VectorBasis, DiagonalStructure>& D,
   IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
               VectorBasis, VectorBasis>& Vh)
{
   typedef IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
      VectorBasis, VectorBasis> IrredT;
   typedef IrredTensor<LinearAlgebra::DiagonalMatrix<double>,
      VectorBasis, VectorBasis, DiagonalStructure> DiagonalT;
   DEBUG_CHECK(is_scalar(m.TransformsAs()));
   DEBUG_CHECK(is_regular_basis(m.Basis1()))(m.Basis1());
   DEBUG_CHECK(is_regular_basis(m.Basis2()))(m.Basis2());
   // make the basis for the D matrix
   // For states that exist in both bases, BasisMap tells us where the state from Basis1
   // occurs in the DBasis
   std::vector<int> BasisMap(m.Basis1().size(), 0);
   VectorBasis DBasis(m.GetSymmetryList());

   // Find the superset of the basis 1 and basis 2.  We need to handle separately the
   // case where a quantum number exists in only one basis
   for (unsigned i = 0; i < m.Basis1().size(); ++i)
   {
      for (unsigned j = 0; j < m.Basis2().size(); ++j)
      {
         if (m.Basis1()[i] != m.Basis2()[j])
            continue;

         BasisMap[i] = DBasis.size();
         DBasis.push_back(m.Basis1()[i], m.Basis1().dim(i));
      }
   }

   int const Mark = DBasis.size(); // all states in DBasis from here only exist in one basis
   // Find states that exist only in Basis1
   for (unsigned i = 0; i < m.Basis1().size(); ++i)
   {
      unsigned j = 0;
      while (j < m.Basis2().size() && m.Basis1()[i] != m.Basis2()[j])
         ++j;
      if (j == m.Basis2().size())
      {
         // state only exists in basis 1
         DBasis.push_back(m.Basis1()[i], m.Basis1().dim(i));
      }
   }

   U = IrredT(m.Basis1(), DBasis);
   D = DiagonalT(DBasis, DBasis);
   Vh = IrredT(DBasis, m.Basis2());

   for (unsigned i = 0; i < m.Basis1().size(); ++i)
   {
      for (unsigned j = 0; j < m.Basis2().size(); ++j)
      {
         if (m.Basis1()[i] != m.Basis2()[j])
            continue;

         if (!iterate_at(m.data(), i, j))
            continue;
         set_element(U.data(), i,BasisMap[i], LinearAlgebra::Matrix<std::complex<double> >());
         set_element(Vh.data(), BasisMap[i], j, LinearAlgebra::Matrix<std::complex<double> >());
         LinearAlgebra::DiagonalMatrix<double> Dvec;
         LinearAlgebra::SingularValueDecompositionFullLeft(m(i,j),
                                                           U(i,BasisMap[i]),
                                                           Dvec,
                                                           Vh(BasisMap[i], j));
         set_element(D, BasisMap[i], BasisMap[i], Dvec);
      }
   }

   // Now for the quantum numbers that only exist in Basis1, set the matrix elements
   int n = Mark;
   for (unsigned i = 0; i < m.Basis1().size(); ++i)
   {
      unsigned j = 0;
      while (j < m.Basis2().size() && m.Basis1()[i] != m.Basis2()[j])
         ++j;
      if (j == m.Basis2().size())
      {
         // state only exists in basis 1
         U(i,n) = LinearAlgebra::DiagonalMatrix<double>(m.Basis1().dim(i), m.Basis1().dim(i), 1.0);
         D(n,n) = LinearAlgebra::DiagonalMatrix<double>(m.Basis1().dim(i), m.Basis1().dim(i), 0.0);
         ++n;
      }
   }

   U.debug_check_structure();
   D.debug_check_structure();
   Vh.debug_check_structure();
}

void
SingularValueDecompositionKeepBasis2Regular(IrredTensor<LinearAlgebra::Matrix<std::complex<double>>,
                                      VectorBasis, VectorBasis> const& m,
                                      IrredTensor<LinearAlgebra::Matrix<std::complex<double>>,
                                      VectorBasis, VectorBasis>& U,
                                      IrredTensor<LinearAlgebra::DiagonalMatrix<double>,
                                      VectorBasis, VectorBasis, DiagonalStructure>& D,
   IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
               VectorBasis, VectorBasis>& Vh)
{
   typedef IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
      VectorBasis, VectorBasis> IrredT;
   typedef IrredTensor<LinearAlgebra::DiagonalMatrix<double>,
      VectorBasis, VectorBasis, DiagonalStructure> DiagonalT;
   DEBUG_CHECK(is_scalar(m.TransformsAs()));
   DEBUG_CHECK(is_regular_basis(m.Basis1()))(m.Basis1());
   DEBUG_CHECK(is_regular_basis(m.Basis2()))(m.Basis2());
   // make the basis for the D matrix
   // For states that exist in both bases, BasisMap tells us where the state from Basis1
   // occurs in the DBasis
   std::vector<int> BasisMap(m.Basis1().size(), 0);
   VectorBasis DBasis(m.GetSymmetryList());

   // Find the superset of the basis 1 and basis 2.  We need to handle separately the
   // case where a quantum number exists in only one basis
   for (unsigned i = 0; i < m.Basis1().size(); ++i)
   {
      for (unsigned j = 0; j < m.Basis2().size(); ++j)
      {
         if (m.Basis1()[i] != m.Basis2()[j])
            continue;

         BasisMap[i] = DBasis.size();
         DBasis.push_back(m.Basis1()[i], m.Basis2().dim(j));
      }
   }

   int const Mark = DBasis.size(); // all states in DBasis from here only exist in one basis
   // Find states that only exist in Basis2
   for (unsigned j = 0; j < m.Basis2().size(); ++j)
   {
      unsigned i = 0;
      while (i < m.Basis1().size() && m.Basis1()[i] != m.Basis2()[j])
         ++i;
      if (i == m.Basis1().size())
      {
         // state only exists in basis 2
         DBasis.push_back(m.Basis2()[j], m.Basis2().dim(j));
      }
   }

   U = IrredT(m.Basis1(), DBasis);
   D = DiagonalT(DBasis, DBasis);
   Vh = IrredT(DBasis, m.Basis2());

   for (unsigned i = 0; i < m.Basis1().size(); ++i)
   {
      for (unsigned j = 0; j < m.Basis2().size(); ++j)
      {
         if (m.Basis1()[i] != m.Basis2()[j])
            continue;

         if (!iterate_at(m.data(), i, j))
            continue;
         set_element(U.data(), i,BasisMap[i], LinearAlgebra::Matrix<std::complex<double> >());
         set_element(Vh.data(), BasisMap[i], j, LinearAlgebra::Matrix<std::complex<double> >());
         LinearAlgebra::DiagonalMatrix<double> Dvec;
         LinearAlgebra::SingularValueDecompositionFullRight(m(i,j),
                                                            U(i,BasisMap[i]),
                                                            Dvec,
                                                            Vh(BasisMap[i], j));
         set_element(D, BasisMap[i], BasisMap[i], Dvec);
      }
   }

   // Now for states that only exist in Basis2
   int n = Mark;
   for (unsigned j = 0; j < m.Basis2().size(); ++j)
   {
      unsigned i = 0;
      while (i < m.Basis1().size() && m.Basis1()[i] != m.Basis2()[j])
         ++i;
      if (i == m.Basis1().size())
      {
         // state only exists in basis 2
         Vh(n,j) = LinearAlgebra::DiagonalMatrix<double>(m.Basis2().dim(j), m.Basis2().dim(j), 1.0);
         D(n,n) = LinearAlgebra::DiagonalMatrix<double>(m.Basis2().dim(j), m.Basis2().dim(j), 0.0);
         ++n;
      }
   }
   U.debug_check_structure();
   D.debug_check_structure();
   Vh.debug_check_structure();
}

void
SingularValueDecompositionFull(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                       VectorBasis, VectorBasis> const& m,
                           IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                       VectorBasis, VectorBasis>& U,
                           IrredTensor<LinearAlgebra::DiagonalMatrix<double>,
                           VectorBasis, VectorBasis, DiagonalStructure>& D,
                           IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                           VectorBasis, VectorBasis>& Vh)
{
   typedef IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                       VectorBasis, VectorBasis> MatrixOperator;

   // If the basis is already regular, don't bother with the extra work
   if (is_regular_basis(m.Basis1()) && is_regular_basis(m.Basis2()))
   {
      SingularValueDecompositionFullRegular(m, U, D, Vh);
      return;
   }
   // else

   Regularizer R1(m.Basis1());
   Regularizer R2(m.Basis2());

   SingularValueDecompositionFullRegular(RegularizeBasis12(R1, m, R2), U, D, Vh);
   U = UnregularizeBasis1(R1, U);
   Vh = UnregularizeBasis2(Vh, R2);
}

void
SingularValueDecompositionKeepBasis1(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                       VectorBasis, VectorBasis> const& m,
                           IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                       VectorBasis, VectorBasis>& U,
                           IrredTensor<LinearAlgebra::DiagonalMatrix<double>,
                           VectorBasis, VectorBasis, DiagonalStructure>& D,
                           IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                           VectorBasis, VectorBasis>& Vh)
{
   typedef IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                       VectorBasis, VectorBasis> MatrixOperator;

   // If the basis is already regular, don't bother with the extra work
   if (is_regular_basis(m.Basis1()) && is_regular_basis(m.Basis2()))
   {
      SingularValueDecompositionKeepBasis1Regular(m, U, D, Vh);
      return;
   }
   // else

   Regularizer R1(m.Basis1());
   Regularizer R2(m.Basis2());

   SingularValueDecompositionKeepBasis1Regular(RegularizeBasis12(R1, m, R2), U, D, Vh);
   U = UnregularizeBasis1(R1, U);
   Vh = UnregularizeBasis2(Vh, R2);
}

void
SingularValueDecompositionKeepBasis2(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                       VectorBasis, VectorBasis> const& m,
                           IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                       VectorBasis, VectorBasis>& U,
                           IrredTensor<LinearAlgebra::DiagonalMatrix<double>,
                           VectorBasis, VectorBasis, DiagonalStructure>& D,
                           IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                           VectorBasis, VectorBasis>& Vh)
{
   typedef IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                       VectorBasis, VectorBasis> MatrixOperator;

   // If the basis is already regular, don't bother with the extra work
   if (is_regular_basis(m.Basis1()) && is_regular_basis(m.Basis2()))
   {
      SingularValueDecompositionKeepBasis2Regular(m, U, D, Vh);
      return;
   }
   // else

   Regularizer R1(m.Basis1());
   Regularizer R2(m.Basis2());

   SingularValueDecompositionKeepBasis1Regular(RegularizeBasis12(R1, m, R2), U, D, Vh);
   U = UnregularizeBasis1(R1, U);
   Vh = UnregularizeBasis2(Vh, R2);
}

IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
InvertDiagonal(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
               VectorBasis, VectorBasis> const& Op, double Cutoff)
{
   IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>
      Result(Op.Basis2(), Op.Basis1());

   for (unsigned i = 0; i < Op.Basis1().size(); ++i)
   {
      //TRACE(i);
      if (iterate_at(Op.data(), i, i))
      {
         Result(i,i) = LinearAlgebra::Matrix<double>(Op.Basis2().dim(i), Op.Basis1().dim(i), 0.0);
         for (int j = 0; j < std::min(Op.Basis1().dim(i), Op.Basis2().dim(i)); ++j)
         {
            std::complex<double> x = Op(i,i)(j,j);
            if (LinearAlgebra::norm_frob(x) < Cutoff)
               x = 0;
            else
               x = 1.0 / x;
            Result(i,i)(j,j) = x;
         }
      }
      //TRACE(Result(i,i));
   }

   return Result;
}

//
// QR Factorization
//

std::pair<IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>,
         IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>>
QR_FactorizeFullRegular(IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis> m)
{
   CHECK(is_scalar(m.TransformsAs()));
   DEBUG_CHECK(is_regular_basis(m.Basis1()));
   DEBUG_CHECK(is_regular_basis(m.Basis2()));

   IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis> Q(m.Basis1(), m.Basis1());
   IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis> R(m.Basis1(), m.Basis2());

   for (int i = 0; i < m.Basis1().size(); ++i)
   {
      int j = m.Basis2().find_first(m.Basis1()[i]);
      // If the quantum number doesn't exist in Basis2(), then the corresponding matrix is n by zero
      if (j < 0)
         continue;

      //CHECK(m.Basis1().dim(i) >= m.Basis2().dim(j))("QR decomposition is not defined for m*n matrix with n > m!");
      auto I = iterate_at(m.data(), i,j);
      if (!I)
      {
         // the component is zero.  We could set the Q matrix to the identity, but probably better to use a random unitary
         Q(i,i) = LinearAlgebra::random_unitary<std::complex<double>>(m.Basis1().dim(i), m.Basis1().dim(i));
         R(i,j) = LinearAlgebra::Matrix<std::complex<double>>(m.Basis1().dim(i), m.Basis2().dim(j), 0.0);
      }
      else
      {
         R(i,j) = std::move(*I);
         Q(i,i) = QR_FactorizeFull(R(i,j));
      }
   }
   return std::make_pair(std::move(Q), std::move(R));
}

std::pair<IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>,
         IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>>
QR_FactorizeFull(IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis> m)
{
   if (is_regular_basis(m.Basis1()) && is_regular_basis(m.Basis2()))
   {
      return QR_FactorizeFullRegular(std::move(m));
   }
   // else

   Regularizer R1(m.Basis1());
   Regularizer R2(m.Basis2());

   auto QR = QR_FactorizeFullRegular(RegularizeBasis12(R1, std::move(m), R2));
   QR.first = UnregularizeBasis1(R1, QR.first);
   QR.second = UnregularizeBasis2(QR.second, R2);
   return QR;
}

std::pair<IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>,
         IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>>
QR_FactorizeRegular(IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis> m)
{
   CHECK(is_scalar(m.TransformsAs()));
   DEBUG_CHECK(is_regular_basis(m.Basis1()));
   DEBUG_CHECK(is_regular_basis(m.Basis2()));

   // Make a basis for L.Basis2() and Q.Basis1().
   std::vector<int> B1Index, B2Index;
   VectorBasis B3;
   std::tie(B3, B1Index, B2Index) = MakeMinimalBasis(m.Basis1(), m.Basis2());
   IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis> Q(m.Basis1(), B3);
   IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis> R(B3, m.Basis2());

   for (int i = 0; i < m.Basis1().size(); ++i)
   {
      int j = m.Basis2().find_first(m.Basis1()[i]);
      // If the quantum number doesn't exist in Basis2(), then the corresponding matrix is n by zero
      if (j < 0)
         continue;

      int k = B1Index[i];
      DEBUG_CHECK_EQUAL(k, B2Index[j]);

      auto I = iterate_at(m.data(), i,j);
      if (!I)
      {
         // the component is zero.  We set the corresponding Q matrix to a random unitary
         Q(i,k) = LinearAlgebra::random_unitary<std::complex<double>>(m.Basis1().dim(i), B3.dim(k));
         R(k,j) = LinearAlgebra::Matrix<std::complex<double>>(B3.dim(k), m.Basis2().dim(j), 0.0);
      }
      else
      {
         std::tie(Q(i,k), R(k,j)) = QR_FactorizeThin(std::move(*I));
      }
   }
   return std::make_pair(std::move(Q), std::move(R));
}

std::pair<IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>,
         IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>>
QR_Factorize(IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis> m)
{
   if (is_regular_basis(m.Basis1()) && is_regular_basis(m.Basis2()))
   {
      return QR_FactorizeRegular(std::move(m));
   }
   // else

   Regularizer R1(m.Basis1());
   Regularizer R2(m.Basis2());

   auto QR = QR_FactorizeRegular(RegularizeBasis12(R1, std::move(m), R2));
   QR.first = UnregularizeBasis1(R1, QR.first);
   QR.second = UnregularizeBasis2(QR.second, R2);
   return QR;
}

//
// LQ Factorization
//

std::pair<IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>,
         IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>>
LQ_FactorizeFullRegular(IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis> m)
{
   CHECK(is_scalar(m.TransformsAs()));
   DEBUG_CHECK(is_regular_basis(m.Basis1()));
   DEBUG_CHECK(is_regular_basis(m.Basis2()));

   IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis> L(m.Basis1(), m.Basis2());
   IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis> Q(m.Basis2(), m.Basis2());

   for (int i = 0; i < m.Basis1().size(); ++i)
   {
      int j = m.Basis2().find_first(m.Basis1()[i]);
      // If the quantum number doesn't exist in Basis2(), then the corresponding matrix is n by zero
      if (j < 0)
         continue;

      //CHECK(m.Basis1().dim(i) >= m.Basis2().dim(j))("QR decomposition is not defined for m*n matrix with n > m!");
      auto I = iterate_at(m.data(), i,j);
      if (!I)
      {
         // the component is zero.  We could set the Q matrix to the identity, but probably better to use a random unitary
         L(i,j) = LinearAlgebra::Matrix<std::complex<double>>(m.Basis1().dim(i), m.Basis2().dim(j), 0.0);
         Q(j,j) = LinearAlgebra::random_unitary<std::complex<double>>(m.Basis2().dim(i), m.Basis1().dim(j));
      }
      else
      {
         L(i,j) = std::move(*I);
         Q(j,j) = LQ_FactorizeFull(L(i,j));
      }
   }
   return std::make_pair(std::move(L), std::move(Q));
}

std::pair<IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>,
         IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>>
LQ_FactorizeFull(IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis> m)
{
   if (is_regular_basis(m.Basis1()) && is_regular_basis(m.Basis2()))
   {
      return LQ_FactorizeFullRegular(std::move(m));
   }
   // else

   Regularizer R1(m.Basis1());
   Regularizer R2(m.Basis2());

   auto LQ = LQ_FactorizeFullRegular(RegularizeBasis12(R1, std::move(m), R2));
   LQ.first = UnregularizeBasis1(R1, LQ.first);
   LQ.second = UnregularizeBasis2(LQ.second, R2);
   return LQ;
}

std::pair<IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>,
         IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>>
LQ_FactorizeRegular(IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis> m)
{
   CHECK(is_scalar(m.TransformsAs()));
   DEBUG_CHECK(is_regular_basis(m.Basis1()));
   DEBUG_CHECK(is_regular_basis(m.Basis2()));

   // Make a basis for L.Basis2() and Q.Basis1().
   std::vector<int> B1Index, B2Index;
   VectorBasis B3;
   std::tie(B3, B1Index, B2Index) = MakeMinimalBasis(m.Basis1(), m.Basis2());
   IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis> L(m.Basis1(), B3);
   IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis> Q(B3, m.Basis2());

   for (int i = 0; i < m.Basis1().size(); ++i)
   {
      int j = m.Basis2().find_first(m.Basis1()[i]);
      // If the quantum number doesn't exist in Basis2(), then the corresponding matrix is n by zero
      if (j < 0)
         continue;

      int k = B1Index[i];
      DEBUG_CHECK_EQUAL(k, B2Index[j]);

      auto I = iterate_at(m.data(), i,j);
      if (!I)
      {
         // the component is zero.  We set the corresponding Q matrix to a random unitary
         L(i,k) = LinearAlgebra::Matrix<std::complex<double>>(m.Basis1().dim(i), B3.dim(k), 0.0);
         Q(k,j) = LinearAlgebra::random_unitary<std::complex<double>>(B3.dim(k), m.Basis2().dim(j));
      }
      else
      {
         std::tie(L(i,k), Q(k,j)) = LQ_FactorizeThin(std::move(*I));
      }
   }
   return std::make_pair(std::move(L), std::move(Q));
}

std::pair<IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>,
         IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>>
LQ_Factorize(IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis> m)
{
   if (is_regular_basis(m.Basis1()) && is_regular_basis(m.Basis2()))
   {
      return LQ_FactorizeRegular(std::move(m));
   }
   // else

   Regularizer R1(m.Basis1());
   Regularizer R2(m.Basis2());

   auto LQ = LQ_FactorizeRegular(RegularizeBasis12(R1, std::move(m), R2));
   LQ.first = UnregularizeBasis1(R1, LQ.first);
   LQ.second = UnregularizeBasis2(LQ.second, R2);
   return LQ;
}

IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
InvertDiagonal(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
               VectorBasis, VectorBasis> const& Op)
{
   return InvertDiagonal(Op, std::sqrt(std::numeric_limits<double>::epsilon()));
}

// IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
// SqrtDiagonal(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
//                VectorBasis, VectorBasis> const& Op, double Tol)
// {
//    IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
//       Result(Op.Basis2(), Op.Basis1());
//
//    for (unsigned i = 0; i < Op.Basis1().size(); ++i)
//    {
//       if (iterate_at(Op.data(), i, i))
//       {
//          Result(i,i) = LinearAlgebra::Matrix<double>(Op.Basis2().dim(i), Op.Basis1().dim(i), 0.0);
//          for (int j = 0; j < std::min(Op.Basis1().dim(i), Op.Basis2().dim(i)); ++j)
//          {
//             std::complex<double> x = Op(i,i)(j,j);
//             CHECK(LinearAlgebra::norm_frob(x.imag()) <= Tol);
//             if (x.real() <= 0)
//             {
//                CHECK(norm_frob(x.real()) <= Tol)(x.real())(Tol)(Op);
//                Result(i,i)(j,j) = 0.0;
//             }
//             else
//                Result(i,i)(j,j) = std::sqrt(x.real());
//          }
//       }
//       //TRACE(Result(i,i));
//    }
//
//    return Result;
// }

IrredTensor<LinearAlgebra::DiagonalMatrix<std::complex<double> >, VectorBasis, VectorBasis, Tensor::DiagonalStructure>
SqrtDiagonal(IrredTensor<LinearAlgebra::DiagonalMatrix<std::complex<double> >,
               VectorBasis, VectorBasis, Tensor::DiagonalStructure> const& Op, double Tol)
{
   IrredTensor<LinearAlgebra::DiagonalMatrix<std::complex<double>>, VectorBasis, VectorBasis, Tensor::DiagonalStructure>
      Result(Op.Basis1(), Op.Basis2());

   for (unsigned i = 0; i < Op.Basis1().size(); ++i)
   {
      Result(i,i) = LinearAlgebra::DiagonalMatrix<double>(Op.Basis1().dim(i), Op.Basis2().dim(i), 0.0);
      for (int j = 0; j < std::min(Op.Basis1().dim(i), Op.Basis2().dim(i)); ++j)
      {
         std::complex<double> x = Op(i,i)(j,j);
         CHECK(LinearAlgebra::norm_frob(x.imag()) <= Tol);
         if (x.real() <= 0)
         {
            CHECK(norm_frob(x.real()) <= Tol)(x.real())(Tol)(Op);
            Result(i,i)(j,j) = 0.0;
         }
         else
            Result(i,i)(j,j) = std::sqrt(x.real());
      }
      //TRACE(Result(i,i));
   }

   return Result;
}

IrredTensor<LinearAlgebra::DiagonalMatrix<double>, VectorBasis, VectorBasis, Tensor::DiagonalStructure>
SqrtDiagonal(IrredTensor<LinearAlgebra::DiagonalMatrix<double>,
               VectorBasis, VectorBasis, Tensor::DiagonalStructure> const& Op, double Tol)
{
   IrredTensor<LinearAlgebra::DiagonalMatrix<double>, VectorBasis, VectorBasis, Tensor::DiagonalStructure>
      Result(Op.Basis1(), Op.Basis2());

   for (unsigned i = 0; i < Op.Basis1().size(); ++i)
   {
      Result(i,i) = LinearAlgebra::DiagonalMatrix<double>(Op.Basis1().dim(i), Op.Basis2().dim(i), 0.0);
      for (int j = 0; j < std::min(Op.Basis1().dim(i), Op.Basis2().dim(i)); ++j)
      {
         double x = Op(i,i)(j,j);
         if (x <= 0)
         {
            CHECK(std::abs(x) <= Tol)(x)(Tol)(Op);
            Result(i,i)(j,j) = 0.0;
         }
         else
            Result(i,i)(j,j) = std::sqrt(x);
      }
      //TRACE(Result(i,i));
   }

   return Result;
}

// CholeskyFactorizeUpper

IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
CholeskyFactorizeUpperRegular(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                              VectorBasis, VectorBasis> const& m)
{
   DEBUG_PRECONDITION_EQUAL(m.Basis1(), m.Basis2());

   typedef IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
      VectorBasis, VectorBasis> TensorType;

   TensorType Result = m;

   for (iterator<TensorType>::type I = iterate(Result); I; ++I)
   {
      for (inner_iterator<TensorType>::type J = iterate(I); J; ++J)
      {
         LinearAlgebra::Matrix<std::complex<double> > x = *J;
         CholeskyFactorizeUpper(x);
         // zero out the lower triangular part
         for (unsigned i = 1; i < size1(x); ++i)
         {
            for (unsigned j = 0; j < i; ++j)
            {
               x(i,j) = 0.0;
            }
         }
         *J = x;
      }
   }
   return Result;
}

IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
CholeskyFactorizeUpper(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                       VectorBasis, VectorBasis> const& m)
{
   DEBUG_PRECONDITION_EQUAL(m.Basis1(), m.Basis2());

   if (is_regular_basis(m.Basis1()))
      return CholeskyFactorizeUpperRegular(m);
   // else
   Regularizer R(m.Basis1());
   return UnregularizeBasis12(R, CholeskyFactorizeUpperRegular(RegularizeBasis12(R, m, R)), R);
}

// CholeskyFactorizeLower

IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
CholeskyFactorizeLowerRegular(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                              VectorBasis, VectorBasis> const& m)
{
   DEBUG_PRECONDITION_EQUAL(m.Basis1(), m.Basis2());

   typedef IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
      VectorBasis, VectorBasis> TensorType;

   TensorType Result = m;

   for (iterator<TensorType>::type I = iterate(Result); I; ++I)
   {
      for (inner_iterator<TensorType>::type J = iterate(I); J; ++J)
      {
         LinearAlgebra::Matrix<std::complex<double> > x = *J;
         CholeskyFactorizeLower(x);
         // zero out the upper triangular part
         unsigned const s2 = size2(x);
         for (unsigned i = 0; i < size1(x); ++i)
         {
            for (unsigned j = i+1; j < s2; ++j)
            {
               x(i,j) = 0.0;
            }
         }
         *J = x;
      }
   }
   return Result;
}

IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
CholeskyFactorizeLower(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                       VectorBasis, VectorBasis> const& m)
{
   DEBUG_PRECONDITION_EQUAL(m.Basis1(), m.Basis2());

   if (is_regular_basis(m.Basis1()))
      return CholeskyFactorizeLowerRegular(m);
   // else
   Regularizer R(m.Basis1());
   return UnregularizeBasis12(R, CholeskyFactorizeLowerRegular(RegularizeBasis12(R, m, R)), R);
}

IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
SingularFactorizeRegular(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                         VectorBasis, VectorBasis> const& m)
{
   DEBUG_PRECONDITION_EQUAL(m.Basis1(), m.Basis2());

   typedef IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
      VectorBasis, VectorBasis> TensorType;

   TensorType Result = m;

   for (iterator<TensorType>::type I = iterate(Result); I; ++I)
   {
      for (inner_iterator<TensorType>::type J = iterate(I); J; ++J)
      {
         *J = SingularFactorize(*J);
      }
   }
   return Result;
}

IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
SingularFactorize(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                  VectorBasis, VectorBasis> const& m)
{
   DEBUG_PRECONDITION_EQUAL(m.Basis1(), m.Basis2());

   if (is_regular_basis(m.Basis1()))
      return SingularFactorizeRegular(m);
   // else
   Regularizer R(m.Basis1());
   return UnregularizeBasis12(R, SingularFactorize(RegularizeBasis12(R, m, R)), R);
}

} // namespace Tensor
