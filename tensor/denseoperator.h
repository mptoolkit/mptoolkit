// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/denseoperator.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#if !defined(DENSEOPERATOR_H_DCKJHFAUOIFHIU4YU58943Y9PFWE)
#define DENSEOPERATOR_H_DCKJHFAUOIFHIU4YU58943Y9PFWE

#include "tensor/tensor.h"
#include "tensor/tensorsum.h"
#include "tensor/tensorproduct.h"
#include "basis.h"

classs DenseOperator
{
   public:
      typedef double value_type;

      typedef IrredTensor<LinearAlgebra::Matrix<value_type> > data_type;

      DenseOperator();

      DenseOperator(VectorBasis const& b1, VectorBasis const& b2, QuantumNumber const& q);

      DenseOperator(VectorBasis const& b1, VectorBasis const& b2, data_type const& d);

      VectorBasis const& Basis1() const { return Basis1_; }
      VectorBasis const& Basis2() const { return Basis2_; }

      data_type& data() { return Data_; }
      data_type const& data() const { return Data_; }

      QuantumNumber const& TransformsAs() const { return Data_.TransformsAs(); }

      DenseOperator& operator+=(DenseOperator const& x);
      DenseOperator& operator-=(DenseOperator const& x);

      DenseOperator& operator*=(value_type a);

   private:
      VectorBasis Basis1_, Basis2_;
      IrredTensor<LinearAlgebra::Matrix<value_type> > Data_;
};

bool operator==(DenseOperator const& x, DenseOperator const& y);
bool operator!=(DenseOperator const& x, DenseOperator const& y);

bool equal(DenseOperator const& x, DenseOperator const& y, 
	   double tol = LinearAlgebra::default_tolerance());

DenseOperator operator-(DenseOperator const& x);

DenseOperator operator+(DenseOperator const& x, DenseOperator const& y);
DenseOperator operator-(DenseOperator const& x, DenseOperator const& y);

DenseOperator::value_type trace(DenseOperator const& x);

double norm_2_sq(DenseOperator const& x);
double norm_2(DenseOperator const& x);

DenseOperator operator*(double a, DenseOperator const& x);
DenseOperator operator*(DenseOperator const& x, double a);

std::ostream& operator<<(std::ostream& out, DenseOperator const& x);

double scalar_prod(DenseOperator const& x, DenseOperator const& y);

DenseOperator prod(DenseOperator const& x, DenseOperator const& y, 
		   QuantumNumber const& q = QuantumNumber());

HermitianProxy<DenseOperator const> herm(DenseOperator const& x);

DenseOperator scalar_product(DenseOperator const& x, HermitianProxy<DenseOperator const> y);
DenseOperator scalar_product(HermitianProxy<DenseOperator const> x, DenseOperator const& y);

DenseOperator adjoint(DenseOperator const& x);
DenseOperator inv_adjoint(DenseOperator const& x);

#if 0
DenseOperator tensor_prod(DenseOperator const& x, DenseOperator const& y,
			  VectorProductBasis const& b1,
			  VectorProductBasis const& b2,
			  QuantumNumber const& q);

DenseOperator tensor_prod(DenseOperator const& x, DenseOperator const& y,
			  VectorProductBasis const& b1,
			  VectorProductBasis const& b2);

DenseOperator tensor_prod(DenseOperator const& x, DenseOperator const& y,
			  VectorProductBasis const& b,
			  QuantumNumber const& q);

DenseOperator tensor_prod(DenseOperator const& x, DenseOperator const& y,
			  VectorProductBasis const& b);

DenseOperator tensor_sum(DenseOperator const& x, DenseOperator const& y,
			 VectorSumBasis const& b1, VectorSumBasis const& b2);

DenseOperator tensor_sum(DenseOperator const& x, DenseOperator const& y,
			 VectorSumBasis const& b);

DenseOperator tensor_row_sum(DenseOperator const& x, DenseOperator const& y,
			     VectorSumBasis const& b1);

DenseOperator tensor_col_sum(DenseOperator const& x, DenseOperator const& y,
			     VectorSumBasis const& b2);
#endif

#include "denseoperator.cc"

#endif
