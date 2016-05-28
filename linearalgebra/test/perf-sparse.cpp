// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/perf-sparse.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include <iostream> 
#include <vector>
#include <boost/timer.hpp>
#include "linearalgebra/sparsematrix.h"
#include "linearalgebra/vector.h"

#include "common/blas2f.h"

using namespace std;
using namespace LinearAlgebra;

int main(int argc, char** argv)
{
  const size_t nrows = 2560;
  const size_t ncols = 2560;
  const size_t ntrials = 100;
  boost::timer::timer t;

  typedef SparseMatrix<float> Mat;
  typedef Vector<float> Vec;
  SparseMatrix<float> m(nrows, ncols);
  Vector<float> v(ncols), r1(nrows);
  std::vector<std::vector<float> > vm(nrows, std::vector<float>(ncols, float(0)));
  std::vector<float> vv(ncols), r2(nrows, 0.0);

  Matrix<float> m3(nrows, ncols, 0);
  Matrix<float> r3(1, nrows);
  Matrix<float> v3(ncols, 1, 0.0);
  
  for (float pctNZ = .1; pctNZ < 1.0; pctNZ += .1) {

    // Initialize matrices
    // Values don't matter except for % NZ
    for (size_t i = 0; i < nrows; ++i) 
      for (size_t j = 0; j < ncols; ++j) 
         if (float(rand()%10)/10.0 < (pctNZ * 0.01)) m(i,j) = 1.0;
    
    // Initialize RHS std::vectors
    // Values don't matter
    for (size_t i = 0; i < ncols; ++i) 
      v[i] = vv[i] = i;
  
    cout << pctNZ << " ";

    // library multiplication
    t.restart();
    for (size_t n = 0; n < ntrials; ++n)
    {
       zero_all(r1);
       for (const_iterator<Mat>::type i = iterate(m); i; ++i)
       {
          for (const_inner_iterator<Mat>::type j = iterate(i); j; ++j)
          {
             r1[j.index1()] += (*j) * v[j.index2()];
          }
       }
    }
    cout << 1000 * (t.elapsed() / ntrials);
    
    // Straight naive multiplication
    t.restart();
    for (size_t n = 0; n < ntrials; ++n)
      for (size_t i = 0; i < nrows; ++i)
	for (size_t j = 0; j < ncols; ++j)
	  r2[i] += vm[i][j] * vv[j];
    cout << " " << 1000 * (t.elapsed() / ntrials);

  // dense
    t.restart();
    for (size_t n = 0; n < ntrials; ++n)
    {
       //       zero_all(r3);
       BLAS::sgemv('N', ncols, nrows, 1.0, data(m3), ncols, data(v3), 1, 0.0, data(r3), 1);
       //       r3 = m3 * v3;
    }
    cout << " " << 1000 * (t.elapsed() / ntrials);

    cout << endl;
  }

  return 0;
}
