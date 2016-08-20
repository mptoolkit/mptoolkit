// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/ublas-test.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include <iostream>
#include <vector>
#include <boost/timer.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

using namespace std;
using namespace boost::numeric::ublas;

int main(int argc, char** argv)
{
  const size_t nrows = 2560;
  const size_t ncols = 2560;
  const size_t ntrials = 10;
  boost::timer::timer t;

  sparse_matrix<float> m(nrows, ncols);
  sparse_vector<float> v(ncols), r1(nrows);
  std::vector<std::vector<float> > vm(nrows, std::vector<float>(ncols, float(0)));
  std::vector<float> vv(ncols), r2(nrows, 0.0);

  for (float pctNZ = .1; pctNZ < 1.0; pctNZ += .1) {

    // Initialize matrices
    // Values don't matter except for % NZ
    for (size_t i = 0; i < nrows; ++i)
      for (size_t j = 0; j < ncols; ++j)
         if (float(rand()%10)/10.0 < pctNZ) m(i,j) = 1.0;

    // Initialize RHS std::vectors
    // Values don't matter
    for (size_t i = 0; i < ncols; ++i)
      v[i] = vv[i] = i;

    cout << pctNZ << " ";

    // Ublas multiplication
    t.restart();
    for (size_t n = 0; n < ntrials; ++n)
      noalias(r1) = prod(m, v);
    cout << 1000 * (t.elapsed() / ntrials);

    // Straight naive multiplication
    t.restart();
    for (size_t n = 0; n < ntrials; ++n)
      for (size_t i = 0; i < nrows; ++i)
        for (size_t j = 0; j < ncols; ++j)
          r2[i] += vm[i][j] * vv[j];
    cout << " " << 1000 * (t.elapsed() / ntrials) << endl;
  }

  return 0;
}
