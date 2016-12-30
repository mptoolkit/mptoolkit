// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// quantumnumbers/coupling.h
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

#if !defined(COUPLING_H_SDFJH3T89EG89J34JIOJ8Q89C89J)
#define COUPLING_H_SDFJH3T89EG89J34JIOJ8Q89C89J

#include "common/halfint.h"
#include "common/gmprational.h"
#include <utility>

//
// Coupling6j
//
// returns the 6j coefficient { j1 j2 j3 }
//                            { j4 j5 j6 }
//

double Coupling6j(half_int j1, half_int j2, half_int j3, half_int j4, half_int j5, half_int j6);

double Racah(half_int ja, half_int jb, half_int jc, half_int jd, half_int je, half_int jf);

// low-level implementation of the Racah symbol, without caching.
double Racah_NoCache(half_int ja, half_int jb, half_int jc, half_int jd, half_int je, half_int jf);


double Coupling9j(half_int j11, half_int j12, half_int j13,
                  half_int j21, half_int j22, half_int j23,
                  half_int j31, half_int j32, half_int j33);

double ClebschGordan(half_int j1, half_int m1,
                     half_int j2, half_int m2,
                     half_int j,  half_int m);

// returns the clebsch-Gordan coefficent in exact form,
// as first * sqrt(second)
std::pair<gmp::rational, gmp::rational>
 ClebschGordanSquared(half_int j1, half_int m1,
                      half_int j2, half_int m2,
                      half_int j,  half_int m);

#endif
