// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// quantumnumbers/dihedral.cpp
//
// Copyright (C) 2016 Ian McCulloch <ian@qusim.net>
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

// Dihedral 3j, CG and 6j symbols
// From: P H Butler and M F Reid, J. Phys. A: Math. Gen., Vol. 12, No. 10, 1979.
// http://iopscience.iop.org/0305-4470/12/10/012

//
// D_infinity
//

double
DihedralInf6j(signed_halfint j1, signed_halfint j2, signed_halfint j3,
              signed_halfint j4, signed_halfint j5, signed_halfint j6)
{

}

double
DihedralInf9j(half_int j11, half_int j12, half_int j13,
              half_int j21, half_int j22, half_int j23,
              half_int j31, half_int j32, half_int j33)
{
}

double
DihedralInfClebschGordan(half_int j1, half_int m1,
                         half_int j2, half_int m2,
                         half_int j,  half_int m)
{
}

//
// Dihedral groups D_n for finite n
//

double
Dihedral6j(int n,
           signed_halfint j1, signed_halfint j2, signed_halfint j3,
           signed_halfint j4, signed_halfint j5, signed_halfint j6)
{

}

double
Dihedral9j(int n,
           half_int j11, half_int j12, half_int j13,
           half_int j21, half_int j22, half_int j23,
           half_int j31, half_int j32, half_int j33)
{
}

double
DihedralClebschGordan(int n,
                      half_int j1, half_int m1,
                      half_int j2, half_int m2,
                      half_int j,  half_int m)
{
}
