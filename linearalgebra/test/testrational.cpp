// -*- C++ -*- $Id$

#include "linearalgebra/vector.h"
#include "linearalgebra/matrix.h"
#include "linearalgebra/scalar.h"
#include "linearalgebra/gmpalgebra.h"

using namespace LinearAlgebra;

int main()
{
   Vector<gmp::rational> Vec1(10, gmp::rational(0.1));
   TRACE(Vec1);

   Vector<gmp::rational> Vec2(10, gmp::rational(2.0));

   Vector<gmp::rational> Vec3 = Vec1 + Vec2;

   TRACE(norm_frob_sq(Vec3));

   TRACE(Vec3);

   TRACE(inner_prod(Vec1, Vec2));
}
