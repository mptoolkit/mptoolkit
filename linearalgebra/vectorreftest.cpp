// -*- C++ -*- $Id$

#include "vector.h"
#include "vectoroperations.h"
#include <iostream>
#include <iomanip>
#include <typeinfo>

using namespace std;

using namespace LinearAlgebra;

int main()
{
  Vector<double> V1(10, 0);

  V1 = Vector<double>(Range(10, 100));

  VectorRef<double> V1Ref(V1);
  TRACE(V1Ref);

  VectorSlice<double> V1Slice(V1, Slice(0, 2, 4));
  TRACE(V1Slice);

  Vector<double> V2 = V1;

  VectorConstRef<double> V2Ref = V2;

  V1[0] = -10;
  V1.resize(5);

  TRACE(V1);
  TRACE(V2);
  TRACE(V1Ref);
  TRACE(V2Ref);
  TRACE(V1Slice);
  TRACE(range(V2, Range(10, 20)) + range(V2, Range(20, 30)));
  TRACE(norm_2(range(V2, Range(10, 20))));
  TRACE(norm_2(range(V2, Range(10, 20)) - range(V2, Range(20, 30))));

  double x = 0;
  for (int i = 20; i < 30; ++i)
  {
    x += i*i;
  }
  std::cout << x << ' ' << std::sqrt(x) << '\n';
}
