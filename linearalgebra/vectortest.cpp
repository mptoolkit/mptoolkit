// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/vectortest.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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

#include "vector.h"
#include "vectoroperations.h"
#include <iostream>
#include <iomanip>
#include <typeinfo>

using namespace std;

using namespace LinearAlgebra;

int main()
{
  Vector<double> MyVec(10, 0);

  for (int i = 0; i < MyVec.size(); ++i)
    MyVec[i] = i;

  cout << "MyVec is " << MyVec << '\n';

  VectorRef<double> RefToMyVec(MyVec);

  cout << "RefToMyVec is " << RefToMyVec << '\n';

  RefToMyVec[3] = 10;
  cout << "After RefToMyVec[3] = 10, MyVec is " << MyVec << '\n';

  VectorConstRef<double> ConstRefToMyVec(RefToMyVec);
  cout << "ConstRefToMyVec is " << ConstRefToMyVec << '\n';

  //  VectorSlice<double> MySlice(MakeSlice(MyVec, Slice(1, 4, 2)));
  VectorSlice<double> MySlice(MyVec.slice(Slice(1, 4, 2)));
  std::cout << "MySlice <- Slice(1,4,2) is " << MySlice << '\n';

  MySlice.fill(1);
  std::cout << "After MySlice.fill(1), MyVec is " << MyVec << '\n';

  MySlice.range(Range(1, 3)).fill(-1);
  std::cout << "After MySlice.range(Range(1, 3)).fill(-1), MyVec is " << MyVec << '\n';

  VectorConstSlice<double> MyConstSlice(MySlice);

  //  MySlice += 2.1 * MySlice;

  MySlice.range(Range(0,2)) = 3.4 * MySlice.range(Range(2,4));
  std::cout << "After MySlice.range(Range(0,2)) = 3.4 * MySlice.range(2,4), MyVec is "
            << MyVec << '\n';

  //  BLAS::daxpy(MySlice.size(), 2.1, MySlice.begin().base(), MySlice.stride(), MySlice.begin().base(), MySlice.stride());

  //  ops::fast_add_scaled(2.1,MyConstSlice.begin(), MyConstSlice.end(), MySlice.begin());

  std::vector<int> i;
  i.push_back(2);
  i.push_back(5);
  i.push_back(7);

  Index I(i.begin(), i.end());

  std::cout << I.begin() << ' ' << I.end() << '\n';

  std::cout << "Index is " << I << '\n';
  std::cout << "MyVec is " << MyVec << '\n';
  std::cout << "Index of MyVec is " << index(MyVec, I) << '\n';
  std::cout << "Range of index is " << range(index(MyVec, I), Range(1, 3)) << '\n';
  std::cout << "Range of index alternate is " << VectorIndex<double>(index(MyVec, I), Range(1,3)) << '\n';
  std::cout << "Index of range is " << index(range(MyVec, Range(1, 9)), I) << '\n';

  std::cout << (5.5 * MySlice) << '\n';
  std::cout << norm_2(5.5 * MySlice) << '\n';

  std::cout << typeid(5.5 * MySlice).name() << '\n';
  std::cout << typeid(3.4 * (5.5 * (MySlice))).name() << '\n';
  std::cout << typeid(3.4 * (5.5 * (-MySlice))).name() << '\n';

  std::cout << ((-100) / 1U) << '\n';
}
