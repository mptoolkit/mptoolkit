// -*- C++ -*- $Id$

#include "linearalgebra/vector.h"
#include "linearalgebra/matrix.h"
#include "linearalgebra/matrix_utility.h"


using namespace std;
using namespace LinearAlgebra;


int main()
{
   typedef Matrix<double,ColMajor> T;
   
   const int D = 20;
   T M = random_matrix<T::value_type>(D,D);
   
   CHECK_EQUAL( Vector<int>(range(3,15))[range(2,7)], range(5,10) );
   CHECK_EQUAL( M( range(3,15), all)( range(2,7), all)( all, range(2,7) ),
                M( range(5,10), range(2,7) ) );
   
   CHECK_EQUAL( Vector<int>(slice(0,10,2))[slice(1,3,3)], slice(2,3,6) );
   CHECK_EQUAL( M( slice(0,10,2), all)( slice(1,3,3), all)( all, range(2,7) ),
                M( slice(2,3,6), range(2,7) ) );
   
   CHECK_EQUAL( Vector<int>(slice(3,8,2))[range(2,7)][slice(1,2,3)], slice(9,2,6) );
   CHECK_EQUAL( M( slice(3,8,2), all)( range(2,7), all)( slice(1,2,3), all),
                M( slice(9,2,6), all ) );
		
   CHECK_EQUAL( M( slice(3,8,2), slice(0,10,2))( range(2,7), slice(1,3,3))( slice(1,2,3), all),
                M( slice(9,2,6), slice(2,3,6) ) );
}
