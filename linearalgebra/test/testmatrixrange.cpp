// -*- C++ -*- $Id$

#include "linearalgebra/vector.h"
#include "linearalgebra/matrix.h"


using namespace std;
using namespace LinearAlgebra;

template <typename T>
typename ProjectMatrix<Matrix<T>,Range,all_t>::result_type
foo(Matrix<T> const& x)
{
   return x(range(1,3), all);
}

template <typename T>
void test()
{
   Matrix<T> M(4,4,2.0);
   Matrix<T> N(4,4,1.0);

   CHECK_EQUAL(inner_prod(N, N), 16.0);
   CHECK_EQUAL( inner_prod( N(all,range(1,2)), M(all,range(1,2)) ), 8.0);
   CHECK_EQUAL( inner_prod( N(range(1,2), all), M(range(1,2), all) ), 8.0);

   N(range(1,2),range(1,2)) = N(range(1,2),all) * M(all,range(1,2));
   CHECK_EQUAL( N(range(1,2),range(1,2)), Matrix<double>(1,1,8) );
}

int main()
{
   Matrix<double> M(4,4,2);
   foo(M);

   Matrix<double> N(4,4,1);

#if 0
   TRACE(N);
   TRACE(N(all, range(1,2)));
   TRACE(N(range(1,2), range(1,2)));
   TRACE(inner_prod(N, N));
   TRACE( inner_prod( N(all,range(1,2)), M(all,range(1,2)) ) );
   
   TRACE( N(range(1,2),range(1,2)) )( N(range(1,2),all) )( M(all,range(1,2)) );
   N(range(1,2),range(1,2)) = N(range(1,2),all) * M(all,range(1,2));
   TRACE( N(range(1,2),range(1,2)) )( N(range(1,2),all) )( M(all,range(1,2)) );
   CHECK_EQUAL( N(range(1,2),range(1,2)), Matrix<double>(1,1,8) );
#endif

   test<double>();
   test<std::complex<double> >();
}
