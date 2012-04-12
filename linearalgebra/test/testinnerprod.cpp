
#include "linearalgebra/vector.h"
#include <complex>
using namespace LinearAlgebra;

typedef std::complex<double> complex;

int main()
{
   Vector<complex> v1 = range(0, 10);
   Vector<complex> v2 = range(20,30);

   TRACE(inner_prod(v1, v2));
   TRACE(inner_prod(herm(v1), v2));
}
