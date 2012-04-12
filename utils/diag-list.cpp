// -*- C++ -*- $Id$

#include "linearalgebra/matrix.h"
#include "linearalgebra/eigen.h"
#include <complex>

using namespace LinearAlgebra;

int main()
{
   std::vector<std::complex<double> > v;
   std::complex<double> z;
   while (std::cin >> z)
      v.push_back(z);

   int n = int(std::sqrt(v.size()));
   Matrix<std::complex<double> > M(n, n);
   std::vector<std::complex<double> >::const_iterator vI = v.begin();
   for (iterator<Matrix<std::complex<double> > >::type I = iterate(M); I; ++I)
      for (inner_iterator<Matrix<std::complex<double> > >::type J = iterate(I); J; ++J)
	 *J = *vI++;

   std::cout.precision(12);
   std::cout << "Size: " << v.size() << '\n';
   std::cout << "Eigenvalues:\n";
   std::cout << EigenvaluesHermitian(M) << '\n';
   DiagonalizeHermitian(M);
   std::cout << "Eigenvectors:\n";
   std::cout << M << std::endl;
}
