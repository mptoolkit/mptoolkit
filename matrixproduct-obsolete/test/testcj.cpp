// -*- C++ -*- $Id$

#include "linearalgebra/matrix.h"
#include "linearalgebra/vector.h"
#include "matrixproduct/conjugategradient.h"
#include "linearalgebra/matrix_utility.h"
#include "linearalgebra/vector_utility.h"

typedef std::complex<double> complex;
typedef LinearAlgebra::Matrix<complex> Matrix;
typedef LinearAlgebra::Vector<complex> Vector;

struct DoMultiply
{
   DoMultiply(Matrix const& M) : M_(M) {}

   Vector operator()(Vector const& v) const
   {
      Vector Result(size1(M_), 0.0);
      for (int i = 0; i < size1(M_); ++i)
      {
	 for (int j = 0; j < size2(M_); ++j)
	 {
	    Result[i] += M_(i,j) * v[j];
	 }
      }
      return Result;
   }

   Matrix M_;
};

int main()
{
   typedef complex num;

   Matrix M = LinearAlgebra::random_matrix<num>(10,10);

   Matrix P = M + herm(M);

   P = P * herm(P);

   Vector v = LinearAlgebra::random_vector<num>(10);

   int MaxIter = 50;
   double Tol = 1E-10;

   Vector Rhs = v;

   TRACE(Rhs);

   ConjugateGradient(v, DoMultiply(P), Rhs, MaxIter, Tol);
   TRACE(MaxIter)(Tol)(v)(Rhs)(DoMultiply(P)(v))(norm_2(Rhs - DoMultiply(P)(v)));

}
