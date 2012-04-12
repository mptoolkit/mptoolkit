
#include "denseoperator.h"

std::ostream& operator<<(std::ostream& out, DenseOperator const& x)
{
   out << "DenseOperator:\nBasis1() = " << x.Basis1()
       << "\nBasis2() = " << x.Basis2()
       << "data() =\n" << x.data();
   return out;
}

   
