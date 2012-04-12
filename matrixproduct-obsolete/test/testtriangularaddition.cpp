
#include "matrixproduct/triangularoperator.h"
#include "models/spin-u1.h"

int main()
{
   SiteBlock Site = CreateU1SpinSite(0.5);

   MpOpTriangular SpSm = TriangularTwoSite(Site["Sp"], Site["Sm"]);
   MpOpTriangular SmSp = TriangularTwoSite(Site["Sm"], Site["Sp"]);
   MpOpTriangular SzSz = TriangularTwoSite(Site["Sz"], Site["Sz"]);

   MpOpTriangular SS = 0.5*(SpSm + SmSp) + SzSz;

   TRACE(SpSm.data())(SS.data());
}
