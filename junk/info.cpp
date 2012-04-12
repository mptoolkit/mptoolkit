// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/u1.h"
#include "pheap/pheap.h"
#include "matrixproduct/copyright.h"
#include <iostream>

void ShowWavefunc(MPWavefunction Psi, std::ostream& out, int Max = -1)
{
   out << "Symmetry list is " << Psi.GetSymmetryList() << '\n';
   out << "State transforms as " << Psi.TransformsAs() << '\n';
   out << "Number of sites = " << Psi.Size() << '\n';
   //   out << "Left-most basis:\n" << Psi.Left().Basis1();
   while (Psi.RightSize() > 1)
   {
      //      out << "Matrices at " << Psi.LeftSize() << '\n';
      //      out << Psi.Left().Basis2() << '\n' << Psi.Right().Basis1() << '\n';
      //      out << mp_prod_left(Psi.Left(), Psi.Left()) << '\n';
      //      out << "trace of norm matrix at " << Psi.LeftSize() << " is "
      //	 out << trace(mp_prod_left(Psi.Left(), Psi.Left())) << '\n';
      DensityMatrix<BlockOperator> DM(prod(Psi.Center(), adjoint(Psi.Center()), 
					   Psi.Center().TransformsAs()));
      out << "\nReduced density matrix at partition (" 
	  << Psi.LeftSize() << "," << Psi.RightSize() << ") :\n";
      DM.DensityMatrixReport(out, Max);
      Psi.RotateRight();
   }

   DensityMatrix<BlockOperator> DM(prod(Psi.Center(), adjoint(Psi.Center()), 
					Psi.Center().TransformsAs()));
   out << "\nReduced density matrix at partition (" 
       << Psi.LeftSize() << "," << Psi.RightSize() << ") :\n";
   DM.DensityMatrixReport(out, Max);
}


int main(int argc, char** argv)
{
   if (argc < 2 || argc > 3)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: mp-info <psi> [max-eigenvalues]\n";
      return 1;
   }

   int Max = -1;
   if (argc == 3) Max = boost::lexical_cast<int>(argv[2]);

   pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(argv[1], 655360, true);

   std::cout.precision(13);
   ShowWavefunc(*Psi, std::cout, Max);

   pheap::Shutdown();
}
