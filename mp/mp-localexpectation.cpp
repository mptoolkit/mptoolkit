// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/centerwavefunction.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include <iomanip>
#include "mp/copyright.h"
#include "common/environment.h"
#include "common/proccontrol.h"

typedef std::complex<double> complex;

int main(int argc, char** argv)
{
   if (argc != 6)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: mp-localexpectation <lattice> <psi> <operator> <first> <last>\n";
      std::cerr << "This uses directly local operators, not lattice operators!\n";
      return 1;
   }

   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pvalue_ptr<OperatorList> System = pheap::OpenPersistent(argv[1], CacheSize, true);
   pvalue_ptr<MPWavefunction> Psi1 = pheap::ImportHeap(argv[2]);

   std::string Op1 = argv[3];
   int FirstSite = boost::lexical_cast<int>(argv[4]);
   int LastSite = boost::lexical_cast<int>(argv[5]);

   Lattice Lat = System->GetLattice();

   CenterWavefunction Psi = *Psi1;
   Psi1 = pvalue_ptr<MPWavefunction>();

   for (int i = 1; i < FirstSite; ++i)
      Psi.RotateRight();

   std::cout.precision(12);

   // do up to LastSite-1 by applying the operator to the left
   for (int i = FirstSite; i < LastSite; ++i)
   {
      SiteBlock::const_iterator I = Lat[i].find(Op1);
      if (I != Lat[i].end())
      {
         SimpleOperator MyOp = I->second;
         MatrixOperator LeftIdentity = MatrixOperator::make_identity(Psi.Left().Basis1());
         MatrixOperator x = operator_prod(herm(MyOp), herm(Psi.Left()), LeftIdentity, Psi.Left());
         complex Res = inner_prod(Psi.Center(), x * Psi.Center());
         std::cout << std::setw(5) << Lat.coordinate_at_site(i) << "   " 
                   << std::setw(18) << Res.real() << "   " 
                   << std::setw(18) << Res.imag() << '\n';
      }
      if (i != LastSite-1)
         Psi.RotateRight();
   }

   // The last site might be at the edge of the lattice, so act on the right side instead
   SiteBlock::const_iterator I = Lat[LastSite].find(Op1);
   if (I != Lat[LastSite].end())
   {
      SimpleOperator MyOp = I->second;
      MatrixOperator Ident = MatrixOperator::make_identity(Psi.Right().Basis2());
      MatrixOperator x = operator_prod(MyOp, Psi.Right(), Ident, herm(Psi.Right()));
      complex Res = inner_prod(Psi.Center() * x, Psi.Center());
      std::cout << std::setw(5) << Lat.coordinate_at_site(LastSite) << "   " 
                << std::setw(18) << Res.real() << "   " 
                << std::setw(18) << Res.imag() << '\n';
   }

   pheap::Shutdown();
}
