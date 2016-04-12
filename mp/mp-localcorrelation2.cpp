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
#include "interface/inittemp.h"
#include "common/proccontrol.h"

int main(int argc, char** argv)
{
   if (argc < 7 || argc > 8)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-localcorrelation <lattice> <psi> <operator1> <first> "
	 "<operator2> <last> [<psi2>]\n";
      std::cerr << "This uses directly local operators, not lattice operators!\n";
      return 1;
   }

   mp_pheap::InitializeTempPHeap();

   pvalue_ptr<OperatorList> System = pheap::ImportHeap(argv[1]);
   pvalue_ptr<MPWavefunction> Psi1Ptr = pheap::ImportHeap(argv[2]);
   pvalue_ptr<MPWavefunction> Psi2Ptr = argc > 7 ? pheap::ImportHeap(argv[2]) : Psi1Ptr;

   std::string Op1 = argv[3];
   int FirstSite = boost::lexical_cast<int>(argv[4]);
   std::string Op2 = argv[5];
   int LastSite = boost::lexical_cast<int>(argv[6]);

   Lattice Lat = System->GetLattice();

   CenterWavefunction Psi1 = *Psi1Ptr;
   CenterWavefunction Psi2 = *Psi2Ptr;
   Psi1Ptr = pvalue_ptr<MPWavefunction>();
   Psi1Ptr = pvalue_ptr<MPWavefunction>();
   
   for (int i = 1; i < FirstSite; ++i)
      Psi.RotateRight();

   std::cout.precision(12);

   // This holds the E matrices for the left system
   typedef std::map<int, pvalue_handle<MatrixOperator> > OpMapType;
   OpMapType OpMap;
   for (int i = FirstSite; i < LastSite; ++i)
   {
      // Update all the existing E matrices
      for (OpMapType::iterator mI = OpMap.begin(); mI != OpMap.end(); ++mI)
      {
         mI->second = new MatrixOperator(operator_prod(herm(Psi.Left()), 
                                                       *mI->second.load(), 
                                                       Psi.Left()));
      }

      // Add another E matrix for the new site, only if the left operator exists in the block
      SiteBlock::const_iterator I = Lat[i].find(Op1);
      if (I != Lat[i].end())
      {
         SimpleOperator MyOp = I->second;
         MatrixOperator LeftIdentity = MatrixOperator::make_identity(Psi.Left().Basis1());
         OpMap[i] = new MatrixOperator(operator_prod(herm(MyOp), 
                                                     herm(Psi.Left()), 
                                                     LeftIdentity, 
                                                     Psi.Left()));
      }

      // For the right operator, construct the F matrix and the expectation value
      I = Lat[i+1].find(Op2);
      if (I != Lat[i+1].end())
      {
         SimpleOperator MyOp = I->second;
         MatrixOperator Ident = MatrixOperator::make_identity(Psi.Right().Basis2());
         MatrixOperator F = operator_prod(MyOp, Psi.Right(), Ident, herm(Psi.Right()));
         for (OpMapType::iterator mI = OpMap.begin(); mI != OpMap.end(); ++mI)
         {
            std::complex<double> Res = inner_prod(Psi.Center(), 
                                                  triple_prod(*mI->second.load(), 
                                                              Psi.Center(), 
                                                              herm(F)));
            std::cout << std::setw(5) << Lat.coordinate_at_site(mI->first) << "   " 
                      << std::setw(5) << Lat.coordinate_at_site(i+1) << "   "
                      << std::setw(18) << Res.real() << "   " 
                      << std::setw(18) << Res.imag() << '\n';
         }
      }

      if (i != LastSite-1)
         Psi.RotateRight();
   }

   pheap::Shutdown();
}
