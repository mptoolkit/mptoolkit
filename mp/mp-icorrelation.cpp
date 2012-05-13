
#include "mps/infinitewavefunction.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include "common/environment.h"
#include "interface/inittemp.h"

#include "tensor/tensor_eigen.h"

#include "models/spin.h"
#include "models/spin-u1.h"
#include "models/spin-su2.h"
#include "models/tj-u1su2.h"
#include "models/tj-u1.h"
#include "models/spinlessfermion-u1.h"
#include "models/kondo-u1su2.h"
#include "models/kondo-u1.h"
#include "models/hubbard-so4.h"
#include "models/bosehubbard-spinless-u1.h"

// For simplicity we don't do cross-correlations here, since that would require
// first calculating the identity operator between the two states.  In fact,
// all cross-correlations are identically zero?

int main(int argc, char** argv)
{
   if (argc < 6 || argc > 6)
   {
      std::cout << "usage: mp-icorrelation <model> <Op1> <Op2> <psi> <length>\n";
      return 1;
   }

   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pvalue_ptr<InfiniteWavefunction> Psi = pheap::OpenPersistent(argv[4], CacheSize, true);
   int Len = boost::lexical_cast<int>(argv[5]);
   std::string Symmetry = argv[1];
   std::string Op1Str = argv[2];
   std::string Op2Str = argv[3];
   int NMax = 3;

   SiteBlock Site;
   if (Symmetry == "sf-u1")
   {
      Site = CreateU1SpinlessFermion();
   }
   else if (Symmetry == "spin")
   {
      Site = CreateSpinSite(0.5);
   }
   else if (Symmetry == "spin1")
   {
      Site = CreateSpinSite(1.0);
   }
   else if (Symmetry == "spin-su2")
   {
      Site = CreateSU2SpinSite(0.5);
   }
   else if (Symmetry == "spin-u1")
   {
      Site = CreateU1SpinSite(0.5);
   }
   else if (Symmetry == "tj-u1")
   {
      Site = CreateU1tJSite();
   }
   else if (Symmetry == "tj-u1su2")
   {
      Site = CreateU1SU2tJSite();
   }
   else if (Symmetry == "klm-u1su2")
   {
      Site = CreateU1SU2KondoSite();
   }
   else if (Symmetry == "klm-u1")
   {
      Site = CreateU1KondoSite();
   }
   else if (Symmetry == "hubbard-so4")
   {
      Site = CreateSO4HubbardSiteA();
   }
   else if (Symmetry == "bh-u1")
   {
      Site = CreateBoseHubbardSpinlessU1Site(NMax);
   }
   else
   {
      PANIC("mp-icorrelation: fatal: model parameter should be one of tj-u1, tj-u1su2, sf-u1, klm-u1su2.");
   }

   SimpleOperator Op1 = Site[Op1Str];
   SimpleOperator Op2 = Site[Op2Str];

   LinearWavefunction p = get_orthogonal_wavefunction(*Psi);
   MatrixOperator Rho = scalar_prod(Psi->C_right, herm(Psi->C_right));
   //   MatrixOperator Rho = scalar_prod(herm(Psi->C_right), Psi->C_right);
   //   MatrixOperator Rho = scalar_prod(herm(Psi->C_old), Psi->C_old);
   //   Rho = delta_shift(Rho, adjoint(Psi->QShift));

   std::cout.precision(14);

   QuantumNumber Ident(Psi->GetSymmetryList());

   // for each site in the unit cell, make the right hand operator
   std::list<MatrixOperator> F;
   LinearWavefunction::const_iterator I = p.end();
   while (I != p.begin())
   {
      --I;
      // add identity operator to each existing operator
      for (std::list<MatrixOperator>::iterator J = F.begin(); J != F.end(); ++J)
      {
         *J = operator_prod(*I, *J, herm(*I));
      }

      // next site
      F.push_back(operator_prod(Op2, *I, herm(*I)));
   }

   // now the left hand operators.   We also do the correlators for the first unit cell
   typedef std::list<std::pair<int, MatrixOperator> > EType;
   typedef std::list<std::pair<std::pair<int, int>, MatrixOperator> > CorrType;
   std::list<std::pair<int, MatrixOperator> > E;
   std::list<std::pair<std::pair<int, int>, MatrixOperator> > Corr;  // first unit cell correlators
   I = p.begin();
   int pLoc = 1;
   while (I != p.end())
   {
      // add identity operator to each existing first unit cell operator
      for (CorrType::iterator J = Corr.begin(); J != Corr.end(); ++J)
      {
         J->second = operator_prod(herm(*I), J->second, *I);
      }
      // add the next site to the first unit cell correlator
      for (EType::const_iterator J = E.begin(); J != E.end(); ++J)
      {
         Corr.push_back(std::make_pair(std::make_pair(J->first, pLoc),
				       operator_prod(herm(Op2), herm(*I), J->second, *I, Ident)));
      }
      // add identity operator to each existing operator
      for (EType::iterator J = E.begin(); J != E.end(); ++J)
      {
         J->second = operator_prod(herm(*I), J->second, *I);
      }
      // next site
      E.push_back(std::make_pair(pLoc, operator_prod(herm(Op1), herm(*I), *I)));
      ++I;
      ++pLoc;
   }

   // output the first unit cell
   for (CorrType::const_iterator J = Corr.begin(); J != Corr.end(); ++J)
   {
      std::complex<double> x = inner_prod(J->second, Rho);
      std::cout << J->first.first << "    " << J->first.second << "   " 
		<< x.real() << "   " << x.imag() << '\n';
   }

   int NCell = 1; // number of cells
   while (NCell * int(p.size()) < Len)
   {
      for (EType::iterator J = E.begin(); J != E.end(); ++J)
      {
	 J->second = delta_shift(J->second, Psi->QShift);
      }

      Corr.clear();
      for (LinearWavefunction::const_iterator I = p.begin(); I != p.end(); ++I)
      {
         // add identity site to operators already in Corr
         for (CorrType::iterator J = Corr.begin(); J != Corr.end(); ++J)
         {
            J->second = operator_prod(herm(*I), J->second, *I);
         }
         // add next site
         for (EType::const_iterator J = E.begin(); J != E.end(); ++J)
         {
            Corr.push_back(std::make_pair(std::make_pair(J->first, pLoc), 
					  operator_prod(herm(Op2), herm(*I), J->second, *I, Ident)));
         }
         // add identity operator to E
         for (EType::iterator J = E.begin(); J != E.end(); ++J)
         {
            J->second = operator_prod(herm(*I), J->second, *I);
         }
	 ++pLoc;
      }

      // output
      for (CorrType::const_iterator J = Corr.begin(); J != Corr.end(); ++J)
      {
	 std::complex<double> x = inner_prod(J->second, Rho);
	 std::cout << J->first.first << "    " << J->first.second << "   " 
		   << x.real() << "   " << x.imag() << '\n';
      }
      ++NCell;
   }

   pheap::Shutdown();
}
