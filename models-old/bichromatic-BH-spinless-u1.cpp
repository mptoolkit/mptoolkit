// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/bosehubbard-spinless-u1.h"

namespace Bichromatic
{
  inline double Potential(int j, double r, double phi)
  {
    return 0.5*(1.0 + cos(2.0*r*math_const::pi*j + 2.0*phi));
  }
}

int main(int argc, char** argv)
{
   if (argc != 11)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: bosehubbard-spinless-u1 <MaxN> <U> <V2> <r> <phi> <Nwell> <K> <L> <PBC> <outfile>\n"
		<< "MaxN - maximum number of bosons per site\n"
		<< "U - on-site coulomb strength\n"
		<< "V2 - Perturbing potential\n"
		<< "r - wave-length ratio\n"
		<< "phi - initial shift\n"
		<< "Nwell - number of wells in the trap, possible is 0,1\n"
		<< "K - trap strength\n"
		<< "L - final number of lattice sites\n"
		<< "PBC - adding a hopping link between 1 and L, 0=>OBC, 1=>PBC, -1=>APBC.\n"
                << "outfile = file name for output lattice.\n";
      return 1;
   }

   int MaxN = boost::lexical_cast<int>(argv[1]);
   double U = boost::lexical_cast<double>(argv[2]);
   double V2 = boost::lexical_cast<double>(argv[3]);
   double r = boost::lexical_cast<double>(argv[4]);
   double phi = boost::lexical_cast<double>(argv[5]);
   int Nwell = boost::lexical_cast<int>(argv[6]);
   double K = boost::lexical_cast<double>(argv[7]);
   int L = boost::lexical_cast<int>(argv[8]);
   int PBC = boost::lexical_cast<int>(argv[9]);

   TRACE(MaxN)(U)(V2)(r)(phi)(K)(L)(PBC);

   // Construct the site block
   SiteBlock Site = CreateBoseHubbardSpinlessU1Site(MaxN);

   // construct a lattice of L copies of Site
   Lattice MyLattice = repeat(Site, L);
   MyLattice.fix_coordinates();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   OperatorAtSite<OperatorList const, int> BH(OpList, "BH");
   OperatorAtSite<OperatorList const, int> B(OpList, "B");
   OperatorAtSite<OperatorList const, int> N2(OpList, "N2");
   OperatorAtSite<OperatorList const, int> N(OpList, "N");

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   MPOperator Hamiltonian;
   // hopping matrix elements
   for (int i = 1; i < L; ++i)
   {
      Hamiltonian += -1.0 * (prod(BH(i), B(i%L+1), Ident) 
			     + prod(B(i), BH(i%L+1), Ident));
      std::cout << "Working.... " << i << "\n";
   }
   if (PBC == 1 || PBC == -1)
   {
     Hamiltonian += -1.0 * double(PBC) * (prod(BH(1), B(L), Ident)
					  + prod(B(1), BH(L), Ident));
     std::string s = ""; if(PBC == -1) s = 'A';
     std::cout << "Working.... "<<s<<"PBC link\n";
   }

   for (int i = 1; i <= L; ++i)
   {
      // Interaction, and trap
      Hamiltonian += (U/2.0) * N2(i);

      Hamiltonian += V2 * Bichromatic::Potential(i,r,phi) * N(i);

      if (Nwell == 1)
	 Hamiltonian += K * ((i - L/2 - 0.5) * (i - L/2 - 0.5)) * N(i);
      else if (Nwell == 2)
	 Hamiltonian += K * ((i - L/2 - 0.5) * (i - L/2 - 0.5)) * (i - 1) * (i - L) * N(i);
      else if (Nwell != 0) PANIC("Nwell parameter invalid")(Nwell);

      std::cout << "Working.... " << i << "\n";
   }

   // save the Hamiltonian into the operator list
   OpList["H"] = Hamiltonian;

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(argv[10], 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
