// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/bosehubbard-spinless-u1.h"

int main(int argc, char** argv)
{
   if (argc != 7)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: bosehubbard-spinless-u1 <MaxN> <U> <Nwell> <K> <L> <outfile>\n"
		<< "MaxN - maximum number of bosons per site\n"
		<< "U  - on-site coulomb strength\n"
		<< "Nwell - number of wells in the trap, possible is 0,1,2\n"
		<< "K - trap strength\n"
		<< "L  - final number of lattice sites\n"
                << "outfile = file name for output lattice.\n";
      return 1;
   }

   int MaxN = boost::lexical_cast<int>(argv[1]);
   double U = boost::lexical_cast<double>(argv[2]);
   int Nwell = boost::lexical_cast<int>(argv[3]);
   double K = boost::lexical_cast<double>(argv[4]);
   int L = boost::lexical_cast<int>(argv[5]);
   std::string FileName = argv[6];

   TRACE(MaxN)(U)(K)(L);

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
   MPOperator& Hamiltonian = OpList["H"];
   MPOperator& HV = OpList["H_V"];
   MPOperator& HVV = OpList["H_VV"];

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   

   // hopping matrix elements
   for (int i = 1; i < L; ++i)
   {
      Hamiltonian += -1.0 * (prod(BH(i), B(i%L+1), Ident) 
			     + prod(B(i), BH(i%L+1), Ident));

      HV += N(i) * N(i+1);
      if (i < L-1)
	 HVV += N(i) * N(i+2);

      std::cout << "Working.... " << i << "\n";
   }

   for (int i = 1; i <= L; ++i)
   {
      // Interaction, and trap
      Hamiltonian += (U/2.0) * N2(i);

      if (Nwell == 1)
	 Hamiltonian += K * ((i - L/2 - 0.5) * (i - L/2 - 0.5)) * N(i);
      else if (Nwell == 2)
	 Hamiltonian += K * ((i - L/2 - 0.5) * (i - L/2 - 0.5)) * (i - 1) * (i - L) * N(i);
      else if (Nwell != 0) PANIC("Nwell parameter invalid")(Nwell);

      std::cout << "Working.... " << i << "\n";
   }

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(FileName, 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
