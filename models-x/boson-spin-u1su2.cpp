// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/bosehubbard-u1su2.h"

int main(int argc, char** argv)
{
   if (argc != 9)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: u1su2bosehubbard <MaxN> <MaxS> <U> <K> <J1> <J2> <L> <outfile>\n"
		<< "MaxN - maximum number of bosons per site\n"
		<< "MaxS - maximum spin on each site\n"
		<< "U  - on-site coulomb strength\n"
		<< "K  - on-site spin^2 strength\n"
		<< "J1 - nearest-neighbor bilinear spin interaction\n"
		<< "J2 - nearest-neighbor biquadratic spin interaction\n"
		<< "L  - final number of lattice sites\n"
                << "outfile = file name for output lattice.\n";
      return 1;
   }

   int MaxN = boost::lexical_cast<int>(argv[1]);
   int MaxS = boost::lexical_cast<int>(argv[2]);
   double U = boost::lexical_cast<double>(argv[3]);
   double K = boost::lexical_cast<double>(argv[4]);
   double J1 = boost::lexical_cast<double>(argv[5]);
   double J2 = boost::lexical_cast<double>(argv[6]);
   int L = boost::lexical_cast<int>(argv[7]);

   TRACE(MaxN)(MaxS)(U)(K)(L);

   // Construct the site block
   SiteBlock Site = CreateU1SU2BoseHubbardSite(MaxN, MaxS);

   // construct a lattice of L copies of Site
   Lattice MyLattice = repeat(Site, L);
   MyLattice.fix_coordinates();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   OperatorAtSite<OperatorList const, int> BH(OpList, "BH");
   OperatorAtSite<OperatorList const, int> B(OpList, "B");
   OperatorAtSite<OperatorList const, int> N2(OpList, "N2");
   OperatorAtSite<OperatorList const, int> Q(OpList, "Q");
   OperatorAtSite<OperatorList const, int> S(OpList, "S");
   OperatorAtSite<OperatorList const, int> S2(OpList, "S2");
   OperatorAtSite<OperatorList, int> Bond(OpList, "Bond");
   OperatorAtSite<OperatorList, int> SS(OpList, "SS");

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   MPOperator Hamiltonian;
   // hopping matrix elements
   for (int i = 1; i < L; ++i)
   {
      MPOperator Term
         = -sqrt(3.0) * (prod(BH(i), B(i%L+1), Ident) + prod(B(i), BH(i%L+1), Ident));
      MPOperator SSn = -sqrt(3.0) * prod(S(i), S(i%L+1), Ident);
      Term += J1 * SSn + J2 * prod(SSn, SSn, Ident);
      Term += U/2.0*(N2(i)+N2(i%L+1));
      Term += K/2.0*(S2(i)+S2(i%L+1));

      if (i == 1) Term += U/2.0*N2(1) + K/2.0*S2(1);
      if (i == L-1) Term += U/2.0*N2(L) + K/2.0*S2(L);

      Bond(i) = Term;
      Hamiltonian = Hamiltonian + Term;
      SS(i) = SSn;

      std::cout << "Working.... " << i << "\n";
   }

   // save the Hamiltonian into the operator list
   OpList["H"] = Hamiltonian;

   // Make some other useful(?) operators
   MPOperator TotalQ, TotalQ2, DimerOrder, DimerOrder2;
   for (int i = 1; i <= L; ++i)
   {
      TotalQ += Q(i);
      if (i != L) 
	 DimerOrder += minus1pow(i) * (-sqrt(3.0)) * prod(S(i), S(i%L+1), Ident);
      std::cout << "Working.... " << i << "\n";
   }
   DimerOrder2 = prod(DimerOrder, DimerOrder, Ident);
   TotalQ2 = sqrt(5.0) * prod(TotalQ, TotalQ, Ident);

   OpList["Q"] = TotalQ;
   OpList["Q2"] = TotalQ2;
   OpList["Dimer"] = DimerOrder;
   OpList["Dimer2"] = DimerOrder2;

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(argv[8], 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
