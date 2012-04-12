// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/spin-u1.h"

typedef std::complex<double> complex;

int main(int argc, char** argv)
{
   if (argc != 7)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: spinchain-zigzag-u1 <L> <EdgeSpin> <BulkSpin> <J1> <J2> <outfile>\n"
                << "L = number of lattice sites\n"
                << "EdgeSpin = spin of the left and right extreme sites\n"
                << "BulkSpin = spin away from the boundary\n"
                << "J1 = nearest neighbor bilinear coupling\n"
                << "J2 = next-nearest neighbor bilinear coupling\n"
                << "outfile = file name for output lattice.\n";
      return 1;
   }

   int L = boost::lexical_cast<int>(argv[1]);
   half_int Edge = boost::lexical_cast<half_int>(argv[2]);
   half_int Bulk = boost::lexical_cast<half_int>(argv[3]);
   double J1 = boost::lexical_cast<double>(argv[4]);
   double J2 = boost::lexical_cast<double>(argv[5]);
   std::string OutName = argv[6];

   TRACE(L)(Edge)(Bulk)(J1)(J2);
   CHECK(L >= 2);
   
   // Construct the site block
   SiteBlock EdgeSite = CreateU1SpinSite(Edge);
   SiteBlock BulkSite = CreateU1SpinSite(Bulk);

   // construct a lattice of L copies of Site
   Lattice MyLattice = join(EdgeSite, repeat(BulkSite, L-2), EdgeSite);
   MyLattice.fix_coordinates();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   OperatorAtSite<OperatorList const, int> Sp(OpList, "Sp");
   OperatorAtSite<OperatorList const, int> Sm(OpList, "Sm");
   OperatorAtSite<OperatorList const, int> Sz(OpList, "Sz");
   OperatorAtSite<OperatorList, int> Kappaz(OpList, "Kappaz");
   OperatorAtSite<OperatorList, int> Kappap(OpList, "Kappap");
   OperatorAtSite<OperatorList, int> Kappam(OpList, "Kappam");
   MPOperator& Hamiltonian = OpList["H"];

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number

   Hamiltonian = J1 * (0.5 * (CreateRepeatedOperator(MyLattice, "Sp", "Sm")
			      + CreateRepeatedOperator(MyLattice, "Sm", "Sp"))
		       + CreateRepeatedOperator(MyLattice, "Sz", "Sz"))
      + J2 * (0.5 * (CreateRepeatedOperator(MyLattice, "Sp", "I", "Sm")
		     + CreateRepeatedOperator(MyLattice, "Sm", "I", "Sp"))
	      + CreateRepeatedOperator(MyLattice, "Sz", "I", "Sz"));

   std::cout << "Constructing Kappa operators...\n";
   for (int i = 1; i < L; ++i)
   {
      Kappaz[i] = complex(0.0, 0.5) * (Sp(i)*Sm(i+1) - Sm(i)*Sp(i+1));
      Kappap[i] = complex(0.0, 1.0) * (Sz(i)*Sp(i+1) - Sp(i)*Sz(i+1));
      Kappam[i] = complex(0.0, 1.0) * (Sm(i)*Sz(i+1) - Sz(i)*Sm(i+1));
   }

   std::cout << "done." << std::endl;

   // make a copy of OpList that exists on the persistent heap
   pheap::Initialize(OutName, 1, 65536, 655360);
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);
   pheap::ShutdownPersistent(OList);
}
