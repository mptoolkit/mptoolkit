// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/spin-su2.h"

int main(int argc, char** argv)
{
   if (argc != 7)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: mp-su2spinchain <L> <EdgeSpin> <BulkSpin> <J1> <J2> <outfile>\n"
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

   TRACE(L)(Edge)(Bulk)(J1)(J2);
   CHECK(L >= 2);

   // initialize the output file
   pheap::Initialize(argv[6], 1, 65536, 655360);

   // Construct the site block
   SiteBlock EdgeSite = CreateSU2SpinSite(Edge);
   SiteBlock BulkSite = CreateSU2SpinSite(Bulk);

   // construct a lattice of L copies of Site
   Lattice MyLattice = join(EdgeSite, repeat(BulkSite, L-2), EdgeSite);
   MyLattice.fix_coordinates();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   OperatorAtSite<OperatorList const, int> S(OpList, "S");
   OperatorAtSite<OperatorList, int> Kappa(OpList, "Kappa");
   MPOperator& Hamiltonian = OpList["H"];
 
   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   QuantumNumber Vec(MyLattice.GetSymmetryList(), "1"); // S=1 quantum number

   // interaction matrix elements.
   // minus sign here because the S operator is skew-hermitian;
   // what we actually intend is S * herm(S) = - S * S
   Hamiltonian = -J1 * CreateRepeatedOperator(MyLattice, "S", "S")
      - J2 * CreateRepeatedOperator(MyLattice, "S", "I", "S");

   std::cout << "Constructing Kappa operators...\n";
   for (int i = 1; i < L; ++i)
   {
      //      Hamiltonian += -J1 * sqrt(3.0) * prod(S(i), S(i%L+1), Ident);
      //if (i < L-1)
      //Hamiltonian += -J2 * sqrt(3.0) * prod(S(i), S((i+1)%L+1), Ident);

      // Chiral operators
      // In principle, there should be a factor -i here,
      // A x B = sqrt(2) * (-i) * (A*B)^{[1]}
      // We omit the imaginary factor as it will cancel out when we take the scalar product
      // anyway.
      Kappa(i) = sqrt(2.0) * prod(S(i), S(i%L+1), Vec);
   }

#if 0
   // This part is quadratic in the length.  Only enable it if we really need
   // these operators.
   std::cout << "Constructing bipartite operators...\n";

   MPOperator& TotalS = OpList["S"];
   MPOperator& TotalS_a = OpList["S_a"];      // spin on A and B bipartite lattices
   MPOperator& TotalS_b = OpList["S_b"];

   // construct the total spin operator, and the bipartite operators
   for (int i = 1; i <= L; ++i)
   {
      TotalS += S(i);
      if (i % 2 == 1)
         TotalS_a += S(i);
      else
         TotalS_b += S(i);
   }

   std::cout << "Constructing scalar products...\n";
   // and S^2
   OpList["S2"] = -sqrt(3.0) * prod(TotalS, TotalS, Ident);
   OpList["S_a^2"] = -sqrt(3.0) * prod(TotalS_a, TotalS_a, Ident);
   OpList["S_b^2"] = -sqrt(3.0) * prod(TotalS_b, TotalS_b, Ident);
   OpList["S_aS_b"] = -sqrt(3.0) * prod(TotalS_a, TotalS_b, Ident);
#endif

   std::cout << "done." << std::endl;

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::ShutdownPersistent(OList);
}
