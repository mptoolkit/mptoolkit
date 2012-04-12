// -*- C++ -*- $Id: bosehubbard-spinless-u1.cpp 920 2008-06-24 18:04:10Z ianmcc $

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/bosehubbard-spinless.h"

int main(int argc, char** argv)
{
   if (argc != 10)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: bosehubbard-spinless <MaxN> <U> <Nwell> <K> <V_L> <a_L> <L> <Ring> <outfile>\n"
                << "MaxN - maximum number of bosons per site\n"
                << "U  - on-site coulomb strength\n"
                << "Nwell - number of wells in the trap, possible is 0,1,2\n"
                << "K - trap strength\n"
                << "V_L - optical lattice strength\n"
                << "a_L - number of computational sites per optical lattice spacing (integer)\n"
                << "L  - final number of lattice sites\n"
                << "Ring  - 0 for fixed boundary conditions, 1 for periodic boundaries\n"
                << "outfile = file name for output lattice.\n";
      return 1;
   }

   int MaxN = boost::lexical_cast<int>(argv[1]);
   double U = boost::lexical_cast<double>(argv[2]);
   int Nwell = boost::lexical_cast<int>(argv[3]);
   double K = boost::lexical_cast<double>(argv[4]);
   double V_L = boost::lexical_cast<double>(argv[5]);
   int a_L = boost::lexical_cast<double>(argv[6]);
   int L = boost::lexical_cast<int>(argv[7]);
   double Ring = boost::lexical_cast<double>(argv[8]);
   std::string FileName = argv[9];

   TRACE(MaxN)(U)(K)(L);

   // Construct the site block
   SiteBlock Site = CreateBoseHubbardSpinlessSite(MaxN);

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
   MPOperator& Hamiltonian_ring = OpList["H_ring"];
   MPOperator& HV = OpList["H_V"];
   MPOperator& HVV = OpList["H_VV"];
   MPOperator& H_mu = OpList["H_mu"];
   MPOperator& H_lattice = OpList["H_lattice"];
   MPOperator& HTT = OpList["H_TT"];
   MPOperator& H_tilt = OpList["H_tilt"];

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number


   // hopping matrix elements
   for (int i = 1; i < L; ++i)
   {
      Hamiltonian += -1.0 * (prod(BH(i), B(i%L+1), Ident)
                           + prod(B(i), BH(i%L+1), Ident));

      HV += N(i) * N(i+1);
      if (i < L-1)
      {
         HVV += N(i) * N(i+2);
         HTT += prod(BH(i), B(i+2), Ident) + prod(B(i), BH(i+2), Ident); // next nearest neighbour hopping
      }

      H_mu += N(i);

      std::cout << "Working.... " << i << "\n";
   }

   if (Ring != 0) Hamiltonian_ring += -Ring * (prod(BH(L), B(1), Ident)
                                             + prod(B(L), BH(1), Ident));

   // additional optical lattice
   if ((a_L != 0) && (V_L != 0.0))
   {
      for (int i = 1; i <= L; ++i)
      {
         H_lattice += 0.5 * V_L * cos(2*M_PI*(i-1)/a_L) * N(i);
      }
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

      H_tilt += (i - L/2.0 - 0.5) * N(i); // simple linear potential

      std::cout << "Working.... " << i << "\n";
   }

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(FileName, 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
