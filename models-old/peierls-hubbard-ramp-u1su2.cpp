// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/hubbard-u1su2.h"

typedef std::complex<double> complex;

int main(int argc, char** argv)
{
   if (argc != 7)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: su2hubbard <L> <t> <U> <outfile>\n"
                << "L = number of lattice sites\n"
                << "t = hopping integral\n"
                << "U = coupling constant\n"
                << "Delta = bond alternation\n"
                << "V = ramp potential\n"
                << "outfile = file name for output lattice.\n";
      return 1;
   }

   int L = boost::lexical_cast<int>(argv[1]);
   std::complex<double> t = boost::lexical_cast<complex>(argv[2]);
   double U = boost::lexical_cast<double>(argv[3]);
   double Delta = boost::lexical_cast<double>(argv[4]);
   double V = boost::lexical_cast<double>(argv[5]);

   TRACE(L)(t)(U)(Delta)(V);

   // Construct the site block
   SiteBlock Site = CreateSU2HubbardSite();

   // construct a lattice of L copies of Site
   Lattice MyLattice = repeat(Site, L);
   MyLattice.fix_coordinates();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   OperatorAtSite<OperatorList const, int> CH(OpList, "CH");
   OperatorAtSite<OperatorList const, int> C(OpList, "C");
   OperatorAtSite<OperatorList const, int> P(OpList, "P");
   OperatorAtSite<OperatorList const, int> N(OpList, "N");
   MPOperator& Hamiltonian = OpList["H"];
   MPOperator& Hop = OpList["Hopping"];
   MPOperator& Coulomb = OpList["Coulomb"];

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   // hopping matrix elements

   for (int i = 1; i < L; ++i)
   {
      MPOperator Hopping 
         = -(1.0 + Delta * minus1pow(i)) * 
         (std::sqrt(2.0) * prod(CH(i), C(i%L+1), Ident) + std::sqrt(2.0) * prod(C(i), CH(i%L+1), Ident));
      Hop += Hopping;
      Hopping *= t;
      Hamiltonian += Hopping;
      std::cout << "Working.... " << i << "\n";
   }
   // coulomb repulsion
   for (int i = 1; i <= L; ++i)
   {
      Coulomb += 0.25 * P(i);
      Hamiltonian = Hamiltonian + (U/4.0) * P(i);

      Hamiltonian += (-V/2 + (i-1) * V / (L-1)) * N(i);

      std::cout << "Working.... " << i << "\n";
   }

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(argv[6], 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
