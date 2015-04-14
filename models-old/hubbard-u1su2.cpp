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
   if (argc != 5)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: hubbard-u1su2 <L> <t> <U> <outfile>\n"
                << "L = number of lattice sites\n"
                << "t = hopping integral\n"
                << "U = coupling constant\n"
                << "outfile = file name for output lattice.\n";
      return 1;
   }

   int L = boost::lexical_cast<int>(argv[1]);
   std::complex<double> t = boost::lexical_cast<complex>(argv[2]);
   double U = boost::lexical_cast<double>(argv[3]);

   TRACE(L)(t)(U);

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
   OperatorAtSite<OperatorList const, int> LocalPg(OpList, "Pg");
   OperatorAtSite<OperatorList, int> Bond(OpList, "Bond");
   MPOperator& Hamiltonian = OpList["H"];
   MPOperator& Hop = OpList["Hopping"];
   MPOperator& Coulomb = OpList["Coulomb"];
   MPOperator& Pg = OpList["Pg"];                // Gutzwiller projector

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   // hopping matrix elements

   complex Sqrt2t = -std::sqrt(2.0) * t;

   for (int i = 1; i < L; ++i)
   {
      MPOperator Hopping 
         = Sqrt2t * prod(CH(i), C(i%L+1), Ident) + conj(Sqrt2t) * prod(C(i), CH(i%L+1), Ident);
      Hop += Hopping;
      Hamiltonian += Hopping;
      Bond(i) = Hopping;
      std::cout << "Working.... " << i << "\n";
   }
   // coulomb repulsion
   for (int i = 1; i <= L; ++i)
   {
      Pg = Pg * LocalPg(i);  // Gutzwiller projector

      Coulomb += 0.25 * P(i);
      Hamiltonian = Hamiltonian + (U/4.0) * P(i);

      // distribute the interaction among the bond terms.
      // There are choices in how to do this.      
      if (i < L) 
         Bond(i) += (U/4.0) * P(i);
      else
         Bond(i-1) += (U/4.0) * P(i);

      std::cout << "Working.... " << i << "\n";
   }

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(argv[4], 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
