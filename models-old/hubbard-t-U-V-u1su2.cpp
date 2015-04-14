
#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/hubbard-u1su2.h"

typedef std::complex<double> complex;

int main(int argc, char** argv)
{
   if (argc != 6)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: hubbard-t-U-V-u1su2 <L> <t> <U> <V> <outfile>\n"
                << "L = number of lattice sites\n"
                << "t = hopping integral\n"
                << "U = coupling constant\n"
                << "V = Coulomb interaction\n"
                << "outfile = file name for output lattice.\n";
      return 1;
   }

   int L = boost::lexical_cast<int>(argv[1]);
   std::complex<double> t = boost::lexical_cast<complex>(argv[2]);
   double U = boost::lexical_cast<double>(argv[3]);
   double V = boost::lexical_cast<double>(argv[4]);

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
   MPOperator& Hamiltonian = OpList["H"];

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   // hopping matrix elements

   complex Sqrt2t = -std::sqrt(2.0) * t;

   for (int i = 1; i < L; ++i)
   {
      MPOperator Hopping = Sqrt2t * prod(CH(i), C(i%L+1), Ident) + conj(Sqrt2t) * prod(C(i), CH(i%L+1), Ident);
      Hamiltonian += Hopping;
      std::cout << "Working.... " << i << "\n";
   }
   // onside coulomb repulsion
   for (int i = 1; i <= L; ++i)
   {
      Hamiltonian += (U/4.0) * P(i);

      std::cout << "Working.... " << i << "\n";
   }
   
   // next-nearest neighbour coulomb repulsion
   for (int i = 1; i < L; ++i)
   {
      Hamiltonian += V * prod(N(i), N(i%L+1), Ident);

      std::cout << "Working.... " << i << "\n";
   }
   

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(argv[5], 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
