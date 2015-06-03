// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "models/hubbard-so4.h"
#include "mp/copyright.h"
#include "pstream/pfilestream.h"

std::string Coord(int i, int j)
{
   std::ostringstream S;
   S << i << ',' << j;
   return S.str();
}

typedef std::complex<double> complex;

int main(int argc, char** argv)
{
   if (argc != 7)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: hubbardladder-so4 <L> <N-legs> <t-long> <t-trans> <U> <outfile>\n"
                << "L = length of ladder (number of rungs)\n"
                << "N-legs = number of legs\n"
                << "t-long = hopping between rungs\n"
                << "t-trans = hopping between ladders\n"
                << "U = onside repulsion\n"
                << "outfile = file name for output lattice.\n";
      return 1;
   }

   int L = boost::lexical_cast<int>(argv[1]);
   int N = boost::lexical_cast<int>(argv[2]);
   double t_long = boost::lexical_cast<double>(argv[3]);
   double t_trans = boost::lexical_cast<double>(argv[4]);
   double U = boost::lexical_cast<double>(argv[5]);

   // initialize the output file	
   pheap::Initialize(argv[6], 1, 65536, 655360);

   Lattice UnitCellAB(CreateSO4HubbardSiteA(), CreateSO4HubbardSiteB());
   Lattice UnitCellBA(CreateSO4HubbardSiteB(), CreateSO4HubbardSiteA());
   Lattice RungAB = repeat(UnitCellAB, N/2);
   Lattice RungBA = repeat(UnitCellBA, N/2);
   if (N%2 == 1)  // if the width is odd, join on an extra site
   {
      RungAB = join(RungAB, Lattice(CreateSO4HubbardSiteA()));
      RungBA = join(RungBA, Lattice(CreateSO4HubbardSiteB()));
   }
   RungAB.fix_coordinates();
   RungBA.fix_coordinates();

   Lattice SuperCell = join(RungAB, RungBA);
   Lattice MyLattice = repeat(SuperCell, L/2);
   if (L%2 == 1)
      MyLattice = join(MyLattice, RungAB);

   MyLattice.fix_coordinates_prepend();  // make the main coordinate the first one

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);
	
   OperatorAtSite<OperatorList const, int, int> CH(OpList, "CH");
   OperatorAtSite<OperatorList const, int, int> C(OpList, "C");
   OperatorAtSite<OperatorList const, int, int> P(OpList, "P");

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   MPOperator& Hamiltonian = OpList["H"];

   // longitudinal matrix elements
   for (int i = 1; i <= L; ++i)
   {
      for (int j = 1; j <= N; ++j)
      {
         Hamiltonian += (-2.) * t_long * prod(C(i%L+1,j), CH(i,j), Ident);
      }
   }

   // transverse matrix elements
   for (int i = 1; i <= L; ++i)
   {
      for (int j = 1; j < N; ++j)
      {
	 Hamiltonian += (-2.) * t_trans * prod(C(i,j), CH(i,j+1), Ident);
      }
   }
   
   //Coulomb repulsion
   for (int i = 1; i <= L; ++i)
   {
      for (int j = 1; j <= N; ++j)
      {
         Hamiltonian += (U / 4.) * P(i,j);
      }
   }

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);
   pheap::ShutdownPersistent(OList);
}
