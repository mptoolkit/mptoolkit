// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/spin-su2.h"

std::string Coord(int i, int j)
{
   std::ostringstream S;
   S << i << ',' << j;
   return S.str();
}

int main(int argc, char** argv)
{
   if (argc != 7)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: su2-heisladder <Spin> <L> <N-legs> <J-long> <J-trans> <outfile>\n"
                << "Spin = spin of the bulk (half-integer)\n"
                << "L = length of ladder (number of rungs)\n"
                << "N-legs = number of legs\n"
                << "J-long = coupling between rungs\n"
                << "J-trans = coupling between ladders\n"
                << "outfile = file name for output lattice.\n";
      return 1;
   }

   half_int S = boost::lexical_cast<half_int>(argv[1]);
   int L = boost::lexical_cast<int>(argv[2]);
   int N = boost::lexical_cast<int>(argv[3]);
   double Jlong = boost::lexical_cast<double>(argv[4]);
   double Jtrans = boost::lexical_cast<double>(argv[5]);

   // Construct the site block
   SiteBlock Site = CreateSU2SpinSite(S);

   // construct the lattice
   Lattice Leg = repeat(Site, N);
   Leg.fix_coordinates();

   Lattice MyLattice = repeat(Leg, L);
   MyLattice.fix_coordinates_prepend();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   MPOperator Hamiltonian;

   // longitudinal matrix elements
   for (int i = 1; i < L; ++i)
   {
      for (int j = 1; j <= N; ++j)
      {
         Hamiltonian = Hamiltonian + (-1) * std::sqrt(3.0) * Jlong *
            prod(OpList["S("+Coord(i,j)+")"], OpList["S("+Coord(i+1,j)+")"], Ident);
      }
   }

   // transverse matrix elements
   for (int i = 1; i <= L; ++i)
   {
      for (int j = 1; j < N; ++j)
      {
         Hamiltonian = Hamiltonian + (-1) * std::sqrt(3.0) * Jtrans *
            prod(OpList["S("+Coord(i,j)+")"], OpList["S("+Coord(i,j+1)+")"], Ident);
      }
   }

   // save the Hamiltonian into the operator list
   OpList["H"] = Hamiltonian;

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(argv[6], 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
