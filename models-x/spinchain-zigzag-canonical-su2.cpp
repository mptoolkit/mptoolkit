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
      std::cerr << "usage: mp-spinchain-zigzag-canonical-su2 <L> <EdgeSpin> <BulkSpin> <J1> <J2> <outfile>\n"
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
   
   SiteBlock AuxEdgeSite = CreateSU2SpinSite(Edge,"Saux");
   SiteBlock AuxBulkSite = CreateSU2SpinSite(Bulk,"Saux");


   // construct a lattice of L copies of Site
   Lattice EdgeUnitCell(SymmetryList("S:SU(2),Saux:SU(2)"), EdgeSite, AuxEdgeSite);
   EdgeUnitCell.fix_coordinates("", "aux");
   Lattice BulkUnitCell(SymmetryList("S:SU(2),Saux:SU(2)"), BulkSite, AuxBulkSite);
   BulkUnitCell.fix_coordinates("", "aux");
  
   Lattice MyLattice = join(EdgeUnitCell, repeat(BulkUnitCell, L-2), EdgeUnitCell);
   MyLattice.fix_coordinates_prepend();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   OperatorAtSite<OperatorList const> S(OpList, "S");
   MPOperator& Hamiltonian = OpList["H"];
   MPOperator& AuxHamiltonian = OpList["AuxH"];
 
   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number

   // interaction matrix elements
   for (int i = 1; i < L; ++i)
      Hamiltonian += - J1 * sqrt(3.) *  prod(S(i), S(i+1), Ident);
   
   for (int i = 1; i < L-1; ++i)
      Hamiltonian += - J2 * sqrt(3.) *  prod(S(i), S(i+2), Ident);
   
   for (int i = 1; i < L; ++i)
      AuxHamiltonian += - J1 * sqrt(3.) *  prod(S(i,"aux"), S(i+1,"aux"), Ident);
   
   for (int i = 1; i < L-1; ++i)
      AuxHamiltonian += - J2 * sqrt(3.) *  prod(S(i,"aux"), S(i+2,"aux"), Ident);

      
         

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::ShutdownPersistent(OList);
}
