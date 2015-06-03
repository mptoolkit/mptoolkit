// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/hubbard-so4.h"

typedef std::complex<double> complex;

int main(int argc, char** argv)
{
   if (argc != 5)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: so4hubbard <L> <t> <U> <outfile>\n"
                << "L = number of lattice sites\n"
                << "t = hopping integral\n"
                << "U = coupling constant\n"
                << "outfile = file name for output lattice.\n";
      return 1;
   }

   int L = boost::lexical_cast<int>(argv[1]);
   complex t = boost::lexical_cast<complex>(argv[2]);
   double U = boost::lexical_cast<double>(argv[3]);

   TRACE(L)(t)(U);

   // Construct the site blocks
   SiteBlock SiteA = CreateSO4HubbardSiteA();
   SiteBlock SiteB = CreateSO4HubbardSiteB();

   Lattice UnitCellAB(SiteA, SiteB);
   UnitCellAB.fix_coordinates("", "aux");

   Lattice UnitCellBA(SiteB, SiteA);
   UnitCellBA.fix_coordinates("", "aux");

   Lattice SuperCell(UnitCellAB, UnitCellBA);

   // construct a lattice of L copies of Site
   Lattice MyLattice = repeat(SuperCell, L/2);
   if (L%2 == 1)
      MyLattice = join(MyLattice, UnitCellAB);
   MyLattice.fix_coordinates();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   OperatorAtSite<OperatorList const, int> C(OpList, "C");
   OperatorAtSite<OperatorList const, int> CH(OpList, "CH");
   OperatorAtSite<OperatorList const, int> P(OpList, "P");
   OperatorAtSite<OperatorList const, int> S(OpList, "S");
   OperatorAtSite<OperatorList const, int> Q(OpList, "Q");
   OperatorAtSite<OperatorList const, std::string, int> AuxC(OpList, "C");
   OperatorAtSite<OperatorList const, std::string, int> AuxCH(OpList, "CH");
   OperatorAtSite<OperatorList const, std::string, int> AuxS(OpList, "S");
   OperatorAtSite<OperatorList const, std::string, int> AuxQ(OpList, "Q");
   MPOperator& Hamiltonian = OpList["H"];
   MPOperator& IG = OpList["IG"];
   MPOperator& IGH = OpList["IGH"];
   MPOperator IGX;

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   // hopping matrix elements

   complex HoppingValue = -2.0 * t;

   for (int i = 1; i < L; ++i)
   {
      MPOperator Hopping = HoppingValue * prod(C(i), CH(i%L+1), Ident);
      Hamiltonian += Hopping;
      IGX += HoppingValue * prod(AuxC("aux",i), AuxCH("aux",i%L+1), Ident);
      std::cout << "Working.... " << i << "\n";
   }
   // coulomb repulsion
   for (int i = 1; i <= L; ++i)
   {
      Hamiltonian = Hamiltonian + (U/4.0) * P(i);
      MPOperator X = prod(C(i), AuxCH("aux", i), Ident);
      IGX += X;
      IG = IG + X - 1.0 * std::sqrt(3.0) * prod(AuxS("aux",i), S(i), Ident)
	 - 1.0 * std::sqrt(3.0) * prod(AuxQ("aux", i), Q(i), Ident);
      std::cout << "Working.... " << i << "\n";
   }

   IGH = IG + 0.1 * (Hamiltonian + IGX);

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(argv[4], 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
