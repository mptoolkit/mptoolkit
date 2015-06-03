// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/hubbard-u1u1.h"
#include "models/site-tensor.h"

using std::string;

int main(int argc, char** argv)
{
   if (argc != 5)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: hubbard-t-u1u1 <L> <t> <U> <outfile>\n"
                << "L = number of lattice sites\n"
                << "t = hopping integral\n"
                << "U = coupling constant\n"
                << "outfile = file name for output lattice.\n";
      return 1;
   }

   int L = boost::lexical_cast<int>(argv[1]);
   double t = boost::lexical_cast<double>(argv[2]);
   double U = boost::lexical_cast<double>(argv[3]);

   TRACE(L)(t)(U);
   
   // Construct the site block
   SiteBlock BasicSite = CreateU1HubbardSite();
   SiteBlock Site = SiteTensorProduct(BasicSite, BasicSite);
   // construct a lattice of L copies of Site
   Lattice MyLattice = repeat(Site, L);
   MyLattice.fix_coordinates();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   OperatorAtSite<OperatorList const, int> CHup(OpList, "CHup*I");
   OperatorAtSite<OperatorList const, int> Cup(OpList, "Cup*I");
   OperatorAtSite<OperatorList const, int> CHdown(OpList, "CHdown*I");
   OperatorAtSite<OperatorList const, int> Cdown(OpList, "Cdown*I");
   OperatorAtSite<OperatorList const, int> Sp(OpList, "Sp*I");
   OperatorAtSite<OperatorList const, int> Sm(OpList, "Sm*I");
   OperatorAtSite<OperatorList const, int> Sz(OpList, "Sz*I");
   OperatorAtSite<OperatorList const, int> Qp(OpList, "Qp*I");
   OperatorAtSite<OperatorList const, int> Qm(OpList, "Qm*I");
   OperatorAtSite<OperatorList const, int> Qz(OpList, "Qz*I");
   OperatorAtSite<OperatorList const, int> P(OpList, "P*I");
   OperatorAtSite<OperatorList const, int> AuxCHup(OpList, "P*CHup");
   OperatorAtSite<OperatorList const, int> AuxCup(OpList, "P*Cup");
   OperatorAtSite<OperatorList const, int> AuxCHdown(OpList, "P*CHdown");
   OperatorAtSite<OperatorList const, int> AuxCdown(OpList, "P*Cdown");
   OperatorAtSite<OperatorList const, int> AuxSp(OpList, "I*Sp");
   OperatorAtSite<OperatorList const, int> AuxSm(OpList, "I*Sm");
   OperatorAtSite<OperatorList const, int> AuxSz(OpList, "I*Sz");
   OperatorAtSite<OperatorList const, int> AuxQp(OpList, "I*Qp");
   OperatorAtSite<OperatorList const, int> AuxQm(OpList, "I*Qm");
   OperatorAtSite<OperatorList const, int> AuxQz(OpList, "I*Qz");

   OperatorAtSite<OperatorList, int> Bond(OpList, "Bond");
   MPOperator& Hamiltonian = OpList["H"];
   MPOperator& IG = OpList["IG"];
   MPOperator& IGH = OpList["IGH"];
   MPOperator IGX;

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   // hopping matrix elements
   for (int i = 1; i < L; ++i)
   {
      MPOperator Hopping = -t * (
				     prod(CHup(i), Cup(i%L+1), Ident)
				   - prod(Cup(i), CHup(i%L+1), Ident)
				   + prod(CHdown(i), Cdown(i%L+1), Ident)
				   - prod(Cdown(i), CHdown(i%L+1), Ident)
				 );
      IGX -= t * (
		  prod(AuxCHup(i), AuxCup(i%L+1), Ident)
		  - prod(AuxCup(i), AuxCHup(i%L+1), Ident)
		  + prod(AuxCHdown(i), AuxCdown(i%L+1), Ident)
		  - prod(AuxCdown(i), AuxCHdown(i%L+1), Ident)
		  );
      Bond(i) = Hopping;
      Hamiltonian += Bond(i);
      std::cout << "Working.... " << i << "\n";
   }
   // coulomb repulsion
   for (int i = 1; i <= L; ++i)
   {
      Hamiltonian += (U/4.0) * P(i);
      if (i < L) 
         Bond(i) += (U/4.0) * P(i);
      else
         Bond(i-1) += (U/4.0) * P(i);
      MPOperator X = -t * (
			   prod(AuxCHup(i), Cup(i), Ident)
			   - prod(AuxCup(i), CHup(i), Ident)
			   + prod(AuxCHdown(i), Cdown(i), Ident)
			   - prod(AuxCdown(i), CHdown(i), Ident)
			   );
      IGX += X;
      IG = IG + X + 0.5 * (prod(AuxSp(i), Sm(i), Ident) + prod(AuxSm(i), Sp(i), Ident))
	 + prod(AuxSz(i), Sz(i), Ident)
	 + 0.5 * (prod(AuxQp(i), Qm(i), Ident) + prod(AuxQm(i), Qp(i), Ident))
	 + prod(AuxQz(i), Qz(i), Ident);
      std::cout << "Working.... " << i << "\n";
   }

   IGH = IG + 0.1 * Hamiltonian;

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(argv[4], 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
