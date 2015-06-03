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
   if (argc != 3)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: hubbard-tj-u1su2 <L> <outfile>\n"
                << "L = number of lattice sites\n"
                << "outfile = file name for output lattice.\n";
      return 1;
   }

   int L = boost::lexical_cast<int>(argv[1]);
   std::string OutFile = argv[2];

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
   OperatorAtSite<OperatorList const, int> S(OpList, "S");
   OperatorAtSite<OperatorList const, int> LocalPg(OpList, "Pg");
   OperatorAtSite<OperatorList, int> Bond(OpList, "Bond");
   OperatorAtSite<OperatorList, int> Cg(OpList, "Cg");
   OperatorAtSite<OperatorList, int> CHg(OpList, "CHg");
   OperatorAtSite<OperatorList, int> Ng(OpList, "Ng");
   MPOperator& Pg = OpList["Pg"];                // Gutzwiller projector
   MPOperator& Hop = OpList["H"];
   MPOperator& Coulomb = OpList["H_U"];
   MPOperator& CoulombV = OpList["H_V"];
   MPOperator& CoulombVg = OpList["H_Vg"];
   MPOperator& Hg = OpList["H_g"];
   MPOperator& HJ = OpList["H_J"];
   OperatorAtSite<OperatorList, int, int> PairHopg(OpList, "PairHopg");
   OperatorAtSite<OperatorList, int, int> SpinHop(OpList, "SpinHop");
   OperatorAtSite<OperatorList, int, int> VgHop(OpList, "VgHop");

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   // hopping matrix elements

   Pg = OpList["I"];
   // Gutswiller projectors
   for (int i = 1; i <= L; ++i)
   {
      Pg = Pg * LocalPg(i);  // Gutzwiller projector
      Cg(i) = LocalPg(i) * C(i) * LocalPg(i);
      CHg(i) = LocalPg(i) * CH(i) * LocalPg(i);
      Ng(i) = LocalPg(i) * N(i); // these operators commute
   }

   for (int i = 1; i < L; ++i)
   {
      Hop -= std::sqrt(2.0) * (prod(CH(i), C(i%L+1), Ident) +  prod(C(i), CH(i%L+1), Ident));
      Hg -= std::sqrt(2.0) * (prod(CHg(i), Cg(i%L+1), Ident) +  prod(Cg(i), CHg(i%L+1), Ident));
      HJ -= std::sqrt(3.0) * (prod(S(i), S(i%L+1), Ident));  // -sqrt(3) is the coupling coefficient
      CoulombV += prod(N(i), N(i%L+1), Ident);
      CoulombVg += prod(Ng(i), Ng(i%L+1), Ident);
      std::cout << "Working.... " << i << "\n";
   }
   // Coulomb
   for (int i = 1; i <= L; ++i)
   {
      Coulomb += 0.25 * P(i);
      std::cout << "Working.... " << i << "\n";
   }

   // pair hoppings
   for (int i = 1; i < L; ++i)
   {                                         
      for (int j = i+1; j <= L; ++j)
      {
         PairHopg(i,j) = -sqrt(2.0)*(prod(CHg(i), Cg(j), Ident) +  prod(Cg(i), CHg(j), Ident));
         SpinHop(i,j) = -sqrt(3.0)*prod(S(i), S(j), Ident);
         VgHop(i,j) = prod(Ng(i), Ng(j), Ident);
      }
   }

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(OutFile, 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
