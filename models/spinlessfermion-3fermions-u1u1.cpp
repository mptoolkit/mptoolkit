// 3 fermion lattice system
// based on -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/spinlessfermion-u1.h"
#include "math.h"

int main(int argc, char** argv)
{
   if (argc != 13)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: spinlessfermion-u1 <L> <t1> <t2> <t3> <U12> <U13> <U23> <m1> <m2> <m3> <k> <outfile>\n"
                << "L = number of lattice sites\n"
                << "t1 = hopping for level 1\n"
                << "t2 = hopping for level 2\n"
                << "t3 = hopping for level 3\n"
                << "U12 = interaction between levels 1,2\n"
                << "U13 = interaction between levels 1,3\n"
                << "U23 = interaction between levels 2,3\n"
                << "m1 = chemical potential for level 1\n"
                << "m2 = chemical potential for level 2\n"
                << "m3 = chemical potential for level 3\n"
                << "k = momentum value for pertubation\n"
                << "outfile = file name for output lattice.\n";
      return 1;
   }

   int L = boost::lexical_cast<int>(argv[1]);
   double t1 = boost::lexical_cast<double>(argv[2]);
   double t2 = boost::lexical_cast<double>(argv[3]);
   double t3 = boost::lexical_cast<double>(argv[4]);
   double U12 = boost::lexical_cast<double>(argv[5]);
   double U13 = boost::lexical_cast<double>(argv[6]);
   double U23 = boost::lexical_cast<double>(argv[7]);
   double m1 = boost::lexical_cast<double>(argv[8]);
   double m2 = boost::lexical_cast<double>(argv[9]);
   double m3 = boost::lexical_cast<double>(argv[10]);
   double k = boost::lexical_cast<double>(argv[11]);

   TRACE(L)(t1)(t2)(t3)(U12)(U13)(U23)(m1)(m2)(m3)(k);

   CHECK(L >= 2);
   // Construct the site block
   SiteBlock SiteA = CreateU1SpinlessFermion("NA");
   SiteBlock SiteBC = CreateU1SpinlessFermion("NBC");
   
   // Create unit cell
   Lattice UnitCell(SymmetryList("NA:U(1),NBC:U(1)"), SiteA, SiteBC, SiteBC);
   UnitCell.fix_coordinates("1", "2", "3");

   // construct a lattice of L copies of Site
   Lattice MyLattice = repeat(UnitCell, L);
   MyLattice.fix_coordinates_prepend();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);
   OperatorAtSite<OperatorList const, int, int> CH(OpList, "CH");
   OperatorAtSite<OperatorList const, int, int> C(OpList, "C");
   OperatorAtSite<OperatorList const, int, int> N(OpList, "N");
   OperatorAtSite<OperatorList, int> SZ(OpList, "SZ");
   MPOperator& Hamiltonian = OpList["H"];
   MPOperator& Perturbation = OpList["P"];
   MPOperator& PerturbationH = OpList["PH"];
   MPOperator& Perturbation1 = OpList["P1"];
   MPOperator& Perturbation1H = OpList["P1H"];

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   // hopping matrix elements
   std::cout << "generating hopping matrix elements\n";
   for (int i = 1; i < L; ++i) {
	  MPOperator Hopping = -t1 * (prod(CH(i,1), C(i%L+1,1), Ident) + prod(CH(i%L+1,1), C(i,1), Ident));
	  Hamiltonian += Hopping;
	  Hamiltonian += -t2 * (prod(CH(i,2), C(i%L+1,2), Ident) + prod(CH(i%L+1,2), C(i,2), Ident));
	  Hamiltonian += -t3 * (prod(CH(i,3), C(i%L+1,3), Ident) + prod(CH(i%L+1,3), C(i,3), Ident));
	  std::cout << "Working.... " << i << "\n";
   }
   
   // interaction
   std::cout << "generating interaction, chemical potential and pertubation matrix elements\n";
   for (int i = 1; i <= L; ++i) {
	   Hamiltonian += U12 * prod(N(i,1),N(i,2),Ident);
	   Hamiltonian += U13 * prod(N(i,1),N(i,3),Ident);
	   Hamiltonian += U23 * prod(N(i,2),N(i,3),Ident);
	   
	   Hamiltonian += m1 * N(i,1);
	   Hamiltonian += m2 * N(i,2);
	   Hamiltonian += m3 * N(i,3);
	   
	   std::complex<double> re(1,0);
	   std::complex<double> im(0,1);
	   MPOperator Perturb = (re*cos(k*i)-im*sin(k*i)) * (prod(CH(i,2),C(i,3),Ident));
	   Perturbation += Perturb;
           MPOperator Perturb1 = (re*cos(k*i)+im*sin(k*i)) * (prod(CH(i,2),C(i,3),Ident));
	   Perturbation1 += Perturb1;

	   std::cout << "Working.... " << i << "\n";
   }

   PerturbationH = adjoint(Perturbation);
   Perturbation1H = adjoint(Perturbation1);

   // build useful operators
   for (int i=1; i<=L; ++i) {
	SZ(i) = N(i,2) - N(i,1);
   }
   
   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(argv[12], 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
