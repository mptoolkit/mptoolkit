// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/hubbard-u1u1.h"
#include "models/spin-u1.h"

typedef std::complex<double> complex;

int main(int argc, char** argv)
{
   if (argc != 6)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: klm-u1u1 <L> <t> <J_K> <Jz_H> <outfile>\n"
                << "L = number of lattice sites\n"
                << "t = hopping integral\n"
                << "J_K = Kondo coupling\n"
                << "Jz_H = Ising coupling of the localized spins\n"
                << "outfile = file name for output lattice.\n";
      return 1;
   }

   int L = boost::lexical_cast<int>(argv[1]);
   complex t = boost::lexical_cast<complex>(argv[2]);
   double Jk = boost::lexical_cast<double>(argv[3]);
   double Jzh = boost::lexical_cast<double>(argv[4]);

   // initialize the output file
   pheap::Initialize(argv[5], 1, 65536, 655360);

   TRACE(L)(t)(Jk)(Jzh);

   // Construct the site blocks
   SiteBlock fSite = CreateU1SpinSite(0.5);
   SiteBlock cSite = CreateU1HubbardSite();

   // construct the lattice
   Lattice UnitCell(SymmetryList("N:U(1),Sz:U(1)"), fSite, cSite);
   UnitCell.fix_coordinates("f", "c");

   Lattice MyLattice = repeat(UnitCell, L);
   MyLattice.fix_coordinates();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   OperatorAtSite<OperatorList const, char, int> UnrolledCup(OpList, "Cup");
   OperatorAtSite<OperatorList const, char, int> UnrolledCHup(OpList, "CHup");
   OperatorAtSite<OperatorList const, char, int> UnrolledCdown(OpList, "Cdown");
   OperatorAtSite<OperatorList const, char, int> UnrolledCHdown(OpList, "CHdown");

   OperatorAtSite<OperatorList const, char, int> UnrolledSp(OpList, "Sp");
   OperatorAtSite<OperatorList const, char, int> UnrolledSm(OpList, "Sm");
   OperatorAtSite<OperatorList const, char, int> UnrolledSz(OpList, "Sz");
   MPOperator& Hamiltonian = OpList["H"];
   MPOperator& H_t = OpList["H_t"];
   MPOperator& H_Jz = OpList["H_Jz"];
   MPOperator& H_Jk = OpList["H_Jk"];
   MPOperator& H_Jzh = OpList["H_Jzh"];

   // Make some combined operators, for convenience
   OperatorAtSite<OperatorList, int> Sp(OpList, "Sp");
   OperatorAtSite<OperatorList, int> Sm(OpList, "Sm");
   OperatorAtSite<OperatorList, int> Sz(OpList, "Sz");

   OperatorAtSite<OperatorList, int> Cup(OpList, "Cup");
   OperatorAtSite<OperatorList, int> CHup(OpList, "CHup");
   OperatorAtSite<OperatorList, int> Cdown(OpList, "Cdown");
   OperatorAtSite<OperatorList, int> CHdown(OpList, "CHdown");

   OperatorAtSite<OperatorList, int> LocalSp(OpList, "fSp");
   OperatorAtSite<OperatorList, int> LocalSm(OpList, "fSm");
   OperatorAtSite<OperatorList, int> LocalSz(OpList, "fSz");

   OperatorAtSite<OperatorList, int> CondSp(OpList, "cSp");
   OperatorAtSite<OperatorList, int> CondSm(OpList, "cSm");
   OperatorAtSite<OperatorList, int> CondSz(OpList, "cSz");

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number

   // Map the original operators into the per-site operators
   for (int i = 1; i <= L; ++i)
   {
      Sp(i) = UnrolledSp('f', i) + UnrolledSp('c', i);
      Sm(i) = UnrolledSm('f', i) + UnrolledSm('c', i);
      Sz(i) = UnrolledSz('f', i) + UnrolledSz('c', i);

      Cup(i) = UnrolledCup('c', i);
      CHup(i) = UnrolledCHup('c', i);
      Cdown(i) = UnrolledCdown('c', i);
      CHdown(i) = UnrolledCHdown('c', i);

      LocalSp(i) = UnrolledSp('f', i);
      LocalSm(i) = UnrolledSm('f', i);
      LocalSz(i) = UnrolledSz('f', i);

      CondSp(i) = UnrolledSp('c', i);
      CondSm(i) = UnrolledSm('c', i);
      CondSz(i) = UnrolledSz('c', i);
   }

   // Construct the Hamiltonian
   for (int i = 1; i <= L; ++i)
   {
      if (i != L)  // open boundary conditions
      {
         // hopping
         Hamiltonian += (-t) * (prod(CHup(i), Cup(i%L+1), Ident) + prod(CHdown(i), Cdown(i%L+1), Ident))
            + conj(-t) * (prod(CHup(i%L+1), Cup(i), Ident) + prod(CHdown(i%L+1), Cdown(i), Ident));

         H_t += (-1) * (CHup(i)*Cup(i%L+1) + CHdown(i)*Cdown(i%L+1)
			+ CHup(i%L+1)*Cup(i) + CHdown(i%L+1)*Cdown(i));

         // Ising interaction between localized spins
	 H_Jzh += LocalSz(i) * LocalSz(i%L+1);
         Hamiltonian += Jzh*prod(LocalSz(i), LocalSz(i%L+1), Ident);
      }

      // Kondo coupling
      Hamiltonian += Jk*(prod(CondSz(i), LocalSz(i), Ident)
         + 0.5 * (prod(CondSp(i), LocalSm(i), Ident) + prod(CondSm(i), LocalSp(i), Ident)));

      H_Jk += CondSz(i) * LocalSz(i) + 0.5*(CondSp(i)*LocalSm(i) + CondSm(i)*LocalSp(i));

      // Jz kondo interaction
      H_Jz += CondSz(i) * LocalSz(i);

      std::cout << "Working.... " << i << "\n";
   }

   // The antiferromagnetic order parameter might be useful
   MPOperator& StaggeredfSz = OpList["StaggeredfSz"];
   for (int i = 1; i <= L; ++i)
   {
      StaggeredfSz = minus1pow(i) * LocalSz(i);
   }
   // operator 'AFM' is the square of the localized staggered Sz
   OpList["AFM"] = prod(StaggeredfSz, StaggeredfSz, Ident);

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);
   pheap::ShutdownPersistent(OList);
}
