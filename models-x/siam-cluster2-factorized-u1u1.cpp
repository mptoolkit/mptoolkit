// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/factorized-u1u1.h"
#include "models/hubbard-u1u1.h"

using std::string;

int main(int argc, char** argv)
{
   if (argc != 4)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: siam-cluster2-factorized-u1u1 <U> <J> <outfile> < couplings\n"
         //                << "L = number of lattice sites\n"
         //                << "t = hopping integral\n"
                << "U = impurity U for impurity 1,2\n"
                << "J = inter-impurity J\n"
                << "outfile = file name for output lattice\n"
                << "couplings = tperp_0 epsilon_0, t_i tcross_i tperp_i epsilon_i for i = 1 .. L.\n";
      return 1;
   }

   double U = boost::lexical_cast<double>(argv[1]);
   double J = boost::lexical_cast<double>(argv[2]);

   std::vector<double> Epsilon;
   std::vector<double> Hop, HopPerp, HopCross;

   double tp, e,t,tc;
   std::cin >> tp >> e;
   HopPerp.push_back(tp);
   Epsilon.push_back(e);
   while (std::cin >> t >> tc >> tp >> e)
   {
      Hop.push_back(t);
      HopCross.push_back(tc);
      HopPerp.push_back(tp);
      Epsilon.push_back(e);
   }

   int L = Hop.size();
   std::cout << "Lattice size (excluding impurity) = " << L << '\n';
   
   // Construct the site block
   SiteBlock DownSite = FactorizedSite(1, -0.5);
   SiteBlock UpSite = FactorizedSite(1, 0.5);
   SiteBlock ImpuritySite = CreateU1HubbardSite();

   // Create the lattice
   Lattice DownL(DownSite);
   DownL.fix_coordinates("down");
   DownL = repeat(DownL, 2);
   DownL.fix_coordinates_prepend();

   Lattice UpL(UpSite);
   UpL.fix_coordinates("up");
   UpL = repeat(UpL, 2);
   UpL.fix_coordinates_prepend();

   Lattice ImpurityL(ImpuritySite);
   ImpurityL.fix_coordinates("impurity");
   ImpurityL = repeat(ImpurityL, 2);
   ImpurityL.fix_coordinates_prepend();

   Lattice DownLattice = repeat(DownL, L);
   DownLattice.fix_coordinates();
   Lattice UpLattice = repeat(UpL, L);
   UpLattice.fix_coordinates();

   Lattice MyLattice = join(DownLattice, ImpurityL, UpLattice);
   MyLattice.fix_coordinates_unique();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number

   // operators on the unwrapped lattice
   OperatorAtSite<OperatorList const, int, string, int> C(OpList, "C");
   OperatorAtSite<OperatorList const, int, string, int> CH(OpList, "CH");
   OperatorAtSite<OperatorList const, int, string, int> LocalSz(OpList, "Sz");
   OperatorAtSite<OperatorList const, int, string, int> LocalN(OpList, "N");
   
   // make operators for the 'wrapped' lattice
   OperatorAtSite<OperatorList, int, int> CHup(OpList, "CHup");
   OperatorAtSite<OperatorList, int, int> Cup(OpList, "Cup");
   OperatorAtSite<OperatorList, int, int> CHdown(OpList, "CHdown");
   OperatorAtSite<OperatorList, int, int> Cdown(OpList, "Cdown");
   OperatorAtSite<OperatorList, int, int> N(OpList, "N");
   OperatorAtSite<OperatorList, int, int> Pdouble(OpList, "Pdouble");
   OperatorAtSite<OperatorList, int, int> Sz(OpList, "Sz");
   OperatorAtSite<OperatorList, int, int> Sp(OpList, "Sp");
   OperatorAtSite<OperatorList, int, int> Sm(OpList, "Sm");

   for (int b = 1; b < 3; ++b)
   {
      for (int i = 1; i <= L; ++i)
      {
         Cup(b,i) = C(b, "up", i);
         CHup(b,i) = CH(b, "up", i);
         Cdown(b,i) = C(b, "down", i);
         CHdown(b,i) = CH(b, "down", i);
         N(b,i) = LocalN(b, "up", i) + LocalN(b, "down", i);
         Pdouble(b,i) = prod(LocalN(b, "up", i), LocalN(b, "down", i), Ident);
         Sz(b,i) = LocalSz(b, "up", i) + LocalSz(b, "down", i);
         Sp(b,i) = prod(CH(b, "up", i), C(b, "down", i), 
                      QuantumNumber(MyLattice.GetSymmetryList(), "0,1"));
         Sm(b,i) = prod(CH(b, "down", i), C(b, "up", i),
                      QuantumNumber(MyLattice.GetSymmetryList(), "0,-1"));
         std::cout << "Working.... " << i << "\n";
      }
      // Site 0 is the impurity site
      Cup(b,0) = OpList["Cup("+boost::lexical_cast<std::string>(b)+",impurity)"];
      CHup(b,0) = OpList["CHup("+boost::lexical_cast<std::string>(b)+",impurity)"];
      Cdown(b,0) = OpList["Cdown("+boost::lexical_cast<std::string>(b)+",impurity)"];
      CHdown(b,0) = OpList["CHdown("+boost::lexical_cast<std::string>(b)+",impurity)"];
      N(b,0) = OpList["N("+boost::lexical_cast<std::string>(b)+",impurity)"];
      Pdouble(b,0) = OpList["Pdouble("+boost::lexical_cast<std::string>(b)+",impurity)"];
      Sz(b,0) = OpList["Sz("+boost::lexical_cast<std::string>(b)+",impurity)"];
      Sp(b,0) = OpList["Sp("+boost::lexical_cast<std::string>(b)+",impurity)"];
      Sm(b,0) = OpList["Sm("+boost::lexical_cast<std::string>(b)+",impurity)"];
   }

   // Now we can refer to the operators as if they were on the original lattice

   MPOperator& Hamiltonian = OpList["H"];
   // hopping matrix elements
   for (int i = 0; i < L; ++i)
   {
      for (int b = 1; b <= 2; ++b)
      {
         MPOperator Hopping 
            = prod(CHup(b,i), Cup(b,i%L+1), Ident)
            - prod(Cup(b,i), CHup(b,i%L+1), Ident)
            + prod(CHdown(b,i), Cdown(b,i%L+1), Ident)
            - prod(Cdown(b,i), CHdown(b,i%L+1), Ident);
         Hamiltonian += -Hop[i] * Hopping + Epsilon[i] * N(b,i);
      }

      MPOperator CrossTerm = prod(CHup(1,i), Cup(2,i+1), Ident)
         + prod(CHdown(1,i), Cdown(2,i+1), Ident)
         - prod(Cup(1,i), CHup(2,i+1), Ident)
         - prod(Cdown(1,i), CHdown(2,i+1), Ident)
         + prod(CHup(2,i), Cup(1,i+1), Ident)
         + prod(CHdown(2,i), Cdown(1,i+1), Ident)
         - prod(Cup(2,i), CHup(1,i+1), Ident)
         - prod(Cdown(2,i), CHdown(1,i+1), Ident);

      Hamiltonian += -HopCross[i] * CrossTerm;

      MPOperator PerpTerm = prod(CHup(1,i), Cup(2,i), Ident)
         + prod(CHdown(1,i), Cdown(2,i), Ident)
         - prod(Cup(1,i), CHup(2,i), Ident)
         - prod(Cdown(1,i), CHdown(2,i), Ident);

      Hamiltonian += -HopPerp[i] * PerpTerm;

      std::cout << "Working.... " << i << "\n";
   }
   // the final epsilon term
   for (int b = 1; b <= 2; ++b)
   {
      Hamiltonian = Hamiltonian + Epsilon[L] * N(b,L);
   }

   // Interchange the bands
   

   // the final HopPerp term

   MPOperator PerpTerm = prod(CHup(1,L), Cup(2,L), Ident)
      + prod(CHdown(1,L), Cdown(2,L), Ident)
      - prod(Cup(1,L), CHup(2,L), Ident)
      - prod(Cdown(1,L), CHdown(2,L), Ident);
   Hamiltonian += -HopPerp[L] * PerpTerm;

   // coulomb term
   Hamiltonian = Hamiltonian + U * (OpList["Pdouble(1,0)"] + OpList["Pdouble(2,0)"]);

   // inter-impurity spin coupling
   Hamiltonian = Hamiltonian + J * (prod(Sz(1,0), Sz(2,0), Ident)
                                    + 0.5 * (prod(Sp(1,0), Sm(2,0), Ident)
                                             + prod(Sm(1,0), Sp(2,0), Ident)));

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(argv[3], 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
