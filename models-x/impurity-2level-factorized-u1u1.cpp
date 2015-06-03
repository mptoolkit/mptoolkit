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
   if (argc != 11)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: impurity-2level-factorized-u1u1 <Bandwidth> <Gamma1> "
         "<EpsilonD1> <Gamma2> <EpsilonD2> <Delta> <Ec> <N> <H> <outfile> < couplings\n"
         //                << "L = number of lattice sites\n"
         //                << "t = hopping integral\n"
                << "Delta = impurity hopping\n"
                << "Ec, N: energy term Ec*(n_d - N)^2"
                << "H = local magnetic field\n"
                << "outfile = file name for output lattice\n"
                << "couplings = epsilon_0, t_i epsilon_i for i = 1 .. L.\n";
      return 1;
   }

   double Bandwidth = boost::lexical_cast<double>(argv[1]);
   double Gamma1 =  boost::lexical_cast<double>(argv[2]);
   double EpsilonD1 =  boost::lexical_cast<double>(argv[3]);
   double Gamma2 =  boost::lexical_cast<double>(argv[4]);
   double EpsilonD2 =  boost::lexical_cast<double>(argv[5]);

   double U = boost::lexical_cast<double>(argv[1]);
   double H0 = boost::lexical_cast<double>(argv[2]);

   std::vector<double> Epsilon;
   std::vector<double> Hop;

   double e,t;
   while (std::cin >> t >> e)
   {
      Hop.push_back(t);
      Epsilon.push_back(e);
   }

   int L = Hop.size();
   std::cout << "Lattice size (excluding impurity) = " << L << '\n';
   
   // Construct the site block
   SiteBlock DownSite = FactorizedSite(1, -0.5);
   SiteBlock UpSite = FactorizedSite(1, 0.5);
   SiteBlock ImpuritySite = CreateU1HubbardSite();

   // Create the lattice
   Lattice MyLattice;

   // Starting with the down spins
   for (int i = 1; i <= L; ++i)
   {
      MyLattice.Append("down", L-i+1, DownSite);
   }
   // add the impurity sites
   MyLattice.Append("impurity1", ImpuritySite);
   MyLattice.Append("impurity2", ImpuritySite);
   // add the up spins
   for (int i = 1; i <= L; ++i)
   {
      MyLattice.Append("up", i, UpSite);
   }

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number

   // operators on the unwrapped lattice
   OperatorAtSite<OperatorList const, string, int> C(OpList, "C");
   OperatorAtSite<OperatorList const, string, int> CH(OpList, "CH");
   OperatorAtSite<OperatorList const, string, int> LocalSz(OpList, "Sz");
   OperatorAtSite<OperatorList const, string, int> LocalN(OpList, "N");
   
   // make operators for the 'wrapped' lattice
   OperatorAtSite<OperatorList, int> CHup(OpList, "CHup");
   OperatorAtSite<OperatorList, int> Cup(OpList, "Cup");
   OperatorAtSite<OperatorList, int> CHdown(OpList, "CHdown");
   OperatorAtSite<OperatorList, int> Cdown(OpList, "Cdown");
   OperatorAtSite<OperatorList, int> N(OpList, "N");
   OperatorAtSite<OperatorList, int> Pdouble(OpList, "Pdouble");
   OperatorAtSite<OperatorList, int> Sz(OpList, "Sz");
   OperatorAtSite<OperatorList, int> Sp(OpList, "Sp");
   OperatorAtSite<OperatorList, int> Sm(OpList, "Sm");

   for (int i = 1; i <= L; ++i)
   {
      Cup(i) = C("up", i);
      CHup(i) = CH("up", i);
      Cdown(i) = C("down", i);
      CHdown(i) = CH("down", i);
      N(i) = LocalN("up", i) + LocalN("down", i);
      Pdouble(i) = prod(LocalN("up", i), LocalN("down", i), Ident);
      Sz(i) = LocalSz("up", i) + LocalSz("down", i);
      Sp(i) = prod(CH("up", i), C("down", i), 
                   QuantumNumber(MyLattice.GetSymmetryList(), "0,1"));
      Sm(i) = prod(CH("down", i), C("up", i),
                   QuantumNumber(MyLattice.GetSymmetryList(), "0,-1"));
      std::cout << "Working.... " << i << "\n";
   }
   // Sites -1 and -2 are the impurities
   Cup(-1) = OpList["Cup(impurity1)"];
   CHup(-1) = OpList["CHup(impurity1)"];
   Cdown(-1) = OpList["Cdown(impurity1)"];
   CHdown(-1) = OpList["CHdown(impurity1)"];
   N(-1) = OpList["N(impurity1)"];
   Pdouble(-1) = OpList["Pdouble(impurity1)"];
   Sz(-1) = OpList["Sz(impurity1)"];
   Sp(-1) = OpList["Sp(impurity1)"];
   Sm(-1) = OpList["Sm(impurity1)"];

   Cup(-2) = OpList["Cup(impurity2)"];
   CHup(-2) = OpList["CHup(impurity2)"];
   Cdown(-2) = OpList["Cdown(impurity2)"];
   CHdown(-2) = OpList["CHdown(impurity2)"];
   N(-2) = OpList["N(impurity2)"];
   Pdouble(-2) = OpList["Pdouble(impurity2)"];
   Sz(-2) = OpList["Sz(impurity2)"];
   Sp(-2) = OpList["Sp(impurity2)"];
   Sm(-2) = OpList["Sm(impurity2)"];

   // Now we can refer to the operators as if they were on the original lattice

   MPOperator& Hamiltonian = OpList["H"];
   // The bath
   for (int i = 1; i <= L; ++i)
   {
      MPOperator Hopping 
         = prod(CHup(i), Cup(i%L+1), Ident)
      	 - prod(Cup(i), CHup(i%L+1), Ident)
         + prod(CHdown(i), Cdown(i%L+1), Ident)
      	 - prod(Cdown(i), CHdown(i%L+1), Ident);
      Hamiltonian += -Hop[i] * Hopping + Epsilon[i] * N(i);

      std::cout << "Working.... " << i << "\n";
   }

   // Impurity 1, energy, hopping
   Hamiltonian += 

   // the final epsilon term
   Hamiltonian = Hamiltonian + Epsilon[L] * N(L);

   // coulomb term
   Hamiltonian = Hamiltonian + U * OpList["Pdouble(0)"];

   // local field
   Hamiltonian = Hamiltonian + H0 * OpList["Sz(0)"];

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(argv[3], 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
