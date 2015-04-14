// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/bosehubbard-spinless-u1.h"

int main(int argc, char** argv)
{
   if (argc != 9)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: bosehubbard-u1u1 <MaxN> <L> <t> <U> <Uint> <mu> <dn> <outfile>\n"
      << "MaxN - maximum number of bosons per site\n"
      << "L  - length of the system\n"
      << "t - hopping strength\n"
      << "U - intratype onside interaction\n"
      << "Uint - intertype onside interaction\n"
      << "mu strength of inital disturbance \n" 
      << "dn width of inital disturbance \n" 
      << "outfile = file name for output lattice.\n";
      return 1;
   }

   int MaxN = boost::lexical_cast<int>(argv[1]);
   int L = boost::lexical_cast<int>(argv[2]);
   double t = boost::lexical_cast<double>(argv[3]);
   double U = boost::lexical_cast<double>(argv[4]);
   double Uint = boost::lexical_cast<double>(argv[5]);
   double mu = boost::lexical_cast<double>(argv[6]);
   double dn = boost::lexical_cast<double>(argv[7]);
   
   if (dn == 0.) {mu = 0.; dn = 1.;}

   //TRACE(MaxN)(L)(t)(U)(Uint)(mu)(dn);
   
   // Construct the site block
   SiteBlock SiteA = CreateBoseHubbardSpinlessU1Site(MaxN,"NA");
   SiteBlock SiteB = CreateBoseHubbardSpinlessU1Site(MaxN,"NB");

   // Construct the unit cell
   Lattice UnitCell(SymmetryList("NA:U(1),NB:U(1)"), SiteA, SiteB);
   UnitCell.fix_coordinates("0", "1");

   // construct a lattice of L copies of the unit cell
   Lattice MyLattice = repeat(UnitCell, L);
   MyLattice.fix_coordinates_prepend();
   
   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   OperatorAtSite<OperatorList const, int, int> BH(OpList, "BH");
   OperatorAtSite<OperatorList const, int, int> B(OpList, "B");
   
   OperatorAtSite<OperatorList const, int, int> N2(OpList, "N2");
   OperatorAtSite<OperatorList const, int, int> N(OpList, "N");
      
   // and some convinient operators
   
   OperatorAtSite<OperatorList, int> N_A(OpList, "N_A");
   OperatorAtSite<OperatorList, int> N_B(OpList, "N_B");
   
   OperatorAtSite<OperatorList, int> N_S(OpList, "N_S");
   OperatorAtSite<OperatorList, int> N_C(OpList, "N_C");

   OperatorAtSite<OperatorList, int> BH_A(OpList, "BH_A");
   OperatorAtSite<OperatorList, int> BH_B(OpList, "BH_B");

   OperatorAtSite<OperatorList, int> B_A(OpList, "B_A");
   OperatorAtSite<OperatorList, int> B_B(OpList, "B_B");

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   MPOperator Hamiltonian;
   MPOperator Hsp, Hsm, Hcm, Hcp; // disturbed operators
   
   // hopping matrix elements
   for (int i = 1; i < L; ++i)
   {
      Hamiltonian += -1.0 * t * ( prod(BH(i,0), B(i+1,0), Ident) + prod(B(i,0), BH(i+1,0), Ident) );
      
      Hamiltonian += -1.0 * t * ( prod(BH(i,1), B(i+1,1), Ident) + prod(B(i,1), BH(i+1,1), Ident) );
   }
   
   // Onside interaction
   for (int i = 1; i <= L; ++i)
   {
      Hamiltonian += (U/2.0) * N2(i,0);
      Hamiltonian += (U/2.0) * N2(i,1);
      Hamiltonian += Uint * prod(N(i,0), N(i,1), Ident);
   }
   
   OpList["H"] = Hamiltonian;
   
   // disturbance
   
   if (mu != 0.) {
      
      for (int i = 1; i <= L; ++i)
      {
         Hcp += ( mu * exp(- double(2*i - 1 - L) * double(2*i - 1 - L) / ( dn * dn * 2.) ) ) * ( N(i,0) + N(i,1) );
         Hcm -=  ( mu * exp(- double(2*i - 1 - L) * double(2*i - 1 - L) / ( dn * dn * 2.) ) ) * ( N(i,0) + N(i,1) );
         
         Hsp += ( mu * exp(- double(2*i - 1 - L) * double(2*i - 1 - L) / ( dn * dn * 2.) ) ) * ( N(i,0) - N(i,1) );
         Hsm -= ( mu * exp(- double(2*i - 1 - L) * double(2*i - 1 - L) / ( dn * dn * 2.) ) ) * ( N(i,0) - N(i,1) );
      }
   }
   
   OpList["H_cp"] = Hamiltonian + Hcp;
   OpList["H_cm"] = Hamiltonian + Hcm;
   
   OpList["H_sp"] = Hamiltonian + Hsp;
   OpList["H_sm"] = Hamiltonian + Hsm;
   
   
   // build convinient operators
   
   for (int i = 1; i<= L; ++i) {
      N_A(i) = N(i,0);
      N_B(i) = N(i,1);
      
      N_C(i) = N(i,0) + N(i,1);
      N_S(i) = N(i,0) - N(i,1);
      
      BH_A(i) = BH(i,0);
      BH_B(i) = BH(i,1);

      B_A(i) = B(i,0);
      B_B(i) = B(i,1);
   }
   
   
   // save the Hamiltonian into the operator list
   OpList["H"] = Hamiltonian;

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(argv[8], 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
