// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/bosehubbard-2bosons-u1u1.h"

int main(int argc, char** argv)
{
   if (argc != 13)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: bosehubbard-u1u1 <MaxN> <L> <tA> <tB> <UA> <UB> <Uint> <mu0_A> <dn_A> <mu0_B> <dn_B> <outfile>\n"
		<< "MaxN - maximum number of bosons per site\n"
		<< "L  - length of the system\n"
		<< "tA - hopping strength for particle type A\n"
		<< "tB - hopping strength for particle type B\n"
		<< "UA - onside interaction for type A\n"
		<< "UB - onside interaction for type B\n"
		<< "Uint - onside interaction between both types \n"
		<< "mu0_+ strength of inital disturbance for type A+B \n" 
		<< "dn0_+ width of inital disturbance for type A+B \n" 
		<< "mu0_- strength of inital disturbance for type A-B \n" 
		<< "dn0_- width of inital disturbance for type A-B \n" 
		<< "outfile = file name for output lattice.\n";
      return 1;
   }

   int MaxN = boost::lexical_cast<int>(argv[1]);
	int L = boost::lexical_cast<int>(argv[2]);
	double tA = boost::lexical_cast<double>(argv[3]);
	double tB = boost::lexical_cast<double>(argv[4]);
   double UA = boost::lexical_cast<double>(argv[5]);
	double UB = boost::lexical_cast<double>(argv[6]);
	double Uint = boost::lexical_cast<double>(argv[7]);
   double mu1 = boost::lexical_cast<double>(argv[8]);
	double dn1 = boost::lexical_cast<double>(argv[9]);
	double mu2 = boost::lexical_cast<double>(argv[10]);
	double dn2 = boost::lexical_cast<double>(argv[11]);
	
	if (dn1 == 0.) {mu1 = 0.; dn1 = 1.;}
	if (dn2 == 0.) {mu2 = 0.; dn2 = 1.;}

   // Construct the site block
   SiteBlock Site = CreateBoseHubbard2BosonsU1U1Site(MaxN);

   // construct a lattice of L copies of Site
   Lattice MyLattice = repeat(Site, L);
   MyLattice.fix_coordinates();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   OperatorAtSite<OperatorList const, int> BH_A(OpList, "BH_A");
   OperatorAtSite<OperatorList const, int> B_A(OpList, "B_A");
   OperatorAtSite<OperatorList const, int> N2_A(OpList, "N2_A");
   OperatorAtSite<OperatorList const, int> N_A(OpList, "N_A");

	OperatorAtSite<OperatorList const, int> BH_B(OpList, "BH_B");
   OperatorAtSite<OperatorList const, int> B_B(OpList, "B_B");
   OperatorAtSite<OperatorList const, int> N2_B(OpList, "N2_B");
   OperatorAtSite<OperatorList const, int> N_B(OpList, "N_B");

   QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
   MPOperator Hamiltonian;
	
	// hopping matrix elements
   for (int i = 1; i < L; ++i)
   {
      Hamiltonian += -1.0 * tA * ( prod(BH_A(i), B_A(i+1), Ident) + prod(B_A(i), BH_A(i+1), Ident) );
      
		Hamiltonian += -1.0 * tB * ( prod(BH_B(i), B_B(i+1), Ident) + prod(B_B(i), BH_B(i+1), Ident) );
   }
	
	// Onside interaction
	for (int i = 1; i <= L; ++i)
   {
      Hamiltonian += (UA/2.0) * N2_A(i);
		Hamiltonian += (UB/2.0) * N2_B(i);
		Hamiltonian += Uint * prod(N_A(i), N_B(i), Ident);
		
   }
	
	OpList["H0"] = Hamiltonian;
	
	// chemical potential
	
	if ((mu1 != 0.) && (mu2 != 0.)) {
		
		for (int i = 1; i <= L; ++i)
   	{
	      Hamiltonian += ( mu1 * exp(- double(2*i - 1 - L) * double(2*i - 1 - L) / ( dn1 * dn1 * 2.) ) +  
			                 mu2 * exp(- double(2*i - 1 - L) * double(2*i - 1 - L) / ( dn2 * dn2 * 2.) )   ) * N_A(i);
		
			Hamiltonian += ( mu1 * exp(- double(2*i - 1 - L) * double(2*i - 1 - L) / ( dn1 * dn1 * 2.) ) -  
		   	              mu2 * exp(- double(2*i - 1 - L) * double(2*i - 1 - L) / ( dn2 * dn2 * 2.) )   ) * N_B(i);
   	}
	}
	
	
   // save the Hamiltonian into the operator list
   OpList["H"] = Hamiltonian;

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(argv[12], 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
