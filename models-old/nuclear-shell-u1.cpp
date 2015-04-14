// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/fermion-m-basis.h"
#include <fstream>

int main(int argc, char** argv)
{
   if (argc != 3)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: nuclear-shell-u1 <4-body-matrix-file> <outfile>\n"
                << "4-body-matrix-file = matrix elements V(i,j,k,l)\n"
                << "outfile = file name for output lattice.\n";
      return 1;
   }

   std::string VFileName = argv[1];
   std::string OutFileName = argv[2];
   
   // Construct the site blocks
   std::map<half_int, SiteBlock> ProtonSites, NeutronSites;
   ProtonSites[0.5] = CreateFermionSite(0.5, "Np", "Sz");
   ProtonSites[1.5] = CreateFermionSite(1.5, "Np", "Sz");
   ProtonSites[2.5] = CreateFermionSite(2.5, "Np", "Sz");
   ProtonSites[3.5] = CreateFermionSite(3.5, "Np", "Sz");
   NeutronSites[0.5] = CreateFermionSite(0.5, "Nn", "Sz");
   NeutronSites[1.5] = CreateFermionSite(1.5, "Nn", "Sz");
   NeutronSites[2.5] = CreateFermionSite(2.5, "Nn", "Sz");
   NeutronSites[3.5] = CreateFermionSite(3.5, "Nn", "Sz");

   SymmetryList QL("Nn:U(1),Np:U(1),Sz:U(1)");
   ProtonSites[0.5].CoerceSymmetryList(QL);
   ProtonSites[1.5].CoerceSymmetryList(QL);
   ProtonSites[2.5].CoerceSymmetryList(QL);
   ProtonSites[3.5].CoerceSymmetryList(QL);
   NeutronSites[0.5].CoerceSymmetryList(QL);
   NeutronSites[1.5].CoerceSymmetryList(QL);
   NeutronSites[2.5].CoerceSymmetryList(QL);
   NeutronSites[3.5].CoerceSymmetryList(QL);

   // Join a proton and a neutron site of the same z-spin
   std::map<half_int, Lattice> zSpinLattice;
   for (half_int i = 0.5; i <= 3.5; ++i)
   {
      zSpinLattice[i] = Lattice(ProtonSites[i], NeutronSites[i]);
      zSpinLattice[i].fix_coordinates();  // '1' for proton, '2' for neutron
   }

   // Join z-spin lattices together to make j-spin site
   std::map<half_int, Lattice> jSpinLattice;
   jSpinLattice[0.5] = zSpinLattice[0.5];    // spin 1/2 case is easy
   for (half_int i = 1.5; i <= 3.5; ++i)
   {
      jSpinLattice[i] = join(jSpinLattice[i-1], zSpinLattice[i]);
   }
   // fix coordinates of the components of the j-site
   // coordinates are 1 for spin 1/2, 2 for spin 3/2 etc
   for (half_int i = 0.5; i <= 3.5; ++i)
      jSpinLattice[i].fix_coordinates();

   // Now assemble the final lattice
   Lattice MyLattice = join(jSpinLattice[3.5], jSpinLattice[1.5],
			    jSpinLattice[0.5], jSpinLattice[2.5]);
   MyLattice.fix_coordinates();

   // construct the operator list for the lattice
   OperatorList OpList(MyLattice);

   OperatorAtSite<OperatorList const> CHup(OpList, "CHup");
   OperatorAtSite<OperatorList const> Cup(OpList, "Cup");
   OperatorAtSite<OperatorList const> CHdown(OpList, "CHdown");
   OperatorAtSite<OperatorList const> Cdown(OpList, "Cdown");
   OperatorAtSite<OperatorList const> N(OpList, "N");
   MPOperator& Hamiltonian = OpList["H"];

   // Assemble an array of operators that matches the convention
   // in the V-file

   // major index 1 = proton, 2 = neutron
   std::vector<std::map<int, MPOperator> > CH(3), C(3);

   for (int i = 1; i <= 2; ++i)
   {
      // Spin 7/2 shell
      CH[i][1] = CHdown(i,1,1);
      CH[i][2] = CHup(i,1,1);
      CH[i][3] = CHdown(i,2,1);
      CH[i][4] = CHup(i,2,1);
      CH[i][5] = CHdown(i,3,1);
      CH[i][6] = CHup(i,3,1);
      CH[i][7] = CHdown(i,4,1);
      CH[i][8] = CHup(i,4,1);
      // Spin 3/2 shell
      CH[i][9] = CHdown(i,1,2);
      CH[i][10] = CHup(i,1,2);
      CH[i][11] = CHdown(i,2,2);
      CH[i][12] = CHup(i,2,2);
      // Spin 1/2 shell
      CH[i][13] = CHdown(i,1,3);
      CH[i][14] = CHup(i,1,3);
      // Spin 5/2 shell
      CH[i][15] = CHdown(i,1,4);
      CH[i][16] = CHup(i,1,4);
      CH[i][17] = CHdown(i,2,3);
      CH[i][18] = CHup(i,2,3);
      CH[i][19] = CHdown(i,3,2);
      CH[i][20] = CHup(i,3,2);
   }

   for (int i = 1; i <= 2; ++i)
   {
      // Spin 7/2 shell
      C[i][1] = Cdown(i,1,1);
      C[i][2] = Cup(i,1,1);
      C[i][3] = Cdown(i,2,1);
      C[i][4] = Cup(i,2,1);
      C[i][5] = Cdown(i,3,1);
      C[i][6] = Cup(i,3,1);
      C[i][7] = Cdown(i,4,1);
      C[i][8] = Cup(i,4,1);
      // Spin 3/2 shell
      C[i][9] = Cdown(i,1,2);
      C[i][10] = Cup(i,1,2);
      C[i][11] = Cdown(i,2,2);
      C[i][12] = Cup(i,2,2);
      // Spin 1/2 shell
      C[i][13] = Cdown(i,1,3);
      C[i][14] = Cup(i,1,3);
      // Spin 5/2 shell
      C[i][15] = Cdown(i,1,4);
      C[i][16] = Cup(i,1,4);
      C[i][17] = Cdown(i,2,3);
      C[i][18] = Cup(i,2,3);
      C[i][19] = Cdown(i,3,2);
      C[i][20] = Cup(i,3,2);
   }

   // read the matrix elements for the 4-body term
   std::vector<std::vector<std::vector<std::vector<double> > > > 
      V1(21, std::vector<std::vector<std::vector<double> > >(21, std::vector<std::vector<double> >(21, std::vector<double>(21, 0.0))));

   std::vector<std::vector<std::vector<std::vector<double> > > > 
      V2(21, std::vector<std::vector<std::vector<double> > >(21, std::vector<std::vector<double> >(21, std::vector<double>(21, 0.0))));

   std::vector<std::vector<std::vector<std::vector<double> > > > 
      V3(21, std::vector<std::vector<std::vector<double> > >(21, std::vector<std::vector<double> >(21, std::vector<double>(21, 0.0))));

   {
      std::cout << "Reading file..." << std::endl;
      std::ifstream VFile(VFileName.c_str());
      int Which, i,j,k,l;
      double x;
      while (VFile >> Which >> i >> j >> k >> l >> x)
      {

	 if (Which == 1) // neutrons
	 {
	    V1[i][j][k][l] = x;
	    //	 Hamiltonian += x * CH[2][i] * CH[2][j] * C[2][l] * C[2][k];
	 }
	 else if (Which == 2) // protons
	 {
	    V2[i][j][k][l] = x;
	    // Hamiltonian += x * CH[1][i] * CH[1][j] * C[1][l] * C[1][k];
	 }
	 else if (Which == 3) // proton/neutron interaction
	 {
	    V3[i][j][k][l] = x;
	    //Hamiltonian += x * CH[1][i] * CH[2][j] * C[1][l] * C[2][k];
	 }
      }
      std::cout << "done.\n";
   }

   for (int i = 1; i <= 20; ++i)
   {
      MPOperator iAccum1, iAccum2, iAccum3;
      for (int j = 1; j <= 20; ++j)
      {
	 MPOperator jAccum1, jAccum2, jAccum3;
	 for (int k = 1; k <= 20; ++k)
	 {
	    std::cout << "Adding matrix element ("<<i<<","<<j<<","<<k<<",*)\n";
	    double* V1i = &V1[i][j][k][0];
	    double* V2i = &V2[i][j][k][0];
	    double* V3i = &V3[i][j][k][0];
	    MPOperator kAccum1, kAccum2, kAccum3;
	    for (int l = 1; l <= 20; ++l)
	    {
	       if (V1i[l] != 0)
		  kAccum1 += V1i[l] * C[2][l];
	       if (V2i[l] != 0)
		  kAccum2 += V2i[l] * C[1][l];
	       if (V3i[l] != 0)
		  kAccum3 += V3i[l] * C[1][l];
	    }
	    std::cout << "Adding matrix element ("<<i<<","<<j<<","<<k<<")\n";
	    jAccum1 += kAccum1 * C[2][k];
	    jAccum2 += kAccum2 * C[1][k];
	    jAccum3 += kAccum3 * C[2][k];
	 }
	 std::cout << "Adding matrix element ("<<i<<","<<j<<")\n";
	 iAccum1 += CH[2][j] * jAccum1;
	 iAccum2 += CH[1][j] * jAccum2;
	 iAccum3 += CH[2][j] * jAccum3;
      }
      std::cout << "Adding matrix element ("<<i<<")\n";
      Hamiltonian += CH[2][i] * iAccum1;
      Hamiltonian += CH[1][i] * iAccum2;
      Hamiltonian += CH[1][i] * iAccum3;
   }

   // 2-body term
   for (int n = 1; n <= 2; ++n)
   {
      Hamiltonian += 2.0 * N(n,1,2);
      Hamiltonian += 2.0 * N(n,2,2);
      Hamiltonian += 4.0 * N(n,1,3);
      Hamiltonian += 6.5 * N(n,1,4);
      Hamiltonian += 6.5 * N(n,2,3);
      Hamiltonian += 6.5 * N(n,3,2);
   }

   // Project the hamiltonian onto the scalar component, but we want to check
   // the remainder term too
   MPOperator& HRem = OpList["HRem"];

   HRem = Hamiltonian;
   Hamiltonian = project(Hamiltonian, QuantumNumber(QL, "0,0,0"));

   // make a copy of OpList that exists on the persistent heap
   pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

   pheap::Initialize(OutFileName, 1, 65536, 655360);
   pheap::ShutdownPersistent(OList);
}
