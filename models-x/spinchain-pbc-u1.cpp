// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/spin-u1.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      double J = 1;
      double Jp = 0;
      double Jz = 0;
      double Beta = 0;
      int L = 0;
      half_int Spin = 0.5;
      std::string LatticeName;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("LatticeSize,L", prog_opt::value(&L), "lattice size [required]")
         ("Spin,S", prog_opt::value(&Spin), "magnitude of the spin in the bulk [default 0.5]")
         ("J", prog_opt::value(&J), "isotropic bilinear coupling [default 1]")
         ("Jp", prog_opt::value(&Jp), "XY bilinear coupling [default J]")
         ("Jz", prog_opt::value(&Jz), "Z bilinear coupling [default J]")
         ("Beta", prog_opt::value(&Beta), "nearest neighbor biquadratic coupling [default 0]")
         ("out,o", prog_opt::value(&LatticeName), "output filename [required]")
         ;
      
      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).run(), vm);
      prog_opt::notify(vm);    
      
      if (vm.count("help") || !vm.count("LatticeSize") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: spinchain-su2 [options]\n";
         std::cerr << desc << '\n';
	 std::cerr << "Alternatively, use the separated components:\n"
		   << "H1p - XY bilinear spin exchange\n"
		   << "H1z - Z bilinear spin exchange\n"
		   << "H2  = biquadratic spin exchange\n"
		   << "\nOr the shortcuts:\n"
		   << "H1 = Jp*H1p + Jz*H1z\n"
		   << "H = H1 + Beta*H2\n";
         return 1;
      }

      // Set the dependent defaults if necessary
      if (!vm.count("Jp"))        Jp = J;
      if (!vm.count("Jz"))        Jz = J;

      TRACE(Spin)(Jp)(Jz)(Beta);
   
      // Construct the site block
      SiteBlock Site = CreateU1SpinSite(Spin);

      // construct a lattice of L copies of Site
      Lattice MyLattice = repeat(Site, L);
     // The optimal numbering of the sites is different for PBC.
      // We want to interlace them as
      // 1  L  2  L-1  3  L-2  ...  L/2  L/2+1
      std::list<std::string> Coords;
      for (int i = 1; i <= L/2; ++i)
      {
	 Coords.push_back(boost::lexical_cast<std::string>(i));
	 Coords.push_back(boost::lexical_cast<std::string>(L-i+1));
      }
      if (L%2 == 1) // if the lattice size is odd, add in the final site
	 Coords.push_back(boost::lexical_cast<std::string>(L/2+1));

      MyLattice.fix_coordinates_from_sequence(Coords.begin(), Coords.end());

      // construct the operator list for the lattice
      OperatorList OpList(MyLattice);
      
      OperatorAtSite<OperatorList const, int> Sp(OpList, "Sp");
      OperatorAtSite<OperatorList const, int> Sm(OpList, "Sm");
      OperatorAtSite<OperatorList const, int> Sz(OpList, "Sz");
      MPOperator& Hamiltonian = OpList["H"];
      MPOperator& TotalSp = OpList["Sp"];
      MPOperator& TotalSm = OpList["Sm"];
      MPOperator& TotalSz = OpList["Sz"];
      MPOperator& TotalS2 = OpList["S2"];
      // Split operators.  H1 is the bilinear term, H2 is biquadratic
      MPOperator& H1p = OpList["H1p"];
      MPOperator& H1z = OpList["H1z"];
      MPOperator& H1 = OpList["H1"];
      MPOperator& H2 = OpList["H2"];

      QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
      // interaction matrix elements
      std::cout << "Working" << std::flush;
      for (int i = 1; i <= L; ++i)
      {
	 MPOperator SSperp = 0.5 * (Sp(i)*Sm(i%L+1) + Sm(i)*Sp(i%L+1));
	 MPOperator SSz = Sz(i)*Sz(i%L+1);
	 MPOperator SS2 = (SSperp+SSz)*(SSperp+SSz);

	 H1p += SSperp;
	 H1z += SSz;
	 H2 += SS2;
         std::cout << '.' << std::flush;
      }
      std::cout << "done!" << std::endl;
 
      H1 = Jp*H1p + Jz*H1z;
      Hamiltonian = H1 + Beta*H2;

      // The total spin operators
      for (int i = 1; i <= L; ++i)
      {
	 TotalSp += Sp(i);
	 TotalSm += Sm(i);
	 TotalSz += Sz(i);
      }

      TotalS2 = 0.5 * (prod(TotalSp, TotalSm, Ident) + prod(TotalSm, TotalSp, Ident))
	 + prod(TotalSz, TotalSz, Ident);

      // make a copy of OpList that exists on the persistent heap
      pvalue_ptr<OperatorList> OList = new OperatorList(OpList);

      pheap::Initialize(LatticeName, 1, 65536, 655360);
      pheap::ShutdownPersistent(OList);
   }
   catch (std::exception& e)
   {
      std::cerr << "Exception: " << e.what() << '\n';
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!\n";
      return 1;
   }
}
