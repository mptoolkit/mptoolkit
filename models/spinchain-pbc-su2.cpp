// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "models/spin-su2.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      double J = 1;
      double Beta = 0;
      int L = 0;
      half_int Spin = 0.5;
      std::string LatticeName;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("LatticeSize,L", prog_opt::value(&L), "lattice size [required]")
         ("Spin,S", prog_opt::value(&Spin), "magnitude of the spin in the bulk [default 0.5]")
         ("J", prog_opt::value(&J), "nearest neighbor bilinear coupling [default 1]")
         ("Beta", prog_opt::value(&Beta), 
	  "nearest neighbor biquadratic coupling [default 0]")
         ("out,o", prog_opt::value(&LatticeName), "output filename [required]")
         ;
      
      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).style(prog_opt::command_line_style::default_style ^
					  prog_opt::command_line_style::allow_guessing).
		      run(), vm);
      prog_opt::notify(vm);    
      
      if (vm.count("help") || !vm.count("LatticeSize") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: spinchain-su2 [options]\n";
         std::cerr << desc << '\n';
	 std::cerr << "\nAlternatively, use shortcut operators\n"
		   << "H1 - bilinear term\n"
		   << "H2 - biquadratic term\n"
		   << "\nH = J*H1 + Neta*H2\n";
         return 1;
      }

      TRACE(Spin)(J)(Beta);

      // Construct the site block
      SiteBlock Site = CreateSU2SpinSite(Spin);

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

      OperatorAtSite<OperatorList const, int> S(OpList, "S");
      //      OperatorAtSite<OperatorList, int> Bond(OpList, "Bond");
      
      MPOperator& Hamiltonian = OpList["H"];

      // Split operators.  H1 is the bilinear term, H2 is biquadratic
      MPOperator& H1 = OpList["H1"];
      MPOperator& H2 = OpList["H2"];

      QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
      // interaction matrix elements
      std::cout << "Working" << std::flush;
      for (int i = 1; i <= L; ++i)
      {
	 MPOperator SS = -sqrt(3.0) * prod(S(i), S(i%L+1), Ident);
         H1 += SS;
         H2 += prod(SS, SS, Ident);
         std::cout << '.' << std::flush;
      }
      std::cout << "done!" << std::endl;

      Hamiltonian = J*H1 + Beta*H2;

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
