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
      double JEdge = 1;
      double BetaEdge = 0;
      //      double J2 = 0;
      //      double Beta2 = 0;
      int L = 0;
      half_int Spin = 0.5;
      half_int SpinEdge = 0;
      std::string LatticeName;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("LatticeSize,L", prog_opt::value(&L), "lattice size [required]")
         ("Spin,S", prog_opt::value(&Spin), "magnitude of the spin in the bulk [default 0.5]")
         ("SpinEdge", prog_opt::value(&SpinEdge), "magnitude of the spin at the edges [default Spin]")
         ("J", prog_opt::value(&J), "nearest neighbor bilinear coupling [default 1]")
         ("JEdge", prog_opt::value(&JEdge), "bilinear coupling at the edge bonds [default J]")
         ("Beta", prog_opt::value(&Beta), 
	  "nearest neighbor biquadratic coupling [default 0]")
         ("BetaEdge", prog_opt::value(&BetaEdge), "biquadratic coupling at the edge bonds [default Beta]")
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
	 std::cerr << "Alternatively, use the separated components\n"
		   << "H1Bulk - bilinear spin exchange for the bulk\n"
		   << "H1Edge - bilinear spin exchange at the edges\n"
		   << "H2Bulk - biquadratic spin exchange for the bulk\n"
		   << "H2Edge - biquadratic spin exchange at the edges\n"
		   << "H1p    - bilinear spin exchange across the boundary (perioidic)\n"
		   << "H2p    - biquadratic spin exchange across the boundary (perioidic)\n"
		   << "H1px   - bilinear spin exchange across the boundary ignoring edge sites\n"
		   << "H2px   - biquadratic spin exchange across the boundary ignoring edge sites\n"
		   << "\nOr the shortcuts\n"
		   << "H1 = JEdge*H1Edge + J*H1Bulk\n"
		   << "H2 = BetaEdge*H2Edge + Beta*H2Bulk\n"
		   << "H = H1 + H2\n";
         return 1;
      }

      // Set the dependent defaults if necessary
      if (!vm.count("SpinEdge"))
         SpinEdge = Spin;
      if (!vm.count("JEdge"))
         JEdge = J;
      if (!vm.count("BetaEdge"))
         BetaEdge = Beta;

      TRACE(Spin)(SpinEdge)(J)(JEdge)(Beta)(BetaEdge);

      // Construct the site block
      SiteBlock EdgeSite = CreateSU2SpinSite(SpinEdge);
      SiteBlock BulkSite = CreateSU2SpinSite(Spin);

      // construct a lattice of L copies of Site
      Lattice MyLattice = join(EdgeSite, repeat(BulkSite, L-2), EdgeSite);
      MyLattice.fix_coordinates();

      // construct the operator list for the lattice
      OperatorList OpList(MyLattice);

      OperatorAtSite<OperatorList const, int> S(OpList, "S");
      //      OperatorAtSite<OperatorList, int> Bond(OpList, "Bond");
      
      MPOperator& Hamiltonian = OpList["H"];

      // Split operators.  H1 is the bilinear term, H2 is biquadratic
      MPOperator& H1 = OpList["H1"];
      MPOperator& H2 = OpList["H2"];
      MPOperator& H1p = OpList["H1p"];
      MPOperator& H2p = OpList["H2p"];
      MPOperator& H1px = OpList["H1px"];
      MPOperator& H2px = OpList["H2px"];
      MPOperator& H1Edge = OpList["H1Edge"];
      MPOperator& H2Edge = OpList["H2Edge"];
      MPOperator& H1Bulk = OpList["H1Bulk"];
      MPOperator& H2Bulk = OpList["H2Bulk"];

      QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
      // interaction matrix elements
      std::cout << "Working" << std::flush;
      for (int i = 1; i < L; ++i)
      {
         //double CurrentJ = (i == 1 || i == L) ? JEdge : J;
         //double CurrentBeta = (i == 1 || i == L) ? BetaEdge : Beta;

         MPOperator SS = -sqrt(3.0) * prod(S(i), S(i%L+1), Ident);
         MPOperator SS2 = prod(SS, SS, Ident);

         if (i == 1 || i == L-1)
         {
            H1Edge += SS;
            H2Edge += SS2;
         }
         else
         {
            H1Bulk += SS;
            H2Bulk += SS2;
         }
         std::cout << '.' << std::flush;
      }
      std::cout << "done!" << std::endl;
      H1p = -sqrt(3.0) * prod(S(1), S(L), Ident);
      H2p = prod(H1p, H1p, Ident);
      H1px = -sqrt(3.0) * prod(S(2), S(L-1), Ident);
      H2px = prod(H1px, H2px, Ident);
      H1 = JEdge*H1Edge + J*H1Bulk;
      H2 = BetaEdge*H2Edge + Beta*H2Bulk;

      Hamiltonian = H1+H2;

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
