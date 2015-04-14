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
      double JpEdge = 0;
      double JzEdge = 0;
      double Beta = 0;
      double BetaEdge = 0;
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
         ("J", prog_opt::value(&J), "isotropic bilinear coupling [default 1]")
         ("Jp", prog_opt::value(&Jp), "XY bilinear coupling [default J]")
         ("Jz", prog_opt::value(&Jz), "Z bilinear coupling [default J]")
         ("JpEdge", prog_opt::value(&JpEdge), "XY bilinear coupling at the edge bonds [default Jp]")
         ("JzEdge", prog_opt::value(&JzEdge), "Z bilinear coupling at the edge bonds [default Jz]")
         ("Beta", prog_opt::value(&Beta), "nearest neighbor biquadratic coupling [default 0]")
         ("BetaEdge", prog_opt::value(&BetaEdge), "biquadratic coupling at the edge bonds [default Beta]")
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
		   << "H1pBulk - XY bilinear spin exchange for the bulk\n"
		   << "H1zBulk - Z bilinear spin exchange for the bulk\n"
		   << "H1pEdge - XY bilinear spin exchange at the edges\n"
		   << "H1zEdge - Z bilinear spin exchange at the edges\n"
		   << "H2Bulk - biquadratic spin exchange for the bulk\n"
		   << "H2Edge - biquadratic spin exchange at the edges\n"
		   << "\nOr the shortcuts:\n"
		   << "H1p = JpEdge*H1pEdge + Jp*H1pBulk\n"
		   << "H1z = JzEdge*H1zEdge + Jz*H1zBulk\n"
		   << "H1 = H1p + H1z\n"
		   << "H2 = BetaEdge*H2Edge + Beta*H2Bulk\n"
		   << "H = H1 + H2\n";
         return 1;
      }

      // Set the dependent defaults if necessary
      if (!vm.count("SpinEdge"))  SpinEdge = Spin;
      if (!vm.count("Jp"))        Jp = J;
      if (!vm.count("Jz"))        Jz = J;
      if (!vm.count("JpEdge"))    JpEdge = Jp;
      if (!vm.count("JzEdge"))    JzEdge = Jz;
      if (!vm.count("BetaEdge"))  BetaEdge = Beta;

      TRACE(Spin)(SpinEdge)(Jp)(Jz)(JpEdge)(JzEdge)(Beta)(BetaEdge);
   
      // Construct the site block
      SiteBlock EdgeSite = CreateU1SpinSite(SpinEdge);
      SiteBlock BulkSite = CreateU1SpinSite(Spin);

      // construct a lattice of L copies of Site
      Lattice MyLattice = join(EdgeSite, repeat(BulkSite, L-2), EdgeSite);
      MyLattice.fix_coordinates();

      // construct the operator list for the lattice
      OperatorList OpList(MyLattice);
      
      OperatorAtSite<OperatorList const, int> Sp(OpList, "Sp");
      OperatorAtSite<OperatorList const, int> Sm(OpList, "Sm");
      OperatorAtSite<OperatorList const, int> Sz(OpList, "Sz");
      OperatorAtSite<OperatorList const, int> mSz(OpList, "mSz");
      OperatorAtSite<OperatorList, int> Bond(OpList, "Bond");
      MPOperator& Hamiltonian = OpList["H"];
      MPOperator& TotalSp = OpList["Sp"];
      MPOperator& TotalSm = OpList["Sm"];
      MPOperator& TotalSz = OpList["Sz"];
      MPOperator& TotalS2 = OpList["S2"];
      MPOperator& TotalmSz = OpList["mSz"]; // staggered magnetization
      MPOperator& M = OpList["M"];
      // Split operators.  H1 is the bilinear term, H2 is biquadratic
      MPOperator& H1pEdge = OpList["H1pEdge"];
      MPOperator& H1zEdge = OpList["H1zEdge"];
      MPOperator& H1pBulk = OpList["H1pBulk"];
      MPOperator& H1zBulk = OpList["H1zBulk"];
      MPOperator& H2Edge = OpList["H2Edge"];
      MPOperator& H2Bulk = OpList["H2Bulk"];
      MPOperator& H1p = OpList["H1p"];
      MPOperator& H1z = OpList["H1z"];
      MPOperator& H1 = OpList["H1"];
      MPOperator& H2 = OpList["H2"];

      QuantumNumber Ident(MyLattice.GetSymmetryList());  // the scalar quantum number
      // interaction matrix elements
      std::cout << "Working" << std::flush;
      for (int i = 1; i < L; ++i)
      {
	 MPOperator SSperp = 0.5 * (Sp(i)*Sm(i%L+1) + Sm(i)*Sp(i%L+1));
	 MPOperator SSz = Sz(i)*Sz(i%L+1);
	 MPOperator SS2 = (SSperp+SSz)*(SSperp+SSz);

         if (i == 1 || i == L-1) // edge interactions
         {
	    H1pEdge += SSperp;
	    H1zEdge += SSz;
	    H2Edge += SS2;
	 }
	 else // bulk
	 {
	    H1pBulk += SSperp;
	    H1zBulk += SSz;
	    H2Bulk += SS2;
	 }
         std::cout << '.' << std::flush;
      }
      std::cout << "done!" << std::endl;
 
      H1p = JpEdge*H1pEdge + Jp*H1pBulk;
      H1z = JzEdge*H1zEdge + Jz*H1zBulk;
      H1 = H1p + H1z;
      H2 = BetaEdge*H2Edge + Beta*H2Bulk;
      Hamiltonian = H1+H2;

      // The total spin operators
      for (int i = 1; i <= L; ++i)
      {
	 TotalSp += Sp(i);
	 TotalSm += Sm(i);
	 TotalSz += Sz(i);
	 TotalmSz += minus1pow(i) * Sz(i);
      }

      TotalS2 = 0.5 * (prod(TotalSp, TotalSm, Ident) + prod(TotalSm, TotalSp, Ident))
	 + prod(TotalSz, TotalSz, Ident);

      // The Marshall sign operator
      M = OpList["I"];
      for (int i = 1; i <= L; i += 2)
      {
         M = M * mSz(i);
      }

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
