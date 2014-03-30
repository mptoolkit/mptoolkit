
#include "mpo/triangular_mpo.h"
#include "mps/infinitewavefunction.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include "common/environment.h"
#include "interface/inittemp.h"
#include "mp-algorithms/gmres.h"
#include "mp-algorithms/arnoldi.h"
#include "common/polynomial.h"
#include "tensor/tensor_eigen.h"
#include "common/prog_opt_accum.h"
#include "mp/copyright.h"
#include "common/prog_options.h"

#include "models/spin.h"
#include "models/spin-u1.h"
#include "models/spin-su2.h"
#include "models/tj-u1su2.h"
#include "models/tj-u1.h"
#include "models/spinlessfermion-u1.h"
#include "models/kondo-u1su2.h"
#include "models/kondo-u1.h"
#include "models/hubbard-so4.h"
#include "models/boson-u1.h"
#include "models/hubbard-u1su2.h"
#include "models/hubbard-u1u1-old.h"
#include "models/hubbard-u1u1.h"

#include "mps/momentum_operations.h"
#include "mp-algorithms/triangular_mpo_solver.h"

namespace prog_opt = boost::program_options;


int main(int argc, char** argv)
{
   std::string FName;
   std::string Operator;

   double t = 0;
   double t2 = 0;
   double tc = 0;
   double U = 0;
   int Power = 1;
   bool Verbose = false;
   int NMax = 3;
   double Spin = 0.5;

   std::cout.precision(getenv_or_default("MP_PRECISION", 14));
 
   try
   {
      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("t", prog_opt::value(&t),
	  FormatDefault("nearest-neighbor hopping (for hubbard etc)", t).c_str())
	 ("t2", prog_opt::value(&t2),
	  FormatDefault("next-nearest-neighbor hopping (for hubbard etc)", t2).c_str())
	 ("tc", prog_opt::value(&tc),
	  FormatDefault("cluster hopping (for triangular cluster)", tc).c_str())
	 ("U", prog_opt::value(&U),
	  FormatDefault("coulomb repulsion", U).c_str())
	 ("spin,s", prog_opt::value(&Spin),
	  FormatDefault("spin", Spin).c_str())
	 ("nmax", prog_opt::value(&NMax),
          FormatDefault("Maximum number of particles (for bose-hubbard model)", NMax).c_str())
	 ("power", prog_opt::value(&Power),
	  FormatDefault("Calculate expectation value of operator to this power", Power).c_str())
         ("verbose,v", prog_opt::bool_switch(&Verbose),
          "extra debug output")
	 ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi", prog_opt::value<std::string>(&FName), "wavefunction")
         ("operator", prog_opt::value<std::string>(&Operator), "operator")
         ;

      prog_opt::positional_options_description p;
      p.add("psi", 1);
      p.add("operator", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);    

      if (vm.count("help") > 0 || vm.count("operator") == 0)
      {
         print_copyright(std::cerr);
         std::cerr << "usage: mp-expectation3 <psi1> <operator>\n";
         std::cerr << desc << '\n';
         return 1;
      }
      
      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pvalue_ptr<InfiniteWavefunction> PsiPtr = pheap::OpenPersistent(FName, CacheSize, true);
      InfiniteWavefunction Psi = *PsiPtr;
      int UnitCellSize = Psi.Psi.size();

      TriangularMPO Op;

      if (Operator == "su2spin")
      {
	 LatticeSite Site = CreateSU2SpinSite(Spin);	 
	 Op = TriangularOneSite(Site["S"]);
	 Op = -sqrt(3.0) * prod(Op, Op, QuantumNumber(Op.GetSymmetryList()));
      }
      else if (Operator == "bh-M")
      {
	 LatticeSite Site = BosonU1(NMax);
	 std::vector<BasisList> Sites(2, Site["I"].Basis());
	 Op = OnePointOperator(Sites, 0, Site["N"]) - OnePointOperator(Sites, 1, Site["N"]);
      }
      else if (Operator == "bh-N")
      {
	 LatticeSite Site = BosonU1(NMax);
	 Op = TriangularOneSite(Site["N"]);
      }
      else if (Operator == "tri-string" || Operator == "tri-string-new"
               || Operator == "tri-string2" || Operator == "tri-string2-new")
      {
	 LatticeSite Site = (Operator == "tri-string" || Operator == "tri-string2") ? CreateU1U1HubbardOldOrderingSite() : CreateU1U1HubbardSite();
	 std::vector<BasisList> Sites(3, Site["I"].Basis());
	 std::vector<SimpleOperator> String(3, Site["ES"]);
	 Op = OnePointStringOperator(Sites, String, 0, Site["Sz"])
	     + OnePointStringOperator(Sites, String, 1, Site["Sz"])
	     + OnePointStringOperator(Sites, String, 2, Site["Sz"]);
	 
	 std::vector<SimpleOperator> Pt0(3, Site["I"]);
	 std::vector<SimpleOperator> Pt1 = Pt0;
	 std::vector<SimpleOperator> Pt2 = Pt0;
	 Pt0[0] = Site["Sz"];
	 Pt1[1] = Site["Sz"];
	 Pt2[2] = Site["Sz"];

	 TriangularMPO OpX = OneCellStringOperator(Sites, String, Pt0)
	    + OneCellStringOperator(Sites, String, Pt1)
	    + OneCellStringOperator(Sites, String, Pt2);
	    
	 std::vector<SimpleOperator> Ptx0(3, Site["ES"]);
	 std::vector<SimpleOperator> Ptx1 = Ptx0;
	 std::vector<SimpleOperator> Ptx2 = Ptx0;
	 Ptx0[0] = Site["ES"] * Site["Sz"];
	 Ptx1[1] = Site["ES"] * Site["Sz"];
	 Ptx2[2] = Site["ES"] * Site["Sz"];

	 TriangularMPO OpA = OneCellStringOperator(Sites, String, Ptx0)
	    + OneCellStringOperator(Sites, String, Ptx1)
	    + OneCellStringOperator(Sites, String, Ptx2);

	 Op = OpA * adjoint(OpX);

         if (Operator == "tri-string2" || Operator == "tri-string2-new")
            Op = Op * adjoint(Op);
      }
      else if (Operator == "tri-su2")
      {
	 LatticeSite Site = CreateU1SU2HubbardSite();
	 double tSqrt2 = (-sqrt(2.0)) * t;  // the -sqrt(2) is an SU(2) factor
	 double tcSqrt2 = (-sqrt(2.0)) * tc;  // the -sqrt(2) is an SU(2) factor
	 // 3-site unit cell
	 std::vector<BasisList> Sites(3, Site["I"].Basis());
	 TriangularMPO Ham;
	 Ham += -tcSqrt2 * (TwoPointStringOperator(Sites, 0, Site["CHP"], Site["P"], 1, Site["C"])
			    + TwoPointStringOperator(Sites, 0, Site["CP"], Site["P"], 1, Site["CH"]));
	 Ham += -tcSqrt2 * (TwoPointStringOperator(Sites, 1, Site["CHP"], Site["P"], 2, Site["C"])
			    + TwoPointStringOperator(Sites, 1, Site["CP"], Site["P"], 2, Site["CH"]));
	 Ham += -tcSqrt2 * (TwoPointStringOperator(Sites, 0, Site["CHP"], Site["P"], 2, Site["C"])
			    + TwoPointStringOperator(Sites, 0, Site["CP"], Site["P"], 2, Site["CH"]));
	 Ham += -tSqrt2 * (TwoPointStringOperator(Sites, 1, Site["CHP"], Site["P"], 4, Site["C"])
			    + TwoPointStringOperator(Sites, 1, Site["CP"], Site["P"], 4, Site["CH"]));
         if (t2 != 0)
         {
	    Ham += - (-sqrt(2.0)*t2) * (TwoPointStringOperator(Sites, 1, Site["CHP"], Site["P"], 7, Site["C"])
			              + TwoPointStringOperator(Sites, 1, Site["CP"], Site["P"], 7, Site["CH"]));
         }
	 Ham += U * OnePointOperator(Sites, 0, Site["Pdouble"]);
	 Ham += U * OnePointOperator(Sites, 1, Site["Pdouble"]);
	 Ham += U * OnePointOperator(Sites, 2, Site["Pdouble"]);
	 Op = Ham;
      }
      else if (Operator == "tri-u1")
      {
	 LatticeSite Site = CreateU1U1HubbardSite();
	 // 3-site unit cell
	 std::vector<BasisList> Sites(3, Site["I"].Basis());
	 TriangularMPO Ham;
	 Ham += -tc * (TwoPointStringOperator(Sites, 0, Site["CHupP"], Site["P"], 1, Site["Cup"])
			    - TwoPointStringOperator(Sites, 0, Site["CupP"], Site["P"], 1, Site["CHup"])
		       + TwoPointStringOperator(Sites, 0, Site["CHdownP"], Site["P"], 1, Site["Cdown"])
			    - TwoPointStringOperator(Sites, 0, Site["CdownP"], Site["P"], 1, Site["CHdown"]));

	 Ham += -tc * (TwoPointStringOperator(Sites, 1, Site["CHupP"], Site["P"], 2, Site["Cup"])
			    - TwoPointStringOperator(Sites, 1, Site["CupP"], Site["P"], 2, Site["CHup"])
		       + TwoPointStringOperator(Sites, 1, Site["CHdownP"], Site["P"], 2, Site["Cdown"])
			    - TwoPointStringOperator(Sites, 1, Site["CdownP"], Site["P"], 2, Site["CHdown"]));

	 Ham += -tc * (TwoPointStringOperator(Sites, 0, Site["CHupP"], Site["P"], 2, Site["Cup"])
			    - TwoPointStringOperator(Sites, 0, Site["CupP"], Site["P"], 2, Site["CHup"])
		       + TwoPointStringOperator(Sites, 0, Site["CHdownP"], Site["P"], 2, Site["Cdown"])
			    - TwoPointStringOperator(Sites, 0, Site["CdownP"], Site["P"], 2, Site["CHdown"]));

	 Ham += -t * (TwoPointStringOperator(Sites, 1, Site["CHupP"], Site["P"], 4, Site["Cup"])
		      - TwoPointStringOperator(Sites, 1, Site["CupP"], Site["P"], 4, Site["CHup"])
		      + TwoPointStringOperator(Sites, 1, Site["CHdownP"], Site["P"], 4, Site["Cdown"])
		      - TwoPointStringOperator(Sites, 1, Site["CdownP"], Site["P"], 4, Site["CHdown"]));

	 if (t2 != 0)
	 {
	    Ham += -t2 * (TwoPointStringOperator(Sites, 1, Site["CHupP"], Site["P"], 7, Site["Cup"])
			  - TwoPointStringOperator(Sites, 1, Site["CupP"], Site["P"], 7, Site["CHup"])
			  + TwoPointStringOperator(Sites, 1, Site["CHdownP"], Site["P"], 7, Site["Cdown"])
			  - TwoPointStringOperator(Sites, 1, Site["CdownP"], Site["P"], 7, Site["CHdown"]));
	 }

	 Ham += U * OnePointOperator(Sites, 0, Site["Pdouble"]);
	 Ham += U * OnePointOperator(Sites, 1, Site["Pdouble"]);
	 Ham += U * OnePointOperator(Sites, 2, Site["Pdouble"]);
	 Op = Ham;
      }
      else
      {
	 std::cerr << "Unknown operator!\n";
	 exit(1);
      }

      // construct the operator to the given power
      TriangularMPO Temp = Op;
      while (Power > 1)
      {
	 Op = Op * Temp;
	 --Power;
      }

      // Make a LinearWavefunction in the symmetric orthogonality constraint
      MatrixOperator Identity = MatrixOperator::make_identity(Psi.Psi.Basis1());
      MatrixOperator Rho = scalar_prod(Psi.C_right, herm(Psi.C_right));
      LinearWavefunction Phi = Psi.Psi; // no need to bugger around with C_old,C_right
 
      Rho = delta_shift(Rho, Psi.QShift);

#if 0
      MatrixOperator LambdaSqrt = SqrtDiagonal(Psi.C_old);
      MatrixOperator LambdaInvSqrt = InvertDiagonal(LambdaSqrt, InverseTol);

      Phi.set_front(prod(LambdaInvSqrt, Phi.get_front()));
      Phi.set_back(prod(Phi.get_back(), delta_shift(LambdaSqrt, adjoint(Psi.QShift))));
      Rho = Psi.C_old;
      Identity = Rho;
#endif

      // Rho and Identity are the same matrix in this representation

      //   TRACE(norm_frob(operator_prod(herm(Phi), Identity, Phi) - Identity));
      //   TRACE(norm_frob(operator_prod(Phi, Rho, herm(Phi)) - Rho));

      // make Op the same size as our unit cell
      Op = repeat(Op, UnitCellSize / Op.size());

      KMatrixPolyType E = SolveMPO_Left(Phi, Psi.QShift, Op, Identity, Rho, Verbose);
      Polynomial<std::complex<double> > aNorm = ExtractOverlap(E[1.0], Rho);

      std::cout << "#degree #real #imag\n";
      for (int i = 0; i <= aNorm.degree(); ++i)
	 {
	    std::cout << i << ' ' << aNorm[i].real() << ' ' << aNorm[i].imag() << '\n';
	 }
      std::cout << std::endl;

      pheap::Shutdown();
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
