// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/mpoperator.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/wavefunc-utils.h"
#include "matrixproduct/triangularoperator.h"
#include "matrixproduct/infinitewavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "mp-algorithms/lanczos.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include "interface/operator-parser.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include <iostream>
#include "common/environment.h"

#include "models/spin-su2.h"
#include "models/spin-u1.h"
#include "models/spin.h"
#include "models/hubbard-u1su2.h"
#include "tensor/tensor_eigen.h"
#include "tensor/regularize.h"

#include "interface/inittemp.h"

namespace prog_opt = boost::program_options;

struct SuperblockMultiply
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator value_type;
   typedef MatrixOperator argument_type;

   SuperblockMultiply(MPStateComponent const& Left_,
		      MPStateComponent const& Right_);

   MatrixOperator operator()(MatrixOperator const& Psi) const
   {
      return operator_prod(Left, Psi, herm(Right));
   }

   MPStateComponent Left, Right;
};

inline
SuperblockMultiply::SuperblockMultiply(MPStateComponent const& Left_,
				       MPStateComponent const& Right_)
   : Left(Left_), Right(Right_)
{
}

bool ExpandL = true, ExpandR = true;

MatrixOperator
DoDMRGSweepLeft(LinearWavefunction& Psi, 
		MatrixOperator const& C_r, 
		SimpleMPOperator const& Ham,
		std::deque<MPStateComponent>& LeftBlockHam,
		MPStateComponent const& IncomingHam,
		StatesInfo const& SInfo)
{
   LinearWavefunction Result;
   std::deque<MPStateComponent> RightBlockHam;

   LinearWavefunction::const_iterator I = Psi.end();
   SimpleMPOperator::const_iterator H = Ham.end();
   --I; --H;

   MPStateComponent R = prod(*I, C_r);
   MatrixOperator C = ExpandBasis1(R);
   RightBlockHam.push_front(operator_prod(*H, R, IncomingHam, herm(R)));
   LeftBlockHam.pop_back();

   while (I != Psi.begin())
   {
      // Expand the left matrix
      --I;
      --H;
      MPStateComponent L = prod(*I, C);
      C = ExpandBasis2(L);

      LeftBlockHam.pop_back();
      LeftBlockHam.push_back(operator_prod(herm(*H), herm(L), LeftBlockHam.back(), L));

      // apply the solver
      int Iterations = 10;
      double Energy = Lanczos(C, SuperblockMultiply(LeftBlockHam.back(), RightBlockHam.front()),
			      Iterations);

      // truncate
      MatrixOperator Rho = scalar_prod(herm(C), C);
      DensityMatrix<MatrixOperator> DM(Rho);
      TruncationInfo Info;
      DensityMatrix<MatrixOperator>::const_iterator DMPivot =
         TruncateFixTruncationErrorAbsolute(DM.begin(), DM.end(),
                                            SInfo,
                                            Info);
      MatrixOperator U = DM.ConstructTruncator(DM.begin(), DMPivot);

      std::cout << "Energy=" << Energy 
		<< " States=" << Info.KeptStates()
		<< " TruncError=" << Info.TruncationError()
		<< " Entropy=" << Info.KeptEntropy() << '\n';

      C = C * herm(U);
      R = prod(U, R);
      RightBlockHam.front() = triple_prod(U, RightBlockHam.front(), herm(U));

      // shift left
      Result.push_front(R);
      R = prod(L, C);
      C = ExpandBasis1(R);
      RightBlockHam.push_front(operator_prod(*H, R, RightBlockHam.front(), herm(R)));
      LeftBlockHam.pop_back();
   }

   // cleanup
   Result.push_front(R);
   LeftBlockHam = RightBlockHam;
   Psi = Result;
   return C;
}

MatrixOperator
DoDMRGSweepRight(MatrixOperator const& C_l, 
		 LinearWavefunction& Psi, 
		 SimpleMPOperator const& Ham,
		 MPStateComponent const& IncomingHam,
		 std::deque<MPStateComponent>& RightBlockHam,
		 StatesInfo const& SInfo)
{
   LinearWavefunction Result;
   std::deque<MPStateComponent> LeftBlockHam;

   LinearWavefunction::const_iterator I = Psi.begin();
   SimpleMPOperator::const_iterator H = Ham.begin();

   MPStateComponent L = prod(C_l, *I);
   MatrixOperator C = ExpandBasis2(L);
   LeftBlockHam.push_back(operator_prod(herm(*H), herm(L), IncomingHam, L));
   RightBlockHam.pop_front();

   ++I; ++H;

   while (I != Psi.end())
   {
      MPStateComponent R = prod(C, *I);
      C = ExpandBasis1(R);

      RightBlockHam.pop_front();
      RightBlockHam.push_front(operator_prod(*H, R, RightBlockHam.back(), herm(R)));

      // apply the solver
      int Iterations = 10;
      double Energy = Lanczos(C, SuperblockMultiply(LeftBlockHam.back(), RightBlockHam.front()),
			      Iterations);

      // truncate
      MatrixOperator Rho = scalar_prod(C, herm(C));
      DensityMatrix<MatrixOperator> DM(Rho);
      TruncationInfo Info;
      DensityMatrix<MatrixOperator>::const_iterator DMPivot =
         TruncateFixTruncationErrorAbsolute(DM.begin(), DM.end(),
                                            SInfo,
                                            Info);
      MatrixOperator U = DM.ConstructTruncator(DM.begin(), DMPivot);

      std::cout << "Energy=" << Energy 
		<< " States=" << Info.KeptStates()
		<< " TruncError=" << Info.TruncationError()
		<< " Entropy=" << Info.KeptEntropy() << '\n';

      C = U * C;
      L = prod(L, herm(U));
      LeftBlockHam.front() = triple_prod(U, LeftBlockHam.front(), herm(U));

      // shift right
      Result.push_back(L);
      L = prod(C, R);
      C = ExpandBasis2(L);
      LeftBlockHam.push_back(operator_prod(herm(*H), herm(L), LeftBlockHam.back(), L));
      RightBlockHam.pop_front();

      ++I;
      ++H;
   }

   // cleanup
   Result.push_back(L);
   RightBlockHam = LeftBlockHam;
   Psi = Result;
   return C;
}

double const InvertEpsilon = 1E-8;

int niter = 0;

void DoIteration(LinearWavefunction& Left, MatrixOperator& C_LR, LinearWavefunction& Right,
		 SimpleMPOperator const& LeftHam, SimpleMPOperator const& RightHam,
		 MatrixOperator& C_RL, 
		 std::deque<MPStateComponent>& LeftBlockHam,
		 std::deque<MPStateComponent>& RightBlockHam,
		 StatesInfo SInfo, int Iterations)
{
   //TRACE(C_LR);
   int Iter = Iterations;

   // These two could be done in parallel
   MPStateComponent LBack = LeftBlockHam.back();
   MatrixOperator C_left = DoDMRGSweepLeft(Left, C_LR, LeftHam, 
					   LeftBlockHam, RightBlockHam.front(), SInfo);

   MatrixOperator C_right = DoDMRGSweepRight(C_LR, Right, RightHam, 
					     LBack, RightBlockHam, SInfo);
   
   //TRACE(C_RL)(InvertDiagonal(C_RL));

   // update C_RL
   //TRACE(C_RL)(InvertDiagonal(C_RL, 1E-7));
   C_RL = C_right * InvertDiagonal(C_RL, InvertEpsilon) * C_left;
   //TRACE(C_RL)(C_right)(C_left);
   
   //TRACE(C_left)(C_LR)(C_RL)(C_right)(C_LR.Basis1())(C_LR.Basis2())(C_left.Basis2())(C_right.Basis1());

   // solve

   double Energy;
   double PsiDifference;
   {
      Iterations = Iter;
      MatrixOperator C_RL_Old = C_RL;
      C_RL_Old *= 1.0 / norm_frob(C_RL_Old);
      C_RL *= 1.0 / norm_frob(C_RL);
      Energy = Lanczos(C_RL, SuperblockMultiply(RightBlockHam.back(), LeftBlockHam.front()),
		       Iterations);
      PsiDifference = 1.0 - norm_frob(inner_prod(C_RL, C_RL_Old));
      // adjust for the energy per site
      RightBlockHam.back().front() -= (0.5 * Energy) * MatrixOperator::make_identity(C_RL.Basis1());
      LeftBlockHam.front().back() -= (0.5 * Energy) * MatrixOperator::make_identity(C_RL.Basis2());
   }

   // truncate
   //SInfo.MaxStates = 40;
   {
      MatrixOperator RhoL = scalar_prod(C_RL, herm(C_RL));
      DensityMatrix<MatrixOperator> DML(RhoL);
      TruncationInfo Info;
      MatrixOperator TruncL = DML.ConstructTruncator(DML.begin(), 
						     TruncateFixTruncationErrorAbsolute(DML.begin(),
											DML.end(),
											SInfo,
											Info));
      std::cout << "A Energy=" << Energy 
		<< " States=" << Info.KeptStates()
		<< " TruncError=" << Info.TruncationError()
		<< " Entropy=" << Info.KeptEntropy() 
                << " PsiDiff=" << PsiDifference
                << '\n';
      //DML.DensityMatrixReport(std::cout);

      C_RL = TruncL * C_RL;
      Right.set_back(prod(Right.get_back(), herm(TruncL)));
      RightBlockHam.back() = triple_prod(TruncL, RightBlockHam.back(), herm(TruncL));

      MatrixOperator RhoR = scalar_prod(herm(C_RL), C_RL);
      DensityMatrix<MatrixOperator> DMR(RhoR);
      TruncationInfo InfoR;
      MatrixOperator TruncR = DMR.ConstructTruncator(DMR.begin(), 
						     TruncateFixTruncationErrorAbsolute(DMR.begin(),
											DMR.end(),
											SInfo,
											InfoR));
      C_RL = C_RL * herm(TruncR);
      Left.set_front(prod(TruncR, Left.get_front()));
      LeftBlockHam.front() = triple_prod(TruncR, LeftBlockHam.front(), herm(TruncR));
   }

   //TRACE(C_RL);

   // sweep back
   MPStateComponent RBack = RightBlockHam.back();
   C_left = DoDMRGSweepLeft(Right, C_RL, RightHam, RightBlockHam, LeftBlockHam.front(), SInfo);
   C_right = DoDMRGSweepRight(C_RL, Left, LeftHam, RBack, LeftBlockHam, SInfo);

   //   TRACE(C_left)(C_LR)(C_right)(C_LR.Basis1())(C_LR.Basis2())(C_right.Basis2())(C_left.Basis1());

   //TRACE(C_LR)(InvertDiagonal(C_LR, 1E-7));
   C_LR = C_right * InvertDiagonal(C_LR, InvertEpsilon) * C_left;
   //TRACE(C_LR)(C_left)(C_right);
   // update C_LR

   // solve
   {
      Iterations = Iter;
      MatrixOperator C_LR_Old = C_LR;
      C_LR_Old *= 1.0 / norm_frob(C_LR_Old);
      C_LR *= 1.0 / norm_frob(C_LR);
      Energy = Lanczos(C_LR, SuperblockMultiply(LeftBlockHam.back(), RightBlockHam.front()),
		       Iterations);
      PsiDifference = 1.0 - norm_frob(inner_prod(C_LR, C_LR_Old));
      // adjust for the energy per site
      LeftBlockHam.back().front() -= (0.5 * Energy) * MatrixOperator::make_identity(C_LR.Basis1());
      RightBlockHam.front().back() -= (0.5 * Energy) * MatrixOperator::make_identity(C_LR.Basis2());
   }

   // truncate
   //SInfo.MaxStates = 40;
   {
      MatrixOperator RhoL = scalar_prod(C_LR, herm(C_LR));
      DensityMatrix<MatrixOperator> DML(RhoL);
      TruncationInfo Info;
      MatrixOperator TruncL = DML.ConstructTruncator(DML.begin(), 
						     TruncateFixTruncationErrorAbsolute(DML.begin(),
											DML.end(),
											SInfo,
											Info));
      std::cout << "B Energy=" << Energy 
		<< " States=" << Info.KeptStates()
		<< " TruncError=" << Info.TruncationError()
		<< " Entropy=" << Info.KeptEntropy() 
                << " PsiDiff=" << PsiDifference
                << '\n';
      //DML.DensityMatrixReport(std::cout);

      C_LR = TruncL * C_LR;
      Left.set_back(prod(Left.get_back(), herm(TruncL)));
      LeftBlockHam.back() = triple_prod(TruncL, LeftBlockHam.back(), herm(TruncL));

      MatrixOperator RhoR = scalar_prod(herm(C_LR), C_LR);
      DensityMatrix<MatrixOperator> DMR(RhoR);
      TruncationInfo InfoR;
      MatrixOperator TruncR = DMR.ConstructTruncator(DMR.begin(), 
						     TruncateFixTruncationErrorAbsolute(DMR.begin(),
											DMR.end(),
											SInfo,
											InfoR));
      C_LR = C_LR * herm(TruncR);
      Right.set_front(prod(TruncR, Right.get_front()));
      RightBlockHam.front() = triple_prod(TruncR, RightBlockHam.front(), herm(TruncR));
   }

}

int main(int argc, char** argv)
{
   try
   {
      int NumIter = 4;
      int MinStates = 1;
      int MaxStates = 100000;
      double MixFactor = 0.01;
      bool TwoSite = false;
      int NumSteps = 10;
      double TruncCutoff = 0;
      double EigenCutoff = -1;
      double Lambda = 1.0;
      double J2 = 0.0;
      double Theta = 0.0;
      half_int Spin = 0.5;
      bool StartFromFixedPoint = false;
      std::string OutName = "test.out";

      std::cout.precision(14);

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("iter,i", prog_opt::value<int>(&NumIter), "Number of Lanczos iterations per step [default 4]")
	 ("max-states,m", prog_opt::value<int>(&MaxStates), 
          "Maximum number of states to keep [default 100000]")
         ("min-states", prog_opt::value<int>(&MinStates), "Minimum number of states to keep [default 1]")
         ("trunc,r", prog_opt::value<double>(&TruncCutoff), 
          "Truncation error cutoff [default 0]")
         ("eigen-cutoff,d", prog_opt::value(&EigenCutoff), 
          ("Cutoff threshold for density matrix eigenvalues (alternative to truncation error) [default "
           +boost::lexical_cast<std::string>(EigenCutoff)+"]").c_str())
	 ("steps,s", prog_opt::value<int>(&NumSteps), "Number of DMRG steps to perform [default 10]")
	 ("out,o", prog_opt::value(&OutName), "Output filename")
	 ("spin", prog_opt::value(&Spin), "spin (for xxx,xxz,xyz hamiltonians) [default 0.5]")
	 ("J2", prog_opt::value(&J2), "next-nearest-neighbor hopping J2 (for xxx) [default 0]")
	 ("theta", prog_opt::value(&Theta), "biquadratic term (biquadratic model, in units of pi) [default 0]")
	 ("lambda", prog_opt::value(&Lambda), "transverse field strength"
	  " (for itf hamiltonian) [default 1.0])")
	 ;

      prog_opt::options_description opt;
      opt.add(desc);

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).run(), vm);
      prog_opt::notify(vm);    

      if (vm.count("help") || vm.count("out") == 0) 
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-idmrg [options]\n";
         std::cerr << desc << '\n';
         return 1;
      }


      int PageSize = getenv_or_default("MP_PAGESIZE", 65536);
      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pheap::Initialize(OutName, 1, PageSize, CacheSize);

      StatesInfo SInfo;
      SInfo.MinStates = MinStates;
      SInfo.MaxStates = MaxStates;
      SInfo.TruncationCutoff = TruncCutoff;
      SInfo.EigenvalueCutoff = EigenCutoff;
      std::cout << SInfo << '\n';

   // Ising model
#if 0
      TRACE(Lambda);
   SiteBlock Boundary = CreateSpinSite(0.5);
   SiteBlock Site = CreateSpinSite(0.5);

   MpOpTriangular Ham = 4.0 * TriangularTwoSite(Site["Sz"], Site["Sz"])
      + Lambda * 2.0 * TriangularOneSite(Site["Sx"]);

   MpOpTriangular BoundaryHam = Ham;
#endif

   // SU(2) Heisenberg model
#if 0
   SiteBlock Boundary = CreateSU2SpinSite(0.5);
   SiteBlock Site = CreateSU2SpinSite(0.5);

   MpOpTriangular Ham = TriangularTwoSite(-sqrt(3.0)*Site["S"], 
   					  Site["S"], 
   					  Site["I"].TransformsAs());

   MpOpTriangular BoundaryHam = TriangularTwoSite(-sqrt(3.0)*Boundary["S"], 
   						  Boundary["S"], 
   						  Boundary["I"].TransformsAs());
#endif

   // SU(2) bilinear-biquadratic Heisenberg model
#if 1
   SiteBlock Boundary = CreateSU2SpinSite(0.5);
   SiteBlock Site = CreateSU2SpinSite(0.5);


   double J = cos(Theta * math_const::pi);
   double B = 1.2*sin(Theta * math_const::pi);
   //   double J=1, B=0.4;

   TRACE(Theta)(J)(B);

   MpOpTriangular Ham = J*TriangularTwoSite(-sqrt(3.0)*Site["S"], 
   					  Site["S"], 
   					  Site["I"].TransformsAs())
      + B*TriangularTwoSite(sqrt(5)*Site["Q"],
                          Site["Q"],
                          Site["I"].TransformsAs());

   MpOpTriangular BoundaryHam = J*TriangularTwoSite(-sqrt(3.0)*Boundary["S"], 
   						  Boundary["S"], 
   						  Boundary["I"].TransformsAs())
      + B*TriangularTwoSite(sqrt(5)*Boundary["Q"],
                          Boundary["Q"],
                          Boundary["I"].TransformsAs());


#endif

   // SU(2) Heisenberg model with PBC
#if 0
   SiteBlock Boundary = CreateSU2SpinSite(0.5);
   SiteBlock Site = CreateSU2SpinSite(0.5);

   MpOpTriangular Ham = TriangularTwoSitePBC(-sqrt(3.0)*Site["S"], 
                                             Site["S"], 
                                             Site["I"].TransformsAs());

   MpOpTriangular BoundaryHam = TriangularTwoSitePBC_Boundary(-sqrt(3.0)*Boundary["S"], 
                                                              Boundary["S"], 
                                                              Boundary["I"].TransformsAs());
#endif

   // Heisenberg model with no symmetries
#if 0
   SiteBlock Boundary = CreateSpinSite(0.5);
   SiteBlock Site = CreateSpinSite(0.5);

   MpOpTriangular Ham = 0.5 * (TriangularTwoSite(Site["Sp"], Site["Sm"]) + 
			       TriangularTwoSite(Site["Sm"], Site["Sp"]))
      + TriangularTwoSite(Site["Sz"], Site["Sz"]);

   MpOpTriangular BoundaryHam = 0.5 * (TriangularTwoSite(Boundary["Sp"], Boundary["Sm"]) + 
				       TriangularTwoSite(Boundary["Sm"], Boundary["Sp"]))
      + TriangularTwoSite(Boundary["Sz"], Boundary["Sz"]);

   TRACE(Ham.data());
#endif

   // SU(2) zig-zag Heisenberg
#if 0
   SiteBlock Boundary = CreateSU2SpinSite(0.5);
   SiteBlock Site = CreateSU2SpinSite(0.5);

   MpOpTriangular Ham = ZigZagChain(-sqrt(3.0)*Site["S"], Site["S"], J1, J2);
   MpOpTriangular BoundaryHam = ZigZagChain(-sqrt(3.0)*Boundary["S"], Boundary["S"], J1, J2);
#endif

   // SU(2) Kagome lattice
#if 0 
   SiteBlock Boundary = CreateSU2SpinSite(0.5);
   SiteBlock Site = CreateSU2SpinSite(0.5);

   MpOpTriangular Ham = Jcross * ();
#endif

   // U(1)xSU(2) Hubbard model
#if 0
   SiteBlock Boundary = CreateSU2HubbardSite();
   SiteBlock Site = CreateSU2HubbardSite();

   MpOpTriangular Ham = TriangularTwoSite(sqrt(2.0)*Site["CHP"], Site["C"], Site["I"].TransformsAs())
      + TriangularTwoSite(sqrt(2.0)*Site["CP"], Site["CH"], Site["I"].TransformsAs())
      + (x/4.0) * TriangularOneSite(Site["P"]);
   MpOpTriangular BoundaryHam = Ham;
#endif

   MPStateComponent E = Initial_E(BoundaryHam);
   MPStateComponent F = Initial_F(BoundaryHam);

   // initial 2-site block
   MatrixOperator Center = MatrixOperator::make_identity(E.Basis1());
   MPStateComponent A1 = ConstructFromLeftBasis(Boundary.Basis2().Basis(), Center.Basis1());
   MPStateComponent B2 = ConstructFromRightBasis(Boundary.Basis1().Basis(), Center.Basis2());
   E = operator_prod(herm(BoundaryHam.data()), herm(A1), E, A1);
   F = operator_prod(BoundaryHam.data(), B2, F, herm(B2));

   Center = MakeRandomMatrixOperator(A1.Basis2(), B2.Basis1());
   int Iterations = NumIter;
   double Energy = Lanczos(Center, 
			   SuperblockMultiply(E, F),
			   Iterations);
   TRACE(Energy);
   // switch to basis where Center is diagonal
   MatrixOperator U, D, Vt;
   SingularValueDecomposition(Center, U, D, Vt);
   Center = D;
   A1 = prod(A1, U);
   B2 = prod(Vt, B2);
   E = triple_prod(herm(U), E, U);
   F = triple_prod(Vt, F, herm(Vt));

   //   double BondEnergy = inner_prod(Center, triple_prod(E[1], Center, herm(F[1]))).real();

   // Add one unit cell away from the boundary
   A1 = ConstructFromLeftBasis(Site.Basis2().Basis(), A1.Basis2());
   B2 = ConstructFromRightBasis(Site.Basis1().Basis(), B2.Basis1());
   E = operator_prod(herm(Ham.data()), herm(A1), E, A1);
   F = operator_prod(Ham.data(), B2, F, herm(B2));

   MatrixOperator OldCenter = Center;
   Center = MakeRandomMatrixOperator(A1.Basis2(), B2.Basis1());

   Iterations = NumIter;
   Energy = Lanczos(Center, 
		    SuperblockMultiply(E, F),
		    Iterations);
   double WavefunctionDifference = 0;
   TRACE(Energy);

   // diagonalize the center matrix again
   SingularValueDecomposition(Center, U, D, Vt);
   Center = D;
   A1 = prod(A1, U);
   B2 = prod(Vt, B2);
   E = triple_prod(herm(U), E, U);
   F = triple_prod(Vt, F, herm(Vt));

   // set up

   LinearWavefunction Left, Right;
   Left.push_back(A1);
   Right.push_front(B2);

   SimpleMPOperator LeftHam, RightHam;
   LeftHam.push_back(Ham.data());
   RightHam.push_front(Ham.data());

   MatrixOperator C_LR = Center;
   MatrixOperator C_RL = OldCenter;

   std::deque<MPStateComponent> LeftBlockHam, RightBlockHam;
   LeftBlockHam.push_back(E);
   RightBlockHam.push_front(F);

   for (int i = 0; i < NumSteps; ++i)
   {
      DoIteration(Left, C_LR, Right, LeftHam, RightHam, C_RL, LeftBlockHam, RightBlockHam, SInfo, NumIter);
   }

      InfiniteWavefunction Psi;
      {
      // convert back to an InfiniteWavefunction
      Psi.C_old = C_RL;
      Psi.Psi = LinearWavefunction();
      for (LinearWavefunction::const_iterator I = Left.begin(); I != Left.end(); ++I)
      {
	 Psi.Psi.push_back(*I);
      }

      MatrixOperator U = C_LR;
      for (LinearWavefunction::const_iterator I = Right.begin(); I != Right.end(); ++I)
      {
	 MPStateComponent x = prod(U, *I);
	 U = TruncateBasis2(x);
	 Psi.Psi.push_back(x);
      }
      Psi.C_right = U;
      Psi.QShift = QuantumNumbers::QuantumNumber(U.GetSymmetryList());
      }

      pvalue_ptr<InfiniteWavefunction> PsiPtr = new InfiniteWavefunction(Psi);
      pheap::ShutdownPersistent(PsiPtr);

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
