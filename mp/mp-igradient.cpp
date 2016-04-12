// -*- C++ -*- $Id$

#include "mpo/triangularoperator.h"
#include "mps/infinitewavefunction.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include "common/environment.h"
#include "interface/inittemp.h"
#include "mp-algorithms/gmres.h"
#include "mp-algorithms/arnoldi.h"
#include "mp-algorithms/lanczos.h"
#include "mp-algorithms/lanczos-nonortho.h"
#include "tensor/tensor_eigen.h"
#include "mpo/mpo_functors.h"
#include "mp/copyright.h"

#include "models/spin-su2.h"
#include "models/spin-u1.h"
#include "models/spin-z2.h"
#include "models/spin.h"
#include "models/tj-u1su2.h"
#include "models/tj-u1.h"
#include "models/kondo-u1su2.h"
#include "models/kondo-u1.h"
#include "models/hubbard-so4.h"
#include "models/hubbard-u1u1.h"
#include "models/spinlessfermion-u1.h"
#include "models/bosehubbard-spinless.h"
#include "models/bosehubbard-spinless-u1.h"
#include "tensor/tensor_eigen.h"
#include "tensor/regularize.h"

#include "interface/inittemp.h"
#include "mp-algorithms/random_wavefunc.h"

namespace prog_opt = boost::program_options;

struct OneMinusTransfer
{
   OneMinusTransfer(double ScaleFactor, StateComponent const& Psi, MatrixOperator const& Rho,
                    MatrixOperator const& Identity)
      : Scale_(ScaleFactor), Psi_(Psi), Rho_(Rho), Identity_(Identity) {}

   MatrixOperator operator()(MatrixOperator const& x) const
   {
      MatrixOperator r = x-Scale_*operator_prod(herm(Psi_), x, Psi_);
      r -= inner_prod(r, Rho_) * Identity_; // orthogonalize to the identity
      return r;
   }

   double Scale_;
   StateComponent const& Psi_;
   MatrixOperator const& Rho_;
   MatrixOperator const& Identity_;
};

template <typename Func>
MatrixOperator
LinearSolve(Func F, MatrixOperator Rhs)
{
   MatrixOperator Guess = Rhs;
   int m = 30;
   int max_iter = 10000;
   double tol = 1e-15;
   GmRes(Guess, F, Rhs, m, max_iter, tol, LinearAlgebra::Identity<MatrixOperator>());
   TRACE(tol);
   return Guess;
}

double const SolverTol = 1E-12;

// Solve a first order MPO, for the symmetric case
// The Psi is passed by reference so that we can normalize it to 1
double
SolveFirstOrderMPO_Symmetric(StateComponent& Psi,
                             MpOpTriangular const& Op, 
                             StateComponent& Eigenvector,
                             int Verbose = 0)
{
   int Dim = Op.Basis1().size();  // dimension of the MPO

   // Solve recursively, starting from the first row, which will be the identity operator
   int Row = 0;


   // symmetrize Psi
   Psi[0] = 0.5 * (Psi[0] + adjoint(Psi[0]));
   Psi[1] = 0.5 * (Psi[1] - adjoint(Psi[1]));

   TRACE(norm_frob(Psi[0] - adjoint(Psi[0])));
   TRACE(norm_frob(Psi[1] + adjoint(Psi[1])));

   TRACE(norm_frob(Psi));

   int Iterations = 30;
   double Tol = SolverTol;
   // use Arnoldi, since we can't yet specify the SolverMode for Lanczos
   std::complex<double> Val = LinearSolvers::Arnoldi(Eigenvector[Row], 
                                                     OperatorProd_MAxBH(Op(Row, Row).scalar(), Psi, Psi),
                                                     Iterations,
                                                     Tol,
                                                     LinearSolvers::LargestAlgebraicReal,
                                                     false,
                                                     Verbose >= 1);
   Eigenvector[0] = 0.5*(Eigenvector[0] + adjoint(Eigenvector[0]));
   while (Iterations == 30)
   {
      TRACE("Restarting Arnoldi")(-Tol);
      Tol = SolverTol;
      Iterations = 30;
      Val =  LinearSolvers::Arnoldi(Eigenvector[Row], 
                                    OperatorProd_MAxBH(Op(Row, Row).scalar(), Psi, Psi),
                                    Iterations,
                                    Tol,
                                    LinearSolvers::LargestAlgebraicReal,
                                    false,
                                    Verbose >= 1);
      Eigenvector[0] = 0.5*(Eigenvector[0] + adjoint(Eigenvector[0]));
   }

   TRACE("Normalization")(Val);
   Psi *= 1.0 / std::sqrt(Val.real()); // normalize the MPS

   // rotate to a basis that diagonalizes Eigenvector[0]
   // if E = herm(U) D U
   // then D is the eigenvector of U Psi herm(U) 

   MatrixOperator D = Eigenvector[0];
   MatrixOperator U = DiagonalizeHermitian(D);

   Psi = triple_prod(U, Psi, herm(U));
   Eigenvector[0] = D;

   // rotate the guess H operator into the new basis
   Eigenvector[Dim-1] = triple_prod(U, Eigenvector[Dim-1], herm(U));

   // Verify that we have a left and right eigenvector
   double N1 = norm_frob(Eigenvector[Row] - OperatorProd_MAxBH(Op(Row, Row).scalar(), Psi, Psi)(Eigenvector[Row]));
   double N2 = norm_frob(Eigenvector[Row] - OperatorProd_MAHxB(Op(Row, Row).scalar(), Psi, Psi)(Eigenvector[Row]));
   CHECK(N1 < SolverTol*10)(N1);
   CHECK(N2 < SolverTol*10)(N2);

   // Make the trace real and positive
   Eigenvector[Row] *= norm_frob(trace(Eigenvector[Row])) / trace(Eigenvector[Row]);

   //   TRACE(SingularValues(Eigenvector[Row]));

   while (Row < Dim-1)
   {
      ++Row;
      Eigenvector[Row] = MatrixOperator();   // clear the next element
      // Assume the diagonal is empty
      for (int j = 0; j < Row; ++j)
      {
         SimpleRedOperator const M = Op(Row, j);
         if (!M.is_null())
         {
            Eigenvector[Row] += operator_prod(M, Psi, Eigenvector[j], herm(Psi), Op.Basis1()[Row]);
         }
      }
   }

   // Solve (1-T)x = Eigenvector[Row]
   MatrixOperator C = Eigenvector[Row];
   // decompose C into parts parallel and perpendicular to the identity
   std::complex<double> Energy = inner_prod(Eigenvector[0], C);
   C -= Energy*Eigenvector[0];

   // a DGKS-style correction, that gets a better orthogonalization of C against the identity part
   std::complex<double> E2 = inner_prod(Eigenvector[0], C);
   C -= E2*Eigenvector[0];
   Energy += E2;

   DEBUG_CHECK(norm_frob(inner_prod(Eigenvector[0], C)) < 1E-14);

   if (Verbose >= 1)
   {
      std::cerr << "Energy=" << Energy << '\n';
   }


   Eigenvector[Row] = LinearSolve(OneMinusTransfer(1.0, Psi, Eigenvector[0], Eigenvector[0]), C);
   TRACE(inner_prod(Eigenvector[0], Eigenvector[Row]));

   return Energy.real();
}

int main(int argc, char** argv)
{
   ProcControl::Initialize(argv[0], 0, 0, true);
   try
   {
      int NumIter = 20;
      int MinStates = 1;
      int MaxStates = 100000;
      //double MixFactor = 0.01;
      //bool TwoSite = false;
      int NumSteps = 10;
      double TruncCutoff = 0;
      double EigenCutoff = -1;
      //bool NoVariance = false;
      //bool UseDGKS = false;
      std::string FName;
      std::string HamStr;
      double Lambda = 1.0;
      double J = 1.0;
      double J2 = 0.0;
      double B = 0;
      double D = 0;
      double U = 0;
      double Jz = 0.0;
      double tprime = 1.0;
      double delta = 0.0;
      double Theta = 0.0;
      double Beta = 0.0;
      half_int Spin = 0.5;
      int NLegs = 1;
      bool TwoSite = true;
      bool OneSite = false;
      double MixFactor = 0.0;
      int Verbose = 10;
      bool NoOrthogonalize = false;

      std::cout.precision(14);

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Hamiltonian,H", prog_opt::value(&HamStr),
          "model Hamiltonian.  Valid choices: itf, itf-z2, xxx-su2, xxx-u1, xxx, tj-zigzag-u1su2, "
          "tj-zigzag-u1, sf-zigzag-u1, klm-u1su2, klm-u1, bh, bh2, bh-u1, bh2-u1")
         ("wavefunction,w", prog_opt::value(&FName),
          "wavefunction to apply DMRG (required)")
	 ("max-states,m", prog_opt::value<int>(&MaxStates),
	  FormatDefault("Maximum number of states to keep", MaxStates).c_str())
         ("min-states", prog_opt::value<int>(&MinStates),
	  FormatDefault("Minimum number of states to keep", MinStates).c_str())
         ("trunc,r", prog_opt::value<double>(&TruncCutoff),
          FormatDefault("Truncation error cutoff", TruncCutoff).c_str())
         ("eigen-cutoff,d", prog_opt::value(&EigenCutoff),
          FormatDefault("Cutoff threshold for density matrix eigenvalues", EigenCutoff).c_str())
	 ("mix-factor,f", prog_opt::value(&MixFactor),
	  FormatDefault("Mixing coefficient for the density matrix", MixFactor).c_str())
	 ("steps,s", prog_opt::value<int>(&NumSteps),
	  FormatDefault("Number of DMRG steps to perform", NumSteps).c_str())
	 ("no-orthogonalize", prog_opt::bool_switch(&NoOrthogonalize),
	  "Don't orthogonalize the wavefunction before saving")
	 ("maxiter", prog_opt::value<int>(&NumIter),
	  FormatDefault("Maximum number of Lanczos iterations per step (Krylov subspace size)", NumIter).c_str())
	 ("spin", prog_opt::value(&Spin),
	  FormatDefault("spin (for xxx,xxz,xyz hamiltonians)", Spin).c_str())
	 ("J", prog_opt::value(&J),
	  FormatDefault("nearest-neighbor exchange J (for xxx,itf, etc)", J).c_str())
	 ("J2", prog_opt::value(&J2),
	  FormatDefault("next-nearest-neighbor exchange J2 (for xxx)", J2).c_str())
	 ("D", prog_opt::value(&D),
	  FormatDefault("single-ion anisotropy (for xxx-u1 and xxx)", D).c_str())
	 ("U", prog_opt::value(&U),
	  FormatDefault("coulomb repulsion", U).c_str())
	 ("B", prog_opt::value(&B),
	  FormatDefault("magnetic field (for xxx)", B).c_str())
	 ("Jz", prog_opt::value(&Jz),
	  FormatDefault("Jz coupling (for Kondo)", Jz).c_str())
	 ("nlegs", prog_opt::value(&NLegs),
	  FormatDefault("Number of legs (for triangular ladder)", NLegs).c_str())
	 ("tprime", prog_opt::value(&tprime),
	  FormatDefault("next-nearest-neighbor hopping t' (for tj-zigzag, sf-zigzag)", tprime).c_str())
	 ("delta", prog_opt::value(&delta),
	  FormatDefault("Zigzag ladder potential imbalance (for tj-zigzag, sf-zigzag)", delta).c_str())
	 ("theta", prog_opt::value(&Theta),
	  FormatDefault("theta (for biquadratic xxx)", Theta).c_str())
	 ("Beta", prog_opt::value(&Beta),
	  FormatDefault("Beta (for biquadratic xxx)", Beta).c_str())
	 ("lambda", prog_opt::value(&Lambda),
	  FormatDefault("transverse field strength (for itf hamiltonian)", Lambda).c_str())
	  ;

      prog_opt::options_description opt;
      opt.add(desc);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("wavefunction") == 0 || HamStr.empty())
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-idmrg [options]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      if (vm.count("one-site"))
	 TwoSite = !OneSite;

      SimpleMPOperator HamMPO;
      if (HamStr == "itf")
      {
	 std::cout << "Hamiltonian is transverse-field Ising, J=" << J << ", Lambda=" << Lambda << "\n";
	 SiteBlock Site = CreateSpinSite(0.5);
	 MpOpTriangular Ham;
	 Ham = J * 4.0 * TriangularTwoSite(Site["Sz"], Site["Sz"])
	    + Lambda * 2.0 * TriangularOneSite(Site["Sx"]);
	 HamMPO.push_back(Ham.data());
      }
      else if (HamStr == "itf-z2")
      {
	 std::cout << "Hamiltonian is transverse-field Ising with Z2, J=" << J << ", Lambda=" << Lambda << "\n";
	 SiteBlock Site = CreateZ2SpinSite(0.5);
	 MpOpTriangular Ham;
	 Ham = J * 4.0 * TriangularTwoSite(Site["Sz"], Site["Sz"])
	    + Lambda * 2.0 * TriangularOneSite(Site["Sx"]);
	 HamMPO.push_back(Ham.data());
      }

      InfiniteWavefunction Psi;
      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pvalue_ptr<InfiniteWavefunction> PsiPtr = pheap::OpenPersistent(FName, CacheSize, true);
      Psi = *PsiPtr;
      
      StateComponent Phi = Psi.Psi.get_front();
      MatrixOperator LambdaSqrt = SqrtDiagonal(Psi.C_old);
      MatrixOperator LambdaInvSqrt = InvertDiagonal(LambdaSqrt, InverseTol);
      Phi = prod(prod(LambdaInvSqrt, Phi), LambdaSqrt);

      OperatorComponent Op = HamMPO.front();

      StateComponent Eigenvector(Op.Basis1(), Phi.Basis1(), Phi.Basis2());
      Eigenvector[0] = Psi.C_old;

      // symmetrize Phi
      Phi[0] = 0.5 * (Phi[0] + adjoint(Phi[0]));
      Phi[1] = 0.5 * (Phi[1] - adjoint(Phi[1]));

      TRACE(norm_frob(Phi[0] - adjoint(Phi[0])));
      TRACE(norm_frob(Phi[1] + adjoint(Phi[1])));

      double E = SolveFirstOrderMPO_Symmetric(Phi, Op, Eigenvector, Verbose);

      double PhiEnergy = E;
      TRACE(E);
      TRACE(trace(Eigenvector[0]));
      TRACE(norm_frob(Phi[0] - adjoint(Phi[0])));
      TRACE(norm_frob(Phi[1] + adjoint(Phi[1])));
      MatrixOperator Normalizer = Eigenvector[0];
      Normalizer = InvertDiagonal(Normalizer, InverseTol);
      //      InvertGeneral(Normalizer);

      // symmetrize Phi
      Phi[0] = 0.5 * (Phi[0] + adjoint(Phi[0]));
      Phi[1] = 0.5 * (Phi[1] - adjoint(Phi[1]));

      StateComponent RightVector(Eigenvector);
      double Timestep = 0.01;

      InverseTol = 1E-6;
      for (int Iter = 0; Iter < 10000; ++Iter)
      {
         for (unsigned i = 0; i < Op.Basis1().size(); ++i)
         {
            RightVector[i] = Eigenvector[Op.Basis1().size()-i-1];
         }
         RightVector[1] *= -0.25;
         Normalizer = adjoint(Eigenvector[0]);
               Normalizer = InvertDiagonal(Normalizer, InverseTol);
      //InvertGeneral(Normalizer);

         TRACE(SingularValues(Eigenvector[0] * Normalizer));


         {
            // This energy calculation appears to be correct!
            StateComponent Next = operator_prod_inner(Op, RightVector, Phi, herm(Eigenvector));
            std::complex<double> CalcEnergy = inner_prod(Phi, Next);
            TRACE(CalcEnergy);
	    E = CalcEnergy.real();
	    TRACE(E);
         }


#if 0
	 for (int k = 0; k < 1; ++k)
	    Phi = Phi - Timestep * triple_prod(Normalizer, operator_prod_inner(Op, RightVector, Phi, herm(Eigenvector)), herm(Normalizer));

#else
         {
	    TRACE(NormFrobNonOrtho<StateComponent>(Eigenvector[0])(Phi));


            int Iterations = 10;
            double Tol = 1E-14;
            StateComponent PhiSave = Phi;
            double Energy = Lanczos(Phi, OperatorProdInner_AxBH_Normal(Op, RightVector, Eigenvector, Normalizer), 
				    InnerProdNonOrtho<StateComponent>(Eigenvector[0]),
				    NormFrobNonOrtho<StateComponent>(Eigenvector[0]),
				    Iterations, Tol);
            //            double Energy = Lanczos(Phi, OperatorProdInner_AxBH_Gen(Op, RightVector, Eigenvector, Eigenvector[0]), Iterations, Tol);
            TRACE(Energy)(Iterations)(Tol);

            if (inner_prod(PhiSave, triple_prod(Normalizer, Phi, herm(Normalizer))).real() < 0)
               Phi *= -1.0;

#if 0
            // shift to delta scheme
            Phi *= -1.0;
            Phi = Phi + PhiSave;
#endif

            double Epsilon = 0.00001;
            StateComponent G1 = PhiSave - Epsilon * Phi;
            StateComponent G2 = PhiSave + Epsilon * Phi;

            StateComponent Eigen1 = Eigenvector;
            StateComponent Eigen2 = Eigenvector;

            double E1 =  SolveFirstOrderMPO_Symmetric(G1, Op, Eigen1, Verbose);
            double E2 =  SolveFirstOrderMPO_Symmetric(G2, Op, Eigen2, Verbose);
            
            // polynomial fit for (E1, PhiEnergy, E2)
            // E(x) = E0 + ax + bx^2
            // where a = (E1 - E2) / (2*epsilon)
            // b = (E1+E2-2*E0)/(epsilon^2)
            // Minimum at x=-alpha/(2*beta), which has energy E0 - a^2/(4*b)

            TRACE(E1)(E2)(PhiEnergy);

            double a = (E2-E1) / (2*Epsilon);
            double b = (E1+E2-2*PhiEnergy) / (Epsilon*Epsilon);
            double x = -0.5 * a / b;
            double PredictedE = PhiEnergy - 0.25*a*a/b;

            TRACE(a)(b)(x)(PredictedE)(PredictedE-PhiEnergy);

            if (b < std::numeric_limits<double>::epsilon())
            {
               TRACE("negative curvature!");

               // make it so that the +ve epsilon direction moves downhill faster
               if (E1 < E2)
               {
                  std::swap(E1, E2);
                  Epsilon = -Epsilon;
               }

               StateComponent G3 = PhiSave + 2*Epsilon * Phi;
               StateComponent Eigen3 = Eigenvector;
               double E3 = SolveFirstOrderMPO_Symmetric(G3, Op, Eigen3, Verbose);

               double Curvature = (PhiEnergy+E3-2*E2) / (Epsilon*Epsilon);
               TRACE(Curvature)(Epsilon)(E3)(E2)(PhiEnergy);
               while (Curvature < 0)
               {
                  TRACE("Still negative!");

                  Epsilon = Epsilon*2;
                  E2=E3;
                  G2=G3;
                  Eigen2=Eigen3;

                  G3 = PhiSave + 2*Epsilon * Phi;
                  E3 = SolveFirstOrderMPO_Symmetric(G3, Op, Eigen3, Verbose);
                  Curvature = (PhiEnergy+E3-2*E2) / (Epsilon*Epsilon);
                  TRACE(Curvature)(Epsilon)(E3)(E2)(PhiEnergy);
               }

               a = (E3-PhiEnergy) / (2*Epsilon);
               b = Curvature;
               PredictedE = E2 - 0.25*a*a/b;
               x = Epsilon - 0.5 * a / b;

               TRACE(a)(b)(x)(PredictedE)(PredictedE-E2);
            }               

            Phi = PhiSave + x*Phi;

#if 0
            // successive over-relaxation step
            double const w = 0.0003;
            Phi = w*Phi + (1.0-w)*PhiSave;

	    E = Energy;
	    TRACE(E);
#endif


      TRACE(norm_frob(Phi[0] - adjoint(Phi[0])));
      TRACE(norm_frob(Phi[1] + adjoint(Phi[1])));

      // symmetrize Phi
      Phi[0] = 0.5 * (Phi[0] + adjoint(Phi[0]));
      Phi[1] = 0.5 * (Phi[1] - adjoint(Phi[1]));


         }
#endif

         E = SolveFirstOrderMPO_Symmetric(Phi, Op, Eigenvector, Verbose);
         PhiEnergy = E;
         TRACE(E);
      TRACE(norm_frob(Phi[0] - adjoint(Phi[0])));
      TRACE(norm_frob(Phi[1] + adjoint(Phi[1])));

#if 0
      Timestep *= 0.999;
      if (Iter > 100)
         Timestep = 0.01;
      TRACE(Timestep);
#endif
      }

      pheap::Shutdown();

      ProcControl::Shutdown();
      return 0;
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
