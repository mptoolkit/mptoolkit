
#include "mps/infinitewavefunction.h"
#include "tensor/tensor_eigen.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include "common/environment.h"
#include "interface/inittemp.h"
#include "mp/copyright.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try 
   {
      int MaxEigenvalues = 10000;
      bool Base2 = false;
      bool ShowDegen = false;
      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("entropy,e", "show the entropy at each partition")
         ("states,s", "show the number of states at each partition")
	 ("basis,a", "show the complete basis at each partition")
         ("localbasis,b", "show the local basis at each site")
         ("density-matrix,d", "show the density matrix eigenvalues")
         ("degen", prog_opt::bool_switch(&ShowDegen),
          "Show degeneracies in the density matrix as repeated eigenvalues (only with -d)")
	 ("casimir,c", "show the values of the casimir invariant operators at each partition")
	 ("trans,t", "calculate left/right eigenvectors of the transfer operator and show how far they deviate from unity")
         ("limit,l", prog_opt::value<int>(&MaxEigenvalues), 
          "limit the density matrix display to N eigenvalues (implies --density-matrix)")
	 ("base2,2", prog_opt::bool_switch(&Base2), "show the entropy using base 2 instead of base e")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("input-wavefunction", prog_opt::value<std::string>(), "input wavefunction (required)")
         ;

      prog_opt::positional_options_description p;
      p.add("input-wavefunction", -1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);    
      
      if (vm.count("help") || vm.count("input-wavefunction") == 0) 
      {
         print_copyright(std::cerr);
         std::cerr << "usage: mp-info [options] input-wavefunction\n";
         std::cerr << desc << "\n";
         return 1;
      }
      
      bool ShowEntropy = vm.count("entropy");
      bool ShowStates = vm.count("states");
      bool ShowLocalBasis = vm.count("localbasis");
      bool ShowBasis = vm.count("basis");
      bool ShowCasimir = vm.count("casimir");
      bool ShowDensity = vm.count("density-matrix") || vm.count("limit");
      bool ShowTrans = vm.count("trans");

      std::string Wavefunc = vm["input-wavefunction"].as<std::string>();

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));

      pvalue_ptr<InfiniteWavefunction> Psi = pheap::OpenPersistent(Wavefunc, mp_pheap::CacheSize(), true);

      std::cout << "Symmetry list=" << Psi->C_right.GetSymmetryList() << '\n';
      std::cout << "Transforms as=" << Psi->shift() << '\n';
      
      std::cout << "Number of states=" << Psi->C_right.Basis1().total_dimension() << '\n';
      std::cout << "Degree=" << Psi->C_right.Basis1().total_degree() << '\n';
      std::cout << "Unit cell size=" << Psi->Psi.size() << '\n';
      
      std::cout << "Orthogonality fidelity=" << (1.0 - orthogonality_fidelity(*Psi)) << '\n';

      MatrixOperator Rho = scalar_prod(Psi->C_right, herm(Psi->C_right));
      DensityMatrix<MatrixOperator> DM(Rho);
      double Entropy1 = DensityEntropy(DM.begin(), DM.end(), DM.EigenSum(), Base2);
      std::cout << "Entropy1=" << Entropy1 << '\n';

      MatrixOperator Rho0 = scalar_prod(Psi->C_old, herm(Psi->C_old));
      DensityMatrix<MatrixOperator> DM0(Rho0);
      double Entropy0 = DensityEntropy(DM0.begin(), DM0.end(), DM0.EigenSum(), Base2);
      std::cout << "Entropy0=" << Entropy0 << '\n';

      // rotate through the unit cell and calculate the entropy
      InfiniteWavefunction P = *Psi;

      if (ShowTrans)
      {
	 MatrixOperator I = MatrixOperator::make_identity(Psi->Psi.Basis1());
	 MatrixOperator J = inject_left(I, Psi->Psi);
	 TRACE(norm_frob(I-J));

	 LinearWavefunction PsiR = Psi->Psi;
	 MatrixOperator Lambda = Psi->C_right;
	 MatrixOperator LambdaInv = InvertDiagonal(Lambda, InverseTol);
	 PsiR.set_front(prod(LambdaInv, PsiR.get_front()));
	 PsiR.set_back(prod(PsiR.get_back(), Lambda));

	 I = MatrixOperator::make_identity(PsiR.Basis2());
	 J = inject_right(I, PsiR);
	 TRACE(norm_frob(I-J));
      }

      if (ShowStates || ShowLocalBasis)
      {
	 std::cout << "quantities in the unit cell:\n";
         int i=0;
         for (LinearWavefunction::const_iterator I = P.Psi.begin();
              I != P.Psi.end(); ++I)
         {
            std::cout << "Site " << i++ << '\n';
            if (ShowStates)
            {
               std::cout << "LeftStates: " << I->Basis1().total_dimension()
                         << " RightStates: " << I->Basis2().total_dimension()
                         << "\n";
            }
            if (ShowLocalBasis)
            {
               std::cout << "LocalBasis: " << I->LocalBasis() << '\n';
            }
         }
      }

      if (ShowCasimir || ShowBasis || ShowEntropy || ShowDensity)
      {

	 std::cout << "quantities in the unit cell:\n";
	 if (ShowCasimir)
	 {
	    SymmetryList const SList = P.GetSymmetryList();
	    int const NumCasimir = SList.NumCasimirOperators();
	    for (int i = 0; i < NumCasimir; ++i)
	    {
	       if (i != 0)
		  std::cout << " ";
	       std::string Name = "#" + SList.CasimirName(i);
	       std::cout << std::setw(20) << Name;
	    }
	    std::cout << '\n';
	 }

	 for (int i = 0; i < P.size(); ++i)
	 {
	    MatrixOperator Rho = scalar_prod(P.C_right, herm(P.C_right));

	    if (ShowEntropy)
	    {
	       DensityMatrix<MatrixOperator> DM(Rho);
	       double Entropy = DensityEntropy(DM.begin(), DM.end(), DM.EigenSum(), Base2);
	       std::cout << Entropy << '\n';
	    }

	    if (ShowBasis)
	    {
	       std::cout << Rho.Basis1() << '\n';
	    }

	    if (ShowCasimir)
	    {
	       SymmetryList const SList = P.GetSymmetryList();
	       int const NumCasimir = SList.NumCasimirOperators();
	       DensityMatrix<MatrixOperator> DM(Rho);
	       for (int i = 0; i < NumCasimir; ++i)
	       {
		  if (i != 0) 
		     std::cout << ' ';
		  std::cout << std::setw(20) << DM.EvaluateCasimir(i);
	       }
	       std::cout << '\n';
	    }

	    if (ShowDensity)
	    {
	       DensityMatrix<MatrixOperator> DM(Rho);
	       DM.DensityMatrixReport(std::cout, MaxEigenvalues, Base2, ShowDegen);
	    }

            // only rotate if we're not on the final iteration
            if (i+1 < P.size())
            {
               P = rotate_left(P, 1);
            }
	 }
      }

   }
   catch(std::exception& e) 
   {
      std::cerr << "error: " << e.what() << "\n";
      return 1;
   }
   catch(...) 
   {
      std::cerr << "Exception of unknown type!\n";
   }
   return 0;
}
