// -*- C++ -*- $Id$

#include "mp/copyright.h"
#include "wavefunction/mpwavefunction.h"
#include <boost/program_options.hpp>
#include "common/environment.h"
#include "interface/inittemp.h"
#include "common/terminal.h"
#include "common/environment.h"
#include "common/prog_options.h"
#include "lattice/unitcell.h"
#include "lattice/infinite-parser.h"
#include "common/statistics.h"
#include "common/formatting.h"
#include "tensor/tensor_eigen.h"

using QuantumNumbers::QuantumNumber;

namespace prog_opt = boost::program_options;

#if 0
InfiniteWavefunctionLeft take_subset(InfiniteWavefunctionLeft const& Psi, 
				     MatrixOperator const& U, int first, int last)
{
   CHECK(first >= 0 && first < last && last <= Psi.size());
   InfiniteWavefunctionLeft Result;

   InfiniteWavefunctionLeft::const_lambda_iterator iL = Psi.lambda_begin() + first;
   InfiniteWavefunctionLeft::const_mps_iterator iA = Psi.begin() + first;
   for (int i = first; i < last-1; ++i, ++iL, ++iA)
   {
      Result.push_back_lambda(*iL);
      Result.push_back(*iA);
   }

   Result.push_back_lambda(*iL); 
   ++iL;

   // for the last site, we wrap around to the first using U
   Result.push_back(prod(*iA, herm(U)));
   // and re-use the first lambda
   Result.push_back_lambda(Result.lambda(0));

   Result.setBasis1(U.Basis1());
   Result.setBasis2(U.Basis1());

   Result.check_structure();

   return Result;
}
#endif

int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      std::string InputFile;
      std::string OutputFile;
      int Iter = 30;
      double Tol = 1E-15;
      int NewUnitCellSize = 0;
      int Div = 2;
      bool Force = false;
      bool Quiet = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("force,f", prog_opt::bool_switch(&Force),
	  "overwrite the output file, if it exists")
	 ("divide,d", prog_opt::value(&Div), "divide the wavefunction into this many parts [default 2]")
	 ("unitcell,u", prog_opt::value(&NewUnitCellSize), "new unit cell size (alternative to --divide)")
         ("tol", prog_opt::value(&Tol),
          FormatDefault("Tolerance of the Arnoldi eigensolver", Tol).c_str())
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose),
          "extra debug output [can be used multiple times]")
         ("quiet", prog_opt::bool_switch(&Quiet), "don't print normal information messages")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("lhs", prog_opt::value<std::string>(&InputFile), "psi1")
         ("rhs", prog_opt::value<std::string>(&OutputFile), "psi2")
         ;

      prog_opt::positional_options_description p;
      p.add("lhs", 1);
      p.add("rhs", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);    

      if (vm.count("help") > 0 || vm.count("rhs") == 0)
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options] <input-psi> <output-psi>\n";
         std::cerr << desc << '\n';
	 std::cerr << "This tool is used to divide the unit cell of an iMPS of unit cell size N,\n"
		   << "in the case where it is translationally invariant under a smaller unit cell size M\n"
		   << "(which will exactly divide N).  This means that the unit cell can be reduced in size\n"
		   << "with no loss of information.  This requires that the wavefunction is exactly translationally\n"
		   << "invariant under an M-site shift; that is, mp-ioverlap --rotate M psi psi gives numerically 1.0.\n";
         return 1;
      }

      // it is an error to specify both --divide and --unitcell
      if (vm.count("divide") && vm.count("unitcell"))
      {
	 std::cerr << basename(argv[0]) << ": fatal: --divide and --unitcell cannot both be specified.\n";
	 return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      // create a new file for output
      pheap::Initialize(OutputFile, 1, mp_pheap::PageSize(), mp_pheap::CacheSize(), false, Force);
      // and load the input wavefunction
      pvalue_ptr<MPWavefunction> InputPsi = pheap::ImportHeap(InputFile);

      InfiniteWavefunctionLeft Psi1 = InputPsi->get<InfiniteWavefunctionLeft>();

      // determine the size of the unit cell, and the new unit cell
      int UnitCellSize = Psi1.size();
      // if --unitcell was specified, then use that number.  Otherwise, use --divide (which has a default value)
      if (vm.count("unitcell"))
      {
	 if (UnitCellSize % NewUnitCellSize != 0)
	 {
	    std::cerr << basename(argv[0]) << ": fatal: the input unit cell is not a multiple of the specified output unit cell.\n";
	    return 1;
	 }
      } 
      else
      {
	 if (UnitCellSize % Div != 0)
	 {
	    std::cerr << basename(argv[0]) << ": fatal: the input unit cell is not a multiple of the number of divisions.\n";
	    return 1;
	 }
	 NewUnitCellSize = UnitCellSize / Div;
      }

      InfiniteWavefunctionLeft Psi2 = Psi1;
      Psi2.rotate_left(NewUnitCellSize);
      // TODO: handle the qshift!


      // the sector is always the identity (**after qshift**!)
      QuantumNumber Ident = QuantumNumber(Psi1.GetSymmetryList());

      std::complex<double> e;
      StateComponent Vec;
      std::tie(e, Vec) = overlap(Psi1, Psi2, Ident, Iter, Tol, Verbose);
      
      if (!Quiet)
      {
	 std::cout << "Overlap eigenvalue is " << format_complex(e) << '\n';
      }

      // Check that e is close to 1.0
      if (norm_frob(e - 1.0) > 0.1)
      {
	 std::cerr <<  basename(argv[0]) << ": warning: overlap eigenvalue " << format_complex(e) 
		   << " is not close to 1.0, results are surely not reliable!\n";
      }
      else if (norm_frob(e - 1.0) > 1E-5)
      {
	 std::cerr <<  basename(argv[0]) << ": warning: overlap eigenvalue " << format_complex(e) 
		   << " is slightly different from 1.0, results might not be reliable.\n";
      }

      // Vec only has one component, in the identity sector.  It should be a unitary matrix.
      MatrixOperator X = Vec[0];

      // X is the eigenvector (with eigenvalue 1) that maps the Basis1() of Psi1 with the Basis1() of Psi2,
      // which is the same as Psi1[NewUnitCellSize].Basis1().  That is, we can use X to 'wrap around' Psi1 and
      // cut it off at NewUnitCellSize sites.

      // Force X to be unitary
      MatrixOperator U, Vh;
      RealDiagonalOperator D;

      SingularValueDecompositionFull(X, U, D, Vh);

      // the unitary part
      X = U*Vh;

      // See how close we were to unitary.  To do this, we check how similar the singular values are -
      // in principle they are all equal (but with some overall normalization). But we don't expect that
      // it will be perfectly unitary for components corresponding to small singular values. 
      // lambda U D^2 U^\dagger \lambda^\dagger
      // The trace of this gives the overall normalization
      // Then once D is normalized, U D^2 U^\dagger should be the identity, so subtract the identity
      // and calculate tr((U D^2 U^\dagger - I)^2 . rho)
      double x = norm_frob_sq(Psi1.lambda(0) * U * D);
      D *= (1.0 / std::sqrt(x));

      MatrixOperator CheckIdent = U*D*D*herm(U);
      CheckIdent = CheckIdent - MatrixOperator::make_identity(CheckIdent.Basis1());
      std::complex<double> diff = norm_frob_sq(Psi1.lambda(0) * CheckIdent);

      TRACE(x);

      TRACE(diff);

      // actually it is much easier than the above - U should map the two lambda matrices
      
      LinearWavefunction PsiNew(Psi1.begin(), Psi1.begin()+NewUnitCellSize);
      PsiNew.set_back(prod(PsiNew.get_back(), herm(X)));

      QuantumNumber QShift = Ident;

      InfiniteWavefunctionLeft PsiOut(PsiNew, QShift);
      Psi1 = InfiniteWavefunctionLeft();
      Psi2 = Psi1;

      MPWavefunction Result;
      Result.AppendHistory(EscapeCommandline(argc, argv));
      Result.Wavefunction() = PsiOut;

      pvalue_ptr<MPWavefunction> OutputPsi = new MPWavefunction(Result);

      pheap::ShutdownPersistent(OutputPsi);

   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << '\n';
      return 1;
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
