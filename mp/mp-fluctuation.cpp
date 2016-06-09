// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-fluctuation.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Reseach publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#include "wavefunction/mpwavefunction.h"
#include "mps/packunpack.h"
#include "interface/inittemp.h"
#include "lattice/latticesite.h"
#include "wavefunction/operator_actions.h"
#include "mp/copyright.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include "common/prog_opt_accum.h"
#include "common/environment.h"
#include "interface/inittemp.h"
#include "tensor/tensor_eigen.h"
#include "mp-algorithms/arnoldi.h"
#include "lattice/unitcell.h"
#include "lattice/infinite-parser.h"
#include "linearalgebra/arpack_wrapper.h"
#include <boost/algorithm/string/predicate.hpp>
#include "mp-algorithms/triangular_mpo_solver.h"
#include "common/statistics.h"
#include <tuple>
#include "parser/matrix-parser.h"

namespace prog_opt = boost::program_options;

template <typename Func>
struct PackApplyFunc
{
   PackApplyFunc(PackStateComponent const& Pack_, Func f_) : Pack(Pack_), f(f_) {}

   void operator()(std::complex<double> const* In, std::complex<double>* Out) const
   {
      StateComponent x = Pack.unpack(In);
      x = f(x);
      Pack.pack(x, Out);
   } 
   PackStateComponent const& Pack;
   Func f;
};

template <typename Func>
PackApplyFunc<Func>
MakePackApplyFunc(PackStateComponent const& Pack_, Func f_)
{
   return PackApplyFunc<Func>(Pack_, f_);
}

std::tuple<std::complex<double>, int, StateComponent>
get_left_eigenvector(LinearWavefunction const& Psi1, QuantumNumber const& QShift1,
		     LinearWavefunction const& Psi2, QuantumNumber const& QShift2, 
                     ProductMPO const& StringOp,
                     double tol = 1E-14, int Verbose = 0)
{
   int ncv = 0;
   int Length = statistics::lcm(Psi1.size(), Psi2.size(), StringOp.size());
   PackStateComponent Pack(StringOp.Basis1(), Psi1.Basis2(), Psi2.Basis2());
   int n = Pack.size();
   //   double tolsave = tol;
   //   int ncvsave = ncv;
   int NumEigen = 1;

   std::vector<std::complex<double> > OutVec;
      LinearAlgebra::Vector<std::complex<double> > LeftEigen = 
         LinearAlgebra::DiagonalizeARPACK(MakePackApplyFunc(Pack,
							    LeftMultiplyOperator(Psi1, QShift1,
										 StringOp, 
										 Psi2, QShift2, Length, Verbose-2)),
					  n, NumEigen, tol, &OutVec, ncv, false, Verbose);

   StateComponent LeftVector = Pack.unpack(&(OutVec[0]));
      
   return std::make_tuple(LeftEigen[0], Length, LeftVector);
}

void PrintFormat(std::complex<double> const& Value, bool ShowRealPart, bool ShowImagPart, 
		 bool ShowMagnitude, bool ShowArgument,
		 bool ShowRadians)
{
   if (ShowRealPart)
   {
      std::cout << std::setw(20) << Value.real() << "    ";
   }
   if (ShowImagPart)
   {
      std::cout << std::setw(20) << Value.imag() << "    ";
   }
   if (ShowMagnitude)
   {
      std::cout << std::setw(20) << LinearAlgebra::norm_frob(Value) << "    ";
   }
   if (ShowArgument)
   {
      double Arg =  std::atan2(Value.imag(), Value.real());
      if (!ShowRadians)
	 Arg *= 180.0 / math_const::pi;
      std::cout << std::setw(20) << Arg << "    ";
   }
}

// we assume that Rho is diagonal here
void
ShowValuesBySector(MatrixOperator const& v, RealDiagonalOperator const& Rho,
		   bool Quiet, bool ShowRealPart, bool ShowImagPart, 
		   bool ShowMagnitude, bool ShowArgument, bool ShowRadians)
{
   if (!Quiet)
   {
      std::cout << "#sector    #weight              ";
      if (ShowRealPart)
	 std::cout << "#real                   ";
      if (ShowImagPart)
	 std::cout << "#imag                   ";
      if (ShowMagnitude)
	 std::cout << "#magnitude              ";
      if (ShowArgument)
	 std::cout << "#argument" << (ShowRadians ? "(rad)" : "(deg)") << "          ";
      std::cout << '\n';
   }

   std::map<QuantumNumber, std::complex<double>> ValueBySector;
   std::map<QuantumNumber, double> WeightBySector;
   for (unsigned i = 0; i < Rho.Basis1().size(); ++i)
   {
      MatrixOperator Component(Rho.Basis1(), Rho.Basis2(), Rho.TransformsAs());
      Component(i,i) = Rho(i,i);
      ValueBySector[Rho.Basis1()[i]] += inner_prod(v, Component);
      // Need to include the degree in the weight by hand, since Rho(i,i) loses the quantum number information
      WeightBySector[Rho.Basis1()[i]] += trace(Rho(i,i)) * degree(Rho.Basis1()[i]);
   }
   for (auto const& x : ValueBySector)
   {
      double Weight = WeightBySector[x.first];
      std::cout << std::setw(10) << std::left << x.first.ToString() << ' ';
      std::cout << std::setw(20) << Weight << ' ';
      PrintFormat(x.second / Weight, 
		  ShowRealPart, ShowImagPart, ShowMagnitude, 
		  ShowArgument, ShowRadians);
      std::cout << '\n';
   }
}

int main(int argc, char** argv)
{
   try
   {
      std::string PsiStr;
      double Tol = 1e-14;
      int Verbose = 0;
      bool ShowRealPart = false;
      bool ShowImagPart = false;
      bool ShowMagnitude = false;
      bool ShowArgument = false;
      bool ShowRadians = false;
      bool ShowPolar = false;
      bool ShowCartesian = false;
      bool Quiet = false;
      bool UseTempFile = false;
      bool NoRandomPhase = false;
      std::vector<std::string> ProductOperators;
      std::vector<std::string> TriangularOperators;
      std::vector<int> PolyDegree;
      int Partition = 0;  // not yet implemented
      double UnityEpsilon = DefaultEigenUnityEpsilon;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("product,p", prog_opt::value(&ProductOperators), "parse a product MPO")
	 ("triangular,t", prog_opt::value(&TriangularOperators), "parse a triangular MPO")
	 ("cart,c", prog_opt::bool_switch(&ShowCartesian),
	  "show the result in cartesian coordinates [equivalent to --real --imag]")
	 ("polar", prog_opt::bool_switch(&ShowPolar),
	  "show the result in polar coodinates [equivalent to --mag --arg]")
	 ("real,r", prog_opt::bool_switch(&ShowRealPart),
	  "display the real part of the result")
	 ("imag,i", prog_opt::bool_switch(&ShowImagPart),
	  "display the imaginary part of the result")
         ("mag,m", prog_opt::bool_switch(&ShowMagnitude),
          "display the magnitude of the result")
         ("arg,a", prog_opt::bool_switch(&ShowArgument),
          "display the argument (angle) of the result")
	 ("radians", prog_opt::bool_switch(&ShowRadians),
	  "display the argument in radians instead of degrees")
	 ("degree", prog_opt::value(&PolyDegree), 
	  "for a TriangularMPO, only show the results for the terms of this degree [can be used more than once]")
         ("tol", prog_opt::value(&Tol),
          FormatDefault("Tolerance of the Arnoldi eigensolver", Tol).c_str())
	 ("unityepsilon", prog_opt::value(&UnityEpsilon),
	  FormatDefault("Epsilon value for testing eigenvalues near unity", UnityEpsilon).c_str())
	 ("norandomize", prog_opt::bool_switch(&NoRandomPhase), "for ProductMPO's, don't randomize the phase of the eigenvectors")
	 ("quiet,q", prog_opt::bool_switch(&Quiet), "suppress headings on the output")
	 ("verbose,v", prog_opt_ext::accum_value(&Verbose), "increase verbosity")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("wavefunction,w", prog_opt::value(&PsiStr), "Wavefunction [required]")
         ;

      prog_opt::positional_options_description p;
      p.add("wavefunction", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);    

      if (vm.count("help") > 0 || vm.count("wavefunction") < 1)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] -w <psi> [-t <triangular-MPO>] [-p <product-MPO] ...\n";
         std::cerr << desc << '\n';
	 return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      // If no output switches are used, default to showing everything
      if (!ShowRealPart && !ShowImagPart && !ShowMagnitude
	  && !ShowCartesian && !ShowPolar && !ShowArgument)
      {
	 ShowCartesian = true;
	 ShowPolar = true;
      }
      
      if (ShowCartesian)
      {
	 ShowRealPart = true;
	 ShowImagPart = true;
      }
      if (ShowPolar)
      {
	 ShowMagnitude = true;
	 ShowArgument = true;
      }
      if (ShowRadians)
	 ShowArgument = true;

      if (!Quiet)
	 print_preamble(std::cout, argc, argv);

      if (!NoRandomPhase)
	 srand(ext::get_unique() % RAND_MAX);

      // if the number of eigenvalues is specified but
      // the cutoff is not, then set the cutoff to zero
      //      if (vm.count("num-eigenvalues") == 1 && vm.count("eigen-cutoff") == 0)
      //      {
      //         EigenCutoff = 0;
      //      }


      // Load the wavefunction
      pvalue_ptr<MPWavefunction> Psi;
      if (UseTempFile)
      {
	  mp_pheap::InitializeTempPHeap(Verbose);
	  Psi = pheap::ImportHeap(PsiStr);
      }
      else
      {
         Psi = pheap::OpenPersistent(PsiStr, mp_pheap::CacheSize(), true);
      }
      InfiniteWavefunctionLeft InfPsi = Psi->get<InfiniteWavefunctionLeft>();

      LinearWavefunction Psi1;
      RealDiagonalOperator D;
      std::tie(Psi1, D) = get_left_canonical(InfPsi);
      RealDiagonalOperator RhoDiag = D*D;
      MatrixOperator Rho = D;
      Rho = scalar_prod(Rho, herm(Rho));
      MatrixOperator Identity = MatrixOperator::make_identity(Psi1.Basis1());
      double Dim = Psi1.Basis1().total_degree();

      bool First = true;

      // Make a set out of the PolyDegree
      std::set<int> ShowPolyDegree(PolyDegree.begin(), PolyDegree.end());

      // process the product operators
      for (unsigned i = 0; i < ProductOperators.size(); ++i)
      {
	 if (!First)
	    std::cout << '\n';
	 First = false;
	 LinearWavefunction* Psi2 = &Psi1;
	 if (!Quiet)
	    std::cout << "#Product Operator " << ProductOperators[i] << '\n';
	 ProductMPO Op;
	 InfiniteLattice Lattice;
	 std::tie(Op, Lattice) = ParseProductOperatorAndLattice(ProductOperators[i]);     

         std::complex<double> e;
         StateComponent v;
	 int n;
         std::tie(e, n, v) = get_left_eigenvector(Psi1, InfPsi.qshift(), *Psi2, InfPsi.qshift(), Op,
						  Tol, Verbose);

	 // Normalization
	 // it might not be unitary, eg anti-unitary.  So we need to take the 4th power
	 std::complex<double> x = inner_prod(scalar_prod(herm(v), operator_prod(herm(v), v, v)), Rho);

	 v *= std::sqrt(std::sqrt(1.0 / x));

	 // randomize phase
	 if (!NoRandomPhase)
	 {
	    v *= std::polar(1.0, LinearAlgebra::random<double>() * 2 * math_const::pi);
	 }

	 ShowValuesBySector(v[0], RhoDiag, Quiet, ShowRealPart, ShowImagPart, ShowMagnitude, ShowArgument, ShowRadians);
      }

      // process the triangular operators

      for (unsigned i = 0; i < TriangularOperators.size(); ++i)
      {
	 if (!First)
	    std::cout << '\n';
	 First = false;
	 if (!Quiet)
	    std::cout << "#Triangular Operator " << TriangularOperators[i] << '\n';
	 TriangularMPO Op;
	 InfiniteLattice Lattice;
	 std::tie(Op, Lattice) = ParseTriangularOperatorAndLattice(TriangularOperators[i]);     

	 Op = repeat(Op, Psi1.size() / Op.size());

	 std::vector<KMatrixPolyType> E;
	 SolveMPO_Left(E, Psi1, InfPsi.qshift(), Op, Identity, Rho, true, Tol, UnityEpsilon, Verbose);
	 
	 Polynomial<MatrixOperator> v = E.back()[1.0];
	 
	 for (int d = 0; d <= v.degree(); ++d)
	 {
	    if (v.has_term(d) && (ShowPolyDegree.empty() || (ShowPolyDegree.find(d) != ShowPolyDegree.end())))
	    {
	       std::cout << "#degree " << d << '\n';
	       ShowValuesBySector(v[d], RhoDiag, Quiet, ShowRealPart, ShowImagPart, ShowMagnitude, ShowArgument, ShowRadians);
	    }
	 }
      }

      pheap::Shutdown();
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
