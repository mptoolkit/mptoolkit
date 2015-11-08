// -*- C++ -*- $Id$

#include "mp/copyright.h"
#include "wavefunction/mpwavefunction.h"
#include "lattice/latticesite.h"
#include "wavefunction/operator_actions.h"
#include "mps/packunpack.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include "common/prog_opt_accum.h"
#include "common/environment.h"
#include "interface/inittemp.h"
#include "tensor/tensor_eigen.h"
#include "mp-algorithms/arnoldi.h"
#include "lattice/siteoperator-parser.h"
#include "linearalgebra/arpack_wrapper.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcell-parser.h"
#include "lattice/infinite-parser.h"

namespace prog_opt = boost::program_options;


void PrintFormat(QuantumNumber const& q, std::complex<double> x, int n, bool ShowRealPart, bool ShowImagPart, 
		 bool ShowCorrLength, bool ShowMagnitude, bool ShowArgument,
		 bool ShowRadians, double ScaleFactor)
{
   std::string SectorStr = boost::lexical_cast<std::string>(q);
   std::complex<double> Value = std::pow(x, ScaleFactor);
   std::cout << std::setw(11) << SectorStr << ' ';
   std::cout << std::setw(4) << n << ' ';
   if (ShowRealPart)
   {
      std::cout << std::setw(20) << Value.real() << "  ";
   }
   if (ShowImagPart)
   {
      std::cout << std::setw(20) << Value.imag() << "  ";
   }
   if (ShowCorrLength)
   {
      std::cout << std::setw(20) << (-1.0/std::log(LinearAlgebra::norm_frob(Value)))
		<< "  ";
   }
   if (ShowMagnitude)
   {
      std::cout << std::setw(20) << LinearAlgebra::norm_frob(Value) << "  ";
   }
   if (ShowArgument)
   {
      double Arg =  std::atan2(Value.imag(), Value.real());
      if (!ShowRadians)
	 Arg *= 180.0 / math_const::pi;
      std::cout << std::setw(20) << Arg << "  ";
   }
}


struct LeftMultiply
{
   typedef MatrixOperator argument_type;
   typedef MatrixOperator result_type;

   LeftMultiply(LinearWavefunction const& L_, QuantumNumber const& QShift_) 
      : L(L_), QShift(QShift_) {}

   result_type operator()(argument_type const& x) const
   {
      result_type r = delta_shift(x, QShift);
      for (LinearWavefunction::const_iterator I = L.begin(); I != L.end(); ++I)
      {
	 r = operator_prod(herm(*I), r, *I);
      }
      return r;
   }

   LinearWavefunction const& L;
   QuantumNumber QShift;
};

struct LeftMultiplyStringSimple
{
   typedef MatrixOperator argument_type;
   typedef MatrixOperator result_type;

   LeftMultiplyStringSimple(LinearWavefunction const& L_, QuantumNumber const& QShift_,
                      ProductMPO const& StringOp_) 
      : L(L_), QShift(QShift_), StringOp(StringOp_) 
   {
      CHECK_EQUAL(int(L.size()), StringOp.size())("The string operator must be the same length as the unit cell");
   }

   result_type operator()(argument_type const& x) const
   {
      result_type r = delta_shift(x, QShift);
      r = inject_left(r, StringOp, L);
      return r;
   }

   LinearWavefunction const& L;
   QuantumNumber QShift;
   ProductMPO StringOp;
};
   
struct RightMultiply
{
   typedef MatrixOperator argument_type;
   typedef MatrixOperator result_type;

   RightMultiply(LinearWavefunction const& L_, QuantumNumber const& QShift_) 
      : L(L_), QShift(QShift_) {}

   result_type operator()(argument_type const& x) const
   {
      result_type r = x; 
      LinearWavefunction::const_iterator I = L.end();
      while (I != L.begin())
      {
         --I;
	 r = operator_prod(*I, r, herm(*I));
      }
      return delta_shift(r, adjoint(QShift));
   }

   LinearWavefunction const& L;
   QuantumNumber QShift;
};

struct RightMultiplyString
{
   typedef MatrixOperator argument_type;
   typedef MatrixOperator result_type;

   RightMultiplyString(LinearWavefunction const& L_, QuantumNumber const& QShift_,
                       ProductMPO const& StringOp_) 
      : L(L_), QShift(QShift_) , StringOp(StringOp_)
   {
      CHECK_EQUAL(int(L.size()), StringOp.size())("The string operator must be the same length as the unit cell");
   }

   result_type operator()(argument_type const& x) const
   {
      return inject_right_qshift(x, StringOp, L, QShift);
   }

   LinearWavefunction const& L;
   QuantumNumber QShift;
   ProductMPO StringOp;
};

struct RightMultiplyOp
{
   typedef MatrixOperator argument_type;
   typedef MatrixOperator result_type;

   RightMultiplyOp(LinearWavefunction const& L_, QuantumNumber const& QShift_, SimpleOperator const& Op_) 
      : L(L_), QShift(QShift_), Op(Op_) {}

   result_type operator()(argument_type const& x) const
   {
      result_type r = x; 
      LinearWavefunction::const_iterator I = L.end();
      while (I != L.begin())
      {
         --I;
	 r = operator_prod(Op, *I, r, herm(*I));
      }
      return delta_shift(x, adjoint(QShift));
   }

   LinearWavefunction const& L;
   QuantumNumber QShift;
   SimpleOperator const& Op;
};

struct MultFunc
{
   MultFunc(LinearWavefunction const& Psi, QuantumNumber const& QShift, 
            QuantumNumbers::QuantumNumber const& q)
      : Mult(Psi, QShift), Pack(Psi.Basis2(), Psi.Basis2(), q) {}

   void operator()(std::complex<double> const* In, std::complex<double>* Out) const
   {
      MatrixOperator x = Pack.unpack(In);
      x = Mult(x);
      Pack.pack(x, Out);
   }

   LeftMultiply Mult;
   PackMatrixOperator Pack;
};

struct MultFuncTrans
{
   MultFuncTrans(LinearWavefunction const& Psi, QuantumNumber const& QShift, 
            QuantumNumbers::QuantumNumber const& q)
      : Mult(Psi, QShift), Pack(Psi.Basis2(), Psi.Basis2(), q) {}

   void operator()(std::complex<double> const* In, std::complex<double>* Out) const
   {
      MatrixOperator x = Pack.unpack(In);
      x = Mult(x);
      Pack.pack(x, Out);
   }

   RightMultiply Mult;
   PackMatrixOperator Pack;
};

struct MultFuncString
{
   MultFuncString(LinearWavefunction const& Psi, QuantumNumber const& QShift, 
                  QuantumNumbers::QuantumNumber const& q,
                  ProductMPO const& StringOp)
      : Mult(Psi, QShift, StringOp), Pack(Psi.Basis2(), Psi.Basis2(), q) { }

   void operator()(std::complex<double> const* In, std::complex<double>* Out) const
   {
      MatrixOperator x = Pack.unpack(In);
      x = Mult(x);
      Pack.pack(x, Out);
   }

   LeftMultiplyStringSimple Mult;
   PackMatrixOperator Pack;
};

struct MultFuncStringTrans
{
   MultFuncStringTrans(LinearWavefunction const& Psi, QuantumNumber const& QShift, 
                       QuantumNumbers::QuantumNumber const& q,
                       ProductMPO const& StringOp)
      : Mult(Psi, QShift, StringOp), Pack(Psi.Basis2(), Psi.Basis2(), q) { }

   void operator()(std::complex<double> const* In, std::complex<double>* Out) const
   {
      MatrixOperator x = Pack.unpack(In);
      x = Mult(x);
      Pack.pack(x, Out);
   }

   RightMultiplyString Mult;
   PackMatrixOperator Pack;
};


// inject_left for a FiniteMPO.  This can have support on multiple wavefunction unit cells
MatrixOperator
inject_left(MatrixOperator const& m, 
            LinearWavefunction const& Psi1,
	    QuantumNumbers::QuantumNumber const& QShift,
            FiniteMPO const& Op, 
            LinearWavefunction const& Psi2)
{
   CHECK_EQUAL(Psi1.size(), Psi2.size());
   DEBUG_CHECK_EQUAL(m.Basis1(), Psi1.Basis2());
   DEBUG_CHECK_EQUAL(m.Basis2(), Psi2.Basis2());
   if (Op.is_null())
      return MatrixOperator();

   // we currently only support simple irreducible operators
   CHECK_EQUAL(Op.Basis1().size(), 1);
   CHECK_EQUAL(Op.Basis2().size(), 1);
   CHECK_EQUAL(Op.Basis1()[0], m.TransformsAs());
   MatrixOperator mm = delta_shift(m, QShift);
   StateComponent E(Op.Basis1(), mm.Basis1(), mm.Basis2());
   E[0] = mm;
   E.debug_check_structure();
   LinearWavefunction::const_iterator I1 = Psi1.begin();
   LinearWavefunction::const_iterator I2 = Psi2.begin();
   FiniteMPO::const_iterator OpIter = Op.begin();
   while (OpIter != Op.end())
   {
      if (I1 == Psi1.end())
      {
	 I1 = Psi1.begin();
	 I2 = Psi2.begin();
	 E = delta_shift(E, QShift);
      }
      E = contract_from_left(*OpIter, herm(*I1), E, *I2);
      ++I1; ++I2; ++OpIter;
   }
   return E[0]; //delta_shift(E[0], QShift);
}

MatrixOperator
inject_right(MatrixOperator const& m, 
             LinearWavefunction const& Psi1,
             QuantumNumbers::QuantumNumber const& QShift,
             GenericMPO const& Op, 
             LinearWavefunction const& Psi2)
{
   PRECONDITION_EQUAL(Psi1.size(), Psi2.size());
   if (Op.is_null())
      return MatrixOperator();

   // we currently only support simple irreducible operators
   CHECK_EQUAL(Op.Basis1().size(), 1);
   CHECK_EQUAL(Op.Basis2().size(), 1);
   StateComponent E(Op.Basis2(), m.Basis1(), m.Basis2());
   E[0] = m;
   MatrixOperator Result = m;
   LinearWavefunction::const_iterator I1 = Psi1.end();
   LinearWavefunction::const_iterator I2 = Psi2.end();
   GenericMPO::const_iterator OpIter = Op.end();
   while (OpIter != Op.begin())
   {
      if (I1 == Psi1.begin())
      {
         I1 = Psi1.end();
         I2 = Psi2.end();
         E = delta_shift(E, adjoint(QShift));
      }
      --I1; --I2; --OpIter;
      E = contract_from_right(herm(*OpIter), *I1, E, herm(*I2));
   }
   return delta_shift(E[0], adjoint(QShift));
}






LinearAlgebra::Vector<std::complex<double> > 
get_spectrum(LinearWavefunction const& Psi, QuantumNumber const& QShift, int NumEigen,
             QuantumNumbers::QuantumNumber const& q, double tol = 1e-10,
             LinearAlgebra::Vector<MatrixOperator>* RightVectors = NULL, 
             LinearAlgebra::Vector<MatrixOperator>* LeftVectors = NULL, 
             int ncv = 0, bool Sort = false, int Verbose = 0)
{
   PackMatrixOperator Pack(Psi.Basis2(), Psi.Basis2(), q);
   int n = Pack.size();
   double tolsave = tol;
   int ncvsave = ncv;

   std::vector<std::complex<double> >* OutVec 
      = RightVectors ? new std::vector<std::complex<double> >() : NULL;
   LinearAlgebra::Vector<std::complex<double> >  RightEigen = 
      LinearAlgebra::DiagonalizeARPACK(MultFunc(Psi, QShift, q),
                                       n, NumEigen, tol, OutVec, ncv, Sort, Verbose);

   if (LeftVectors)
   {
      tol = tolsave;
      ncv = ncvsave;
      std::vector<std::complex<double> >* OutLeftVec = new std::vector<std::complex<double> >();
      LinearAlgebra::Vector<std::complex<double> >  LeftEigen = 
         LinearAlgebra::DiagonalizeARPACK(MultFuncTrans(Psi, QShift, q),
                                          n, NumEigen, tol, OutLeftVec, ncv, Sort, Verbose);
      // The left vectors are the hermitian conjugate, not the transpose.  So conjugate eigenvalues
      LeftEigen = conj(LeftEigen);
      *LeftVectors = LinearAlgebra::Vector<MatrixOperator>(LeftEigen.size());

      // If we have both left & right eigenvectors, then match the corresponding eigenvalues      
      if (RightVectors)
         LinearAlgebra::MatchEigenvectors(n, LeftEigen, *OutLeftVec, RightEigen, *OutVec, tolsave*10);

      // Unpack the eigenvectors into the output array.  We also conjugate them here
      for (unsigned i = 0; i < LeftEigen.size(); ++i)
      {
         (*LeftVectors)[i] = Pack.unpack(&((*OutLeftVec)[n*i]));
      }
   }

   // eigenvectors
   if (RightVectors)
   {
      *RightVectors = LinearAlgebra::Vector<MatrixOperator>(RightEigen.size());
      for (unsigned i = 0; i < RightEigen.size(); ++i)
      {
         (*RightVectors)[i] = Pack.unpack(&((*OutVec)[n*i]));
      }
   }

   return RightEigen;
}

LinearAlgebra::Vector<std::complex<double> > 
get_spectrum_string(LinearWavefunction const& Psi, QuantumNumber const& QShift, 
                    ProductMPO const& StringOp,
                    int NumEigen,
                    QuantumNumbers::QuantumNumber const& q, double tol = 1e-10,
                    LinearAlgebra::Vector<MatrixOperator>* RightVectors = NULL, 
                    LinearAlgebra::Vector<MatrixOperator>* LeftVectors = NULL, 
                    int ncv = 0, bool Sort = false, int Verbose = 0)
{
   PackMatrixOperator Pack(Psi.Basis2(), Psi.Basis2(), q);
   int n = Pack.size();
   double tolsave = tol;
   int ncvsave = ncv;

   if (Verbose >= 1)
   {
      std::cerr << "Calculating right eigenvalues\n";
   }

   std::vector<std::complex<double> >* OutVec 
      = RightVectors ? new std::vector<std::complex<double> >() : NULL;
   LinearAlgebra::Vector<std::complex<double> >  RightEigen = 
      LinearAlgebra::DiagonalizeARPACK(MultFuncString(Psi, QShift, q, StringOp),
                                       n, NumEigen, tol, OutVec, ncv, Sort, Verbose);

   if (LeftVectors)
   {
      if (Verbose >= 1)
      {
         std::cerr << "Calculating left eigenvalues\n";
      }
      tol = tolsave;
      ncv = ncvsave;
      std::vector<std::complex<double> >* OutLeftVec = new std::vector<std::complex<double> >();
      LinearAlgebra::Vector<std::complex<double> >  LeftEigen = 
         LinearAlgebra::DiagonalizeARPACK(MultFuncStringTrans(Psi, QShift, q, StringOp),
                                          n, NumEigen, tol, OutLeftVec, ncv, Sort, Verbose);
      // The left vectors are the hermitian conjugate, not the transpose.  So conjugate eigenvalues
      LeftEigen = conj(LeftEigen);
      *LeftVectors = LinearAlgebra::Vector<MatrixOperator>(LeftEigen.size());

      // If we have both left & right eigenvectors, then match the corresponding eigenvalues      
      if (RightVectors)
         LinearAlgebra::MatchEigenvectors(n, LeftEigen, *OutLeftVec, RightEigen, *OutVec, tolsave*10);

      // Unpack the eigenvectors into the output array.  We also conjugate them here
      for (unsigned i = 0; i < LeftEigen.size(); ++i)
      {
         (*LeftVectors)[i] = Pack.unpack(&((*OutLeftVec)[n*i]));
      }
   }

   // eigenvectors
   if (RightVectors)
   {
      *RightVectors = LinearAlgebra::Vector<MatrixOperator>(RightEigen.size());
      for (unsigned i = 0; i < RightEigen.size(); ++i)
      {
         (*RightVectors)[i] = Pack.unpack(&((*OutVec)[n*i]));
      }
   }

   return RightEigen;
}

// returns the left/right eigenvector pair corresponding to the eigenvalue 1
// TODO: handle the degenerate case where there is more than one eigenvalue = 1
std::pair<MatrixOperator, MatrixOperator>
get_principal_eigenpair(LinearWavefunction const& Psi, QuantumNumber const& QShift, 
                        double TolIn = 1E-14, int Verbose = 0,
			MatrixOperator LeftGuess = MatrixOperator(),
			MatrixOperator RightGuess = MatrixOperator())
{
   // get the left eigenvector

   LeftMultiply LeftMult(Psi, QShift);
   MatrixOperator LeftEigen = LeftGuess;
   if (LeftEigen.is_null())
      LeftEigen = MakeRandomMatrixOperator(Psi.Basis2(), Psi.Basis2());

   int Iterations = 20;
   double Tol = TolIn;
   LeftEigen = 0.5 * (LeftEigen + adjoint(LeftEigen)); // make the eigenvector symmetric
   std::complex<double> EtaL = LinearSolvers::Arnoldi(LeftEigen, LeftMult, 
                                                      Iterations, Tol, LinearSolvers::LargestAlgebraicReal, false,
						      Verbose);
   while (Iterations == 20)
   {
      Iterations = 20; Tol = 1E-14;
      LeftEigen = 0.5 * (LeftEigen + adjoint(LeftEigen)); // make the eigenvector symmetric
      EtaL = LinearSolvers::Arnoldi(LeftEigen, LeftMult, Iterations, Tol, LinearSolvers::LargestAlgebraicReal, 
				    false, Verbose);
   }

   DEBUG_TRACE(EtaL);

   // get the right eigenvector

   RightMultiply RightMult(Psi, QShift);
   MatrixOperator RightEigen = RightGuess;
   if (RightEigen.is_null())
      RightEigen = MakeRandomMatrixOperator(Psi.Basis2(), Psi.Basis2());

   Iterations = 20; Tol = TolIn;
   RightEigen = 0.5 * (RightEigen + adjoint(RightEigen));
   std::complex<double> EtaR = LinearSolvers::Arnoldi(RightEigen, RightMult, 
                                                      Iterations, Tol, LinearSolvers::LargestAlgebraicReal, 
						      false, Verbose);
   //   DEBUG_TRACE(norm_frob(RightEigen - adjoint(RightEigen)));
   while (Iterations == 20)
   {
      Iterations = 20; Tol = 1E-14;
      RightEigen = 0.5 * (RightEigen + adjoint(RightEigen));
      EtaR = LinearSolvers::Arnoldi(RightEigen, RightMult, Iterations, Tol, LinearSolvers::LargestAlgebraicReal, 
				    false, Verbose);
   }
   DEBUG_TRACE(EtaR);

   RightEigen *= 1.0 / inner_prod(LeftEigen, RightEigen);

   return std::make_pair(LeftEigen, RightEigen);
}

int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      bool ShowRealPart = false, ShowImagPart = false, ShowMagnitude = false;
      bool ShowCartesian = false, ShowPolar = false, ShowArgument = false;
      bool ShowRadians = false, ShowCorrLength = false;
      int UnitCellSize = 0;
      std::string PsiStr;
      std::vector<std::string> Sector;
      double Tol = 1E-15;
      int Iter = 30;
      bool Sort = false;
      bool Quiet = false;
      bool Print = false;
      bool Symmetric = false;
      int KrylovLength = 0;
      std::string String;
      int MaxEigen = 10;
      std::vector<std::string> LeftOpStr;
      std::vector<std::string> RightOpStr;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("cart,c", prog_opt::bool_switch(&ShowCartesian),
	  "show the result in cartesian coordinates [equivalent to --real --imag]")
	 ("polar,p", prog_opt::bool_switch(&ShowPolar),
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
	 ("corr,x", prog_opt::bool_switch(&ShowCorrLength),
	  "display the equivalent correlation length")
	 ("unitcell,u", prog_opt::value(&UnitCellSize),
	  "scale the results to use this unit cell size [default wavefunction unit cell]")
         ("numeigen,n", prog_opt::value(&MaxEigen),
          FormatDefault("Number of eigenvalues to calculate in each sector", MaxEigen).c_str())
	 ("left", prog_opt::value(&LeftOpStr), 
	  "Calculate the expansion coefficients of this operator acting on the left")
	 ("right", prog_opt::value(&RightOpStr), 
	  "Calculate the expansion coefficients of this operator acting on the right")
         ("string", prog_opt::value(&String),
          "use this product operator as a string operator")
         ("quantumnumber,q", prog_opt::value(&Sector),
          "calculate the overlap only in this quantum number sector, "
          "can be used multiple times [default is to calculate all sectors]")
	 ("symmetric", prog_opt::bool_switch(&Symmetric),
	  "transform into the symmetric canonical form")
         ("sort,s", prog_opt::bool_switch(&Sort),
          "sort the eigenvalues by magnitude (not yet implemented)")
         ("tol", prog_opt::value(&Tol),
          FormatDefault("Tolerance of the Arnoldi eigensolver", Tol).c_str())
         //("iter", prog_opt::value(&Iter),
	 //FormatDefault("Maximum subspace size in the Arnoldi basis", Iter).c_str())
	 ("krylov,k", prog_opt::value(&KrylovLength), 
	  "Length of the Krylov sequence [default 2*num-eigenvalues]")
         ("quiet", prog_opt::bool_switch(&Quiet),
          "don't show the column headings")
         ("print", prog_opt::bool_switch(&Print), "with --string, Print the MPO to standard output")
	 //         ("overlaps", prog_opt::bool_switch(&ShowEigenvectorOverlaps),
	 //"Write the matrix of overlaps of the left/right eigenvectors to cerr (for debugging)")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose),
          "extra debug output [can be used multiple times]")
	 ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi", prog_opt::value(&PsiStr), "psi1")
         ;

      prog_opt::positional_options_description p;
      p.add("psi", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);    

      if (vm.count("help") > 0 || vm.count("psi") < 1)
      {
         print_copyright(std::cerr);
         std::cerr << "usage: mp-ispectrum [options] <psi>\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      // If no output switches are used, default to showing everything
      if (!ShowRealPart && !ShowImagPart && !ShowMagnitude
	  && !ShowCartesian && !ShowPolar && !ShowArgument
	  && !ShowCorrLength)
      {
	 ShowCartesian = true;
	 ShowPolar = true;
	 ShowCorrLength = true;
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

      SimpleOperator MyOpL, MyOpR, MyOpL2;
      LatticeSite Site;

      if (Verbose)
         std::cout << "Loading wavefunction...\n";

      pvalue_ptr<MPWavefunction> Psi1 
         = pheap::OpenPersistent(PsiStr, mp_pheap::CacheSize(), true);

      // firstly, get the LinearWavefunction
      InfiniteWavefunctionLeft InfPsi = Psi1->get<InfiniteWavefunctionLeft>();

      LinearWavefunction Psi;
      QuantumNumber QShift = InfPsi.qshift();

      RealDiagonalOperator D;
      boost::tie(Psi, D) = get_left_canonical(InfPsi);

      if (Symmetric)
      {
         MatrixOperator LambdaSqrt = D;
	 LambdaSqrt = SqrtDiagonal(LambdaSqrt);
         MatrixOperator LambdaInvSqrt = InvertDiagonal(LambdaSqrt, InverseTol);
         Psi.set_back(prod(Psi.get_back(), LambdaSqrt));
         Psi.set_front(prod(delta_shift(LambdaInvSqrt, QShift), Psi.get_front()));
      }


      // now get the principal eigenpair
      MatrixOperator LeftIdent, RightIdent;
      if (Symmetric)
      {
         MatrixOperator R = D;

	 if (Verbose)
	    std::cout << "Solving principal eigenpair...\n";
	 boost::tie(LeftIdent, RightIdent) = get_principal_eigenpair(Psi, QShift, Tol, Verbose,
								     R, R);
      }
      else
      {
	 // if we're in the left canonical basis, then we know what the eigenpair is.
	 if (Verbose)
	    std::cout << "In left-canonical basis: principal eigenpair is already known.\n";
	 LeftIdent = MatrixOperator::make_identity(Psi.Basis2());
	 RightIdent = D;
	 RightIdent = scalar_prod(RightIdent, herm(RightIdent));

	 boost::tie(LeftIdent, RightIdent) = get_principal_eigenpair(Psi, QShift, Tol, Verbose,
								     LeftIdent, RightIdent);
      }

      std::complex<double> IdentNormalizationFactor = inner_prod(LeftIdent, RightIdent);

      // Get the string operator
      ProductMPO StringOp;
      if (vm.count("string"))
      {
	 InfiniteLattice Lattice;
	 boost::tie(StringOp, Lattice) = ParseProductOperatorAndLattice(String);
	 if (Print)
	 {
	    std::cout << "String MPO is:\n" << StringOp << '\n';
	 }
	 CHECK(Psi.size() % StringOp.size() == 0)
	    ("Wavefunction size must be a multiple of the string operator size")
	    (Psi.size())(StringOp.size());
	 StringOp = repeat(StringOp, Psi.size() / StringOp.size());
      }
      else
      {
         StringOp = ProductMPO::make_identity(ExtractLocalBasis(Psi));
      }

      // The default UnitCellSize for output is the wavefunction size
      if (UnitCellSize == 0)
      {
	 UnitCellSize = Psi.size();
      }
      double ScaleFactor = double(UnitCellSize) / double(Psi.size());

      // get the set of quantum numbers to show
      typedef std::set<QuantumNumbers::QuantumNumber> QSetType;
      QSetType QL;
      if (vm.count("quantumnumber") != 0)
      {
	 while (!Sector.empty())
	 {
	    QL.insert(QuantumNumbers::QuantumNumber(Psi.GetSymmetryList(), Sector.back()));
	    Sector.pop_back();
	 }
      }
      else
      {
         // Assemble the list of all possible quantum numbers
         QSetType Q1 = QuantumNumbersInBasis(Psi.Basis2());
         for (QSetType::const_iterator I = Q1.begin(); I != Q1.end(); ++I)
         {
            for (QSetType::const_iterator J = Q1.begin(); J != Q1.end(); ++J)
            {
               QuantumNumberList NextQ = inverse_transform_targets(*I, *J);
               QL.insert(NextQ.begin(), NextQ.end());
            }
         }
      }

      // show the title
      if (!Quiet)
      {
	 std::cout << "#" << argv[0];
	 for (int i = 1; i < argc; ++i)
	    std::cout << ' ' << argv[i];
	 std::cout << "\n#quantities are calculated per unit cell size of " << UnitCellSize 
		   << (UnitCellSize == 1 ? " site\n" : " sites\n");
         std::cout << "#sector     #n   ";
         if (ShowRealPart)
            std::cout << "#real                 ";
         if (ShowImagPart)
            std::cout << "#imag                 ";
         if (ShowCorrLength)
            std::cout << "#corr_length          ";
         if (ShowMagnitude)
            std::cout << "#magnitude            ";
         if (ShowArgument)
            std::cout << "#argument" << (ShowRadians ? "(rad)" : "(deg)") << "        ";

	 // titles for the overlaps
	 for (unsigned i = 0; i < LeftOpStr.size(); ++i)
	 {
	    for (unsigned j = 0; j < RightOpStr.size(); ++j)
	    {
	       std::string Title = "#overlap_" + boost::lexical_cast<std::string>(i) + '_'
		  + boost::lexical_cast<std::string>(j); 
	       if (ShowRealPart)
		  std::cout << std::setw(20) << std::left << (Title+"_real") << "  ";
	       if (ShowImagPart)
		  std::cout << std::setw(20) << std::left << (Title+"_imag") << "  ";
	       if (ShowMagnitude)
		  std::cout << std::setw(20) << std::left << (Title+"_mag") << "  ";
	       if (ShowArgument)
		  std::cout << std::setw(20) << std::left << (Title+"_arg" + (ShowRadians ? "(rad)" : "(deg)")) << "  ";
	    }
	 }
         std::cout << '\n';
         std::cout << std::left;
      }

      // initialize the finite operators for the observables
      std::vector<MatrixOperator> LeftOp;
      std::vector<MatrixOperator> RightOp;

      //LeftIdent = delta_shift(LeftIdent, QShift);

      for (unsigned i = 0; i < LeftOpStr.size(); ++i)
      {
	 UnitCellMPO Op = ParseUnitCellOperatorAndLattice(LeftOpStr[i]).first;
         Op.ExtendToCoverUnitCell(Psi.size());
         // Adjust the MPO so that the left basis is the identity
         FiniteMPO Mpo = Op.MPO() * FiniteMPO::make_identity(Op.MPO().LocalBasis2List(), 
                                                              adjoint(Op.TransformsAs()));
         Mpo = project(Mpo, QuantumNumbers::QuantumNumber(Op.GetSymmetryList()));
	 LeftOp.push_back(inject_left(LeftIdent, Psi, QShift, Mpo, Psi));
      }

      for (unsigned i = 0; i < RightOpStr.size(); ++i)
      {
	 UnitCellMPO Op = ParseUnitCellOperatorAndLattice(RightOpStr[i]).first;
         Op.ExtendToCoverUnitCell(Psi.size());
	 RightOp.push_back(inject_right(RightIdent, Psi, QShift, Op.MPO(), Psi));
      }

      // iterate over the relevant quantum number sectors
      for (QSetType::const_iterator qI = QL.begin(); qI != QL.end(); ++qI)
      {
         LinearAlgebra::Vector<MatrixOperator> RightEigenvectors;
         LinearAlgebra::Vector<MatrixOperator> LeftEigenvectors;
         LinearAlgebra::Vector<std::complex<double> > EValues;

         LinearAlgebra::Vector<MatrixOperator>* RightEigenvectorsPtr = NULL;
         LinearAlgebra::Vector<MatrixOperator>* LeftEigenvectorsPtr = NULL;

	 if (!RightOp.empty())
         {
            RightEigenvectorsPtr = &RightEigenvectors;
            LeftEigenvectorsPtr = &LeftEigenvectors;
         }

         // determine the spectrum
         EValues = get_spectrum_string(Psi, QShift, StringOp * ProductMPO::make_identity(StringOp.LocalBasis2List(), *qI),
				       MaxEigen, *qI, Tol, RightEigenvectorsPtr, 
                                       LeftEigenvectorsPtr, KrylovLength, true, Verbose);

#if 0
         if (!RightOp.is_null() && RightOp.TransformsAs() == *qI && ShowEigenvectorOverlaps)
         {
            LinearAlgebra::Matrix<double> EigenOverlaps(size(RightEigenvectors), size(LeftEigenvectors));
            for (unsigned i = 0; i < size(RightEigenvectors); ++i)
               for (unsigned j = 0; j < size(LeftEigenvectors); ++j)
                  EigenOverlaps(i,j) = norm_frob(inner_prod(LeftEigenvectors[j], RightEigenvectors[i]));
            std::cerr << "Eigenvector overlap matrix:\n" << EigenOverlaps << '\n';
         }
#endif

         for (int i = 0; i < int(size(EValues)); ++i)
         {
            PrintFormat(*qI, EValues[i], i, ShowRealPart, ShowImagPart, ShowCorrLength, ShowMagnitude,
			ShowArgument, ShowRadians, ScaleFactor);

            // show eigenvector info?
	    for (unsigned iL = 0; iL < LeftOp.size(); ++iL)
	    {
	       for (unsigned iR = 0; iR < RightOp.size(); ++iR)
	       {
		  if (LeftOp[iL].TransformsAs() == LeftEigenvectors[i].TransformsAs()
		      && RightOp[iR].TransformsAs() == RightEigenvectors[i].TransformsAs())
		  {
#if 1
		     std::complex<double> Overlap = inner_prod(LeftOp[iL], LeftEigenvectors[i])
			* inner_prod(RightEigenvectors[i], RightOp[iR])
			/ (inner_prod(RightEigenvectors[i], LeftEigenvectors[i]) * IdentNormalizationFactor);
		     //TRACE(inner_prod(LeftOp[iL], LeftEigenvectors[i]))(inner_prod(RightEigenvectors[i], RightOp[iR]))
		     //(inner_prod(RightEigenvectors[i], LeftEigenvectors[i]));
#else
		     std::complex<double> Overlap = inner_prod(LeftOp[iL], RightEigenvectors[i])
			* inner_prod(LeftEigenvectors[i], RightOp[iR])
			/ (inner_prod(LeftEigenvectors[i], RightEigenvectors[i]) * IdentNormalizationFactor);
#endif
		     if (ShowRealPart)
			std::cout << std::setw(20) << Overlap.real() << "  ";
		     if (ShowImagPart)
			std::cout << std::setw(20) << Overlap.imag() << "  ";
		     double Magnitude = norm_frob(Overlap);
		     double Arg = std::atan2(Overlap.imag(), Overlap.real());
		     if (!ShowRadians)
			Arg *= 180.0 / math_const::pi;
		     if (ShowMagnitude)
			std::cout << std::setw(20) << Magnitude << "  ";
		     if (ShowArgument)
			std::cout << std::setw(20) << Arg << "  ";
		  }
	       }
	    }
            std::cout << '\n';

         }// for i

      } // for qI

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
