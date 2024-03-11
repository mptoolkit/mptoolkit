// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/old/mp-ispectrum2.cpp
//
// Copyright (C) 2012-2017 Ian McCulloch <ian@qusim.net>
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

//
// Version of mp-ispectrum where we use the UnitCell structure
//

#include "mps/infinitewavefunction.h"
#include "mps/packunpack.h"
#include "lattice/latticesite.h"
#include "mps/operator_actions.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include "common/prog_opt_accum.h"
#include "common/environment.h"
#include "interface/inittemp.h"
#include "tensor/tensor_eigen.h"
#include "mp-algorithms/arnoldi.h"
#include "lattice/unitcell-parser.h"

#include "models/spin.h"
#include "models/spin-u1.h"
#include "models/spin-u1u1.h"
#include "models/spin-z2.h"
#include "models/spin-su2.h"
#include "models/tj-u1su2.h"
#include "models/tj-u1.h"
#include "models/spinlessfermion-u1.h"
#include "models/kondo-u1su2.h"
#include "models/kondo-u1.h"
#include "models/kondo-u1u1.h"
#include "models/kondo-so4.h"
#include "models/boson-u1.h"
#include "models/boson-2component-u1z2.h"
#include "models/hubbard-u1u1-old.h"
#include "models/hubbard-u1u1.h"
#include "models/hubbard-u1su2.h"
#include "models/hubbard-so4.h"

#include "linearalgebra/arpack_wrapper.h"

//
// TODO: the decomposition into left/right operators is bogus, since the magnitudes of the left/right operators
// is also unspecified.  We should get rid of that, and add one (pair) column for each left/right pair that
// makes physical sense (ie the quantum numbers match)
//

namespace prog_opt = boost::program_options;

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

struct LeftMultiplyString
{
   typedef MatrixOperator argument_type;
   typedef MatrixOperator result_type;

   LeftMultiplyString(LinearWavefunction const& L_, QuantumNumber const& QShift_,
                      GenericMPO const& StringOp_)
      : L(L_), QShift(QShift_), StringOp(StringOp_)
   {
      CHECK_EQUAL(L.size(), StringOp.size())("The string operator must be the same length as the unit cell");
   }

   result_type operator()(argument_type const& x) const
   {
      result_type r = delta_shift(x, QShift);
      r = inject_left(r, L, StringOp, L);
      return r;
   }

   LinearWavefunction const& L;
   QuantumNumber QShift;
   GenericMPO StringOp;
};

struct LeftMultiplyOp
{
   typedef MatrixOperator argument_type;
   typedef MatrixOperator result_type;

   LeftMultiplyOp(LinearWavefunction const& L_, QuantumNumber const& QShift_, SimpleOperator const& Op_)
      : L(L_), QShift(QShift_), Op(Op_) {}

   result_type operator()(argument_type const& x) const
   {
      result_type r = delta_shift(x, QShift);
      for (LinearWavefunction::const_iterator I = L.begin(); I != L.end(); ++I)
      {
         r = operator_prod(herm(Op), herm(*I), r, *I);
      }
      return r;
   }

   LinearWavefunction const& L;
   QuantumNumber QShift;
   SimpleOperator const& Op;
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
                       GenericMPO const& StringOp_)
      : L(L_), QShift(QShift_) , StringOp(StringOp_)
   {
      CHECK_EQUAL(L.size(), StringOp.size())("The string operator must be the same length as the unit cell");
   }

   result_type operator()(argument_type const& x) const
   {
      result_type r = inject_right(x, L, StringOp, L);
      return delta_shift(r, adjoint(QShift));
   }

   LinearWavefunction const& L;
   QuantumNumber QShift;
   GenericMPO StringOp;
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
                  BasicFiniteMPO const& StringOp)
      : Mult(Psi, QShift, StringOp), Pack(Psi.Basis2(), Psi.Basis2(), q) { }

   void operator()(std::complex<double> const* In, std::complex<double>* Out) const
   {
      MatrixOperator x = Pack.unpack(In);
      x = Mult(x);
      Pack.pack(x, Out);
   }

   LeftMultiplyString Mult;
   PackMatrixOperator Pack;
};

struct MultFuncStringTrans
{
   MultFuncStringTrans(LinearWavefunction const& Psi, QuantumNumber const& QShift,
                       QuantumNumbers::QuantumNumber const& q,
                       GenericMPO const& StringOp)
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
                    BasicFiniteMPO const& StringOp,
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
                        double TolIn = 1E-14, int Verbose = 0)
{
   // get the left eigenvector

   LeftMultiply LeftMult(Psi, QShift);
   MatrixOperator LeftEigen = MakeRandomMatrixOperator(Psi.Basis2(), Psi.Basis2());

   int Iterations = 20;
   double Tol = TolIn;
   LeftEigen = 0.5 * (LeftEigen + adjoint(LeftEigen)); // make the eigenvector symmetric
   std::complex<double> EtaL = LinearSolvers::Arnoldi(LeftEigen, LeftMult,
                                                      Iterations, Tol, LinearSolvers::LargestAlgebraicReal, false);
   while (Iterations == 20)
   {
      Iterations = 20; Tol = 1E-14;
      LeftEigen = 0.5 * (LeftEigen + adjoint(LeftEigen)); // make the eigenvector symmetric
      EtaL = LinearSolvers::Arnoldi(LeftEigen, LeftMult, Iterations, Tol, LinearSolvers::LargestAlgebraicReal, false);
   }

   DEBUG_TRACE(EtaL);

   // get the right eigenvector

   RightMultiply RightMult(Psi, QShift);
   MatrixOperator RightEigen = MakeRandomMatrixOperator(Psi.Basis2(), Psi.Basis2());

   Iterations = 20; Tol = TolIn;
   RightEigen = 0.5 * (RightEigen + adjoint(RightEigen));
   std::complex<double> EtaR = LinearSolvers::Arnoldi(RightEigen, RightMult,
                                                      Iterations, Tol, LinearSolvers::LargestAlgebraicReal, false);
   //   DEBUG_TRACE(norm_frob(RightEigen - adjoint(RightEigen)));
   while (Iterations == 20)
   {
      Iterations = 20; Tol = 1E-14;
      RightEigen = 0.5 * (RightEigen + adjoint(RightEigen));
      EtaR = LinearSolvers::Arnoldi(RightEigen, RightMult, Iterations, Tol, LinearSolvers::LargestAlgebraicReal, false);
   }
   DEBUG_TRACE(EtaR);

   RightEigen *= 1.0 / inner_prod(LeftEigen, RightEigen);

   return std::make_pair(LeftEigen, RightEigen);
}



void PrintEigenvalue(std::ostream& out, std::complex<double> x, bool ShowCart, bool ShowRadians)
{
   if (ShowCart)
   {
      out << std::setw(21) << std::scientific << x.real() << "  "
          << std::setw(21) << std::scientific << x.imag() << "  ";
   }
   else
   {
      double Magnitude = norm_frob(x);
      double Arg = std::atan2(x.imag(), x.real());
      if (!ShowRadians)
         Arg *= 180.0 / math_const::pi;
      out << std::setw(21) << std::scientific << Magnitude << "  "
          << std::setw(21) << std::fixed << Arg << "  ";
   }
}

void PrintEmptyEigenvalue(std::ostream& out, bool ShowCart, bool ShowRadians)
{
   out << std::setw(21) << ' ' << "  "
       << std::setw(21) << ' ' << "  ";
}


int main(int argc, char** argv)
{
   try
   {
      //bool Verbose = false;
      std::string PsiStr;
      //      double EigenCutoff = 0.1;
      int MaxEigen = 10;
      std::string Target;
      std::vector<std::string> OpL, OpR;
      double Tol = 1e-10;
      int Verbose = 0;
      int KrylovLength = 0;
      bool Symmetric = false;
      std::string Model;
      double Spin = 0.5;
      int NMax = 3;
      bool ShowEigenvaluesCart = false;      // show the eigenvalues in (real,imag) units
      bool ShowRadians = false;              // show the angle in radians instead of degrees
      bool ShowCorrelationLength = false;
      int UnitCellSize = 1;
      bool ShowEigenvectorOverlaps = false;
      bool Debug = false;
      std::string StringOp = "I";

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         //         ("eigen-cutoff,d", prog_opt::value(&EigenCutoff),
         //          ("Cutoff threshold for eigenvalues [default "
         //           +boost::lexical_cast<std::string>(EigenCutoff)+"]").c_str())
         ("num-eigenvalues,n", prog_opt::value(&MaxEigen),
          FormatDefault("Number of eigenvalues to calculate", MaxEigen).c_str())
         ("quantumnumber,q", prog_opt::value(&Target), "Calculate spectrum only in this symmetry sector")
         ("tol,t", prog_opt::value(&Tol), FormatDefault("Tolerance of eigensolver", Tol).c_str())
         ("model", prog_opt::value(&Model), "use this model for the operator to examine "
          "(tri-u1)")
         ("spin", prog_opt::value(&Spin), "for spin models, the value of the spin [default 0.5]")
         ("nmax", prog_opt::value(&NMax), "for Bose-Hubbard model, the max number of bosons per site")
         ("lcoefficients,l", prog_opt::value(&OpL),
          "Calculate the expansion coefficients of this operator acting on the left (can be used more than once)")
         ("rcoefficients,r", prog_opt::value(&OpR),
          "Calculate the expansion coefficients of this operator acting on the right (can be used more than once)")
         ("symmetric,s", prog_opt::bool_switch(&Symmetric),
          "Transform the state to the approximate symmetric orthonormalization constraint")
         ("krylov,k", prog_opt::value(&KrylovLength),
          "Length of the Krylov sequence [default 2*num-eigenvalues]")
         ("cart", prog_opt::bool_switch(&ShowEigenvaluesCart),
          "Show the real and imaginary parts of the eigenvalue, instead of (magnitude,argument)")
         ("radians", prog_opt::bool_switch(&ShowRadians),
          "Show the argument of the eigenvalues in radians, instead of degrees")
         ("corr", prog_opt::bool_switch(&ShowCorrelationLength),
          "Show the correlation length of the eigenvalue")
         ("string", prog_opt::value(&StringOp),
          FormatDefault("String operator", StringOp).c_str())
         ("overlaps", prog_opt::bool_switch(&ShowEigenvectorOverlaps),
          "Write the matrix of overlaps of the left/right eigenvectors to cerr (for debugging)")
         ("verbose,v", prog_opt_ext::accum_value(&Verbose), "increase verbosity")
         ("debug", prog_opt::bool_switch(&Debug), "debug mode")
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
         std::cerr << "usage: mp-ispectrum [options] <psi>\n";
         std::cerr << desc << '\n';
         return 1;
      }

      // if the number of eigenvalues is specified but
      // the cutoff is not, then set the cutoff to zero
      //      if (vm.count("num-eigenvalues") == 1 && vm.count("eigen-cutoff") == 0)
      //      {
      //         EigenCutoff = 0;
      //      }

      UnitCell Cell;
      if (Model == "tri-u1")
      {
         LatticeSite Site = CreateU1U1HubbardSite();
         Cell = UnitCell(3, Site);
      }
      else if (Model != "")
      {
         PANIC("Unknown model");
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      pvalue_ptr<InfiniteWavefunction> Psi1
         = pheap::OpenPersistent(PsiStr, mp_pheap::CacheSize(), true);

      // firstly, get the LinearWavefunction
      QuantumNumber QShift = Psi1->QShift;
      MatrixOperator x_unit =  Psi1->C_right * delta_shift(InvertDiagonal(Psi1->C_old, InverseTol),
                                                           adjoint(QShift));
      LinearWavefunction Psi = Psi1->Psi;
      Psi.set_back(prod(Psi.get_back(), x_unit));
      if (Symmetric)
      {
         MatrixOperator LambdaSqrt = SqrtDiagonal(Psi1->C_old);
         MatrixOperator LambdaInvSqrt = InvertDiagonal(LambdaSqrt, InverseTol);
         Psi.set_back(prod(Psi.get_back(), delta_shift(LambdaSqrt, adjoint(QShift))));
         Psi.set_front(prod(LambdaInvSqrt, Psi.get_front()));
      }

      if (Debug)
         std::cerr << Psi.get_front() << '\n';

      // now get the principal eigenpair
      MatrixOperator LeftIdent, RightIdent;
      std::tie(LeftIdent, RightIdent) = get_principal_eigenpair(Psi, QShift);
      std::complex<double> IdentNormalizationFactor = inner_prod(LeftIdent, RightIdent);

      // Assemble the string operator (which may be the identity)
      BasicFiniteMPO MyStringOp = ParseUnitCellOperator(Cell, StringOp);
      MyStringOp = repeat(MyStringOp, Psi.size()/MyStringOp.size());

      // get the set of quantum numbers to show
      typedef std::set<QuantumNumbers::QuantumNumber> QSetType;
      QSetType QL;
      if (vm.count("quantumnumber") != 0)
      {
         QL.insert(QuantumNumbers::QuantumNumber(Psi.GetSymmetryList(), Target));
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

      // the wavefunction size should be a multiple of the primitive unit cell size,
      // but there is no bad consequences here if it isn't, except that the correlation length
      // will be a strange number
      if (ShowCorrelationLength && (Psi.size() % UnitCellSize != 0))
      {
         std::cerr << "mp-ispectrum: warning: the wavefunction unit cell is "
            "not a multiple of the primitive unit cell!\n";
      }

      // show the help info
      if (Verbose >= 0)
      {
         if (ShowCorrelationLength)
            std::cout << "#Unit cell for the correlation lengths is " << UnitCellSize << "sites.\n";
         if (!OpL.empty() || !OpR.empty())
         {
            std::cout << "#Absolute phases for left and right observables is arbitrary.\n"
                      << "#Relative phases of operators between the same eigenvalue are well-defined.\n"
                      << "#The product of a left operator with a right operator has a well-defined phase.\n";
         }
      }

      // show the title
      if (ShowEigenvaluesCart)
         std::cout << "#sector      #n  #e-value_real          #e-value_imag          ";
      else
         std::cout << "#sector      #n   #e-value_mag            #e-value_arg         ";
      if (ShowCorrelationLength)
         std::cout << " #xi                   ";

      if (!OpL.empty())
         std::cout << "#L ";
      for (unsigned j = 0; j < OpL.size(); ++j)
      {
         std::cout << "#" << std::left << std::setw(20)
                   << (OpL[j] + (ShowEigenvaluesCart ? " (real)" : " (mag)")) << "  ";
         std::cout << "#" << std::left << std::setw(20)
                   << (OpL[j] + (ShowEigenvaluesCart ? " (imag)" : " (arg)")) << "  ";
      }
      if (!OpR.empty())
         std::cout << "#R ";
      for (unsigned j = 0; j < OpR.size(); ++j)
      {
         std::cout << "#" << std::left << std::setw(20)
                   << (OpR[j] + (ShowEigenvaluesCart ? " (real)" : " (mag)")) << "  ";
         std::cout << "#" << std::left << std::setw(20)
                   << (OpR[j] + (ShowEigenvaluesCart ? " (imag)" : " (arg)")) << "  ";
      }
      std::cout << '\n' << std::right;

      // initialize the linear operators for the observables
      std::vector<MatrixOperator> LeftObservables, RightObservables;

      for (unsigned i = 0; i < OpL.size(); ++i)
      {
         GenericMPO Op = ParseUnitCellOperator(Cell, OpL[i]);
         local_basis_compatible_or_abort(Psi, Op);
         LeftObservables.push_back(inject_left_qshift(LeftIdent, Op, Psi, QShift));
      }
      for (unsigned i = 0; i < OpR.size(); ++i)
      {
         GenericMPO Op = ParseUnitCellOperator(Cell, OpR[i]);
         local_basis_compatible_or_abort(Psi, Op);
         RightObservables.push_back(inject_right_qshift(RightIdent, Op, Psi, QShift));
      }

      // iterate over the relevant quantum number sectors
      for (QSetType::const_iterator qI = QL.begin(); qI != QL.end(); ++qI)
      {
         LinearAlgebra::Vector<MatrixOperator> RightEigenvectors;
         LinearAlgebra::Vector<MatrixOperator> LeftEigenvectors;
         LinearAlgebra::Vector<std::complex<double> > EValues;

         // determine the spectrum
         EValues = get_spectrum_string(Psi, QShift, MyStringOp*identity_mpo(Cell, *qI),
                                       MaxEigen, *qI, Tol, &RightEigenvectors,
                                       &LeftEigenvectors, KrylovLength, true, Verbose);

         if (ShowEigenvectorOverlaps)
         {
            LinearAlgebra::Matrix<double> EigenOverlaps(size(RightEigenvectors), size(LeftEigenvectors));
            for (unsigned i = 0; i < size(RightEigenvectors); ++i)
               for (unsigned j = 0; j < size(LeftEigenvectors); ++j)
                  EigenOverlaps(i,j) = norm_frob(inner_prod(LeftEigenvectors[j], RightEigenvectors[i]));
            std::cerr << "Eigenvector overlap matrix:\n" << EigenOverlaps << '\n';
         }

         for (int i = 0; i < int(size(EValues)); ++i)
         {
            std::map<int, std::complex<double> > LeftOverlaps;
            std::map<int, std::complex<double> > RightOverlaps;

            std::cout << std::setw(7) << boost::lexical_cast<std::string>(*qI) << "  "
                      << std::setw(6) << i << "  ";

            PrintEigenvalue(std::cout, EValues[i], ShowEigenvaluesCart, ShowRadians);

            if (ShowCorrelationLength)
            {
               double CorrLen = -1.0 / std::log(norm_frob(EValues[i]));
               CorrLen *= double(Psi.size()) / UnitCellSize;
               std::cout << std::setw(21) << std::scientific << CorrLen << "  ";
            }

            // The phase of each component is arbitrary, only the relative phase counts.
            // Also the phase of a left observable multplied by a right observable is well-defined.
            if (!LeftObservables.empty())
               std::cout << "   "; // for the "#L " title
            for (unsigned j = 0; j < LeftObservables.size(); ++j)
            {
               if (LeftObservables[j].TransformsAs() == *qI)
               {
                  std::complex<double> E = inner_prod(LeftObservables[j], LeftEigenvectors[i]);
                  PrintEigenvalue(std::cout, E, ShowEigenvaluesCart, ShowRadians);
               }
               else
               {
                  PrintEmptyEigenvalue(std::cout, ShowEigenvaluesCart, ShowRadians);
               }
            }

            if (!RightObservables.empty())
               std::cout << "   "; // for the "#R " title
            for (unsigned j = 0; j < RightObservables.size(); ++j)
            {
               if (RightObservables[j].TransformsAs() == *qI)
               {
                  std::complex<double> E = inner_prod(RightEigenvectors[i], RightObservables[j])
                     / (inner_prod(RightEigenvectors[i], LeftEigenvectors[i])*IdentNormalizationFactor);
                  PrintEigenvalue(std::cout, E, ShowEigenvaluesCart, ShowRadians);
               }
               else
               {
                  PrintEmptyEigenvalue(std::cout, ShowEigenvaluesCart, ShowRadians);
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
