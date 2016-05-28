// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-ispectrum-old.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "mps/infinitewavefunction.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include "common/prog_opt_accum.h"
#include "common/environment.h"
#include "interface/inittemp.h"
#include "tensor/tensor_eigen.h"

#include "models/spin.h"
#include "models/spin-u1.h"
#include "models/spin-z2.h"
#include "models/spin-su2.h"
#include "models/tj-u1su2.h"
#include "models/tj-u1.h"
#include "models/spinlessfermion-u1.h"
#include "models/kondo-u1su2.h"
#include "models/kondo-u1.h"
#include "models/bosehubbard-spinless-u1.h"
#include "models/hubbard-u1u1.h"

#include "mps/packunpack.h"
#include "linearalgebra/arpack_wrapper.h"

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
      return delta_shift(x, QShift);
   }

   LinearWavefunction const& L;
   QuantumNumber QShift;
   SimpleOperator const& Op;
};

struct MultFunc
{
   MultFunc(LinearWavefunction const& Psi, QuantumNumber const& QShift, 
            QuantumNumbers::QuantumNumber const& q, SimpleOperator const& UnitOp)
      : Mult(Psi, QShift, UnitOp), Pack(Psi.Basis2(), Psi.Basis2(), q) {}

   void operator()(std::complex<double> const* In, std::complex<double>* Out) const
   {
      MatrixOperator x = Pack.unpack(In);
      x = Mult(x);
      Pack.pack(x, Out);
   }

   LeftMultiplyOp Mult;
   PackMatrixOperator Pack;
};

struct MultFuncTrans
{
   MultFuncTrans(LinearWavefunction const& Psi, QuantumNumber const& QShift, 
            QuantumNumbers::QuantumNumber const& q, SimpleOperator const& UnitOp)
      : Mult(Psi, QShift, UnitOp), Pack(Psi.Basis2(), Psi.Basis2(), q) {}

   void operator()(std::complex<double> const* In, std::complex<double>* Out) const
   {
      MatrixOperator x = Pack.unpack(In);
      x = Mult(x);
      Pack.pack(x, Out);
   }

   RightMultiplyOp Mult;
   PackMatrixOperator Pack;
};

LinearAlgebra::Vector<std::complex<double> > 
get_spectrum(LinearWavefunction const& Psi, QuantumNumber const& QShift, int NumEigen,
             QuantumNumbers::QuantumNumber const& q, SimpleOperator const& UnitOp, double tol = 1e-10,
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
      LinearAlgebra::DiagonalizeARPACK(MultFunc(Psi, QShift, q, UnitOp),
                                       n, NumEigen, tol, OutVec, ncv, Sort, Verbose);

   int nev = RightEigen.size();

   if (LeftVectors)
   {
      tol = tolsave;
      ncv = ncvsave;
      std::vector<std::complex<double> >* OutLeftVec 
         = RightVectors ? new std::vector<std::complex<double> >() : NULL;
      LinearAlgebra::Vector<std::complex<double> >  LeftEigen = 
         LinearAlgebra::DiagonalizeARPACK(MultFuncTrans(Psi, QShift, q, UnitOp),
                                          n, NumEigen, tol, OutLeftVec, ncv, Sort, Verbose);
      
      if (RightVectors)
         LinearAlgebra::MatchEigenvectors(n, LeftEigen, *OutLeftVec, RightEigen, *OutVec);

      for (unsigned i = 0; i < LeftEigen.size(); ++i)
      {
         (*LeftVectors)[i] = Pack.unpack(&((*OutLeftVec)[n*i]));
      }
   }

   // eigenvectors
   if (RightVectors)
   {
      *RightVectors = LinearAlgebra::Vector<MatrixOperator>(nev);
      for (unsigned i = 0; i < RightEigen.size(); ++i)
      {
         (*RightVectors)[i] = Pack.unpack(&((*OutVec)[n*i]));
      }
   }

   return RightEigen;
}

LinearAlgebra::Vector<std::complex<double> > 
make_spectrum(InfiniteWavefunction const& x, QuantumNumbers::QuantumNumber const& q,
              int MaxEigen, SimpleOperator const& UnitOp, double Tol, int ncv = 0, bool Symmetric = false,
              LinearAlgebra::Vector<MatrixOperator>* RightVectors = NULL,
              LinearAlgebra::Vector<MatrixOperator>* LeftVectors = NULL,
              int Verbose = 0)
{
   MatrixOperator x_unit =  x.C_right * delta_shift(InvertDiagonal(x.C_old, InverseTol), adjoint(x.QShift));
   LinearWavefunction xPsi = x.Psi;
   xPsi.set_back(prod(xPsi.get_back(), x_unit));
   if (Symmetric)
   {
      MatrixOperator LambdaSqrt = SqrtDiagonal(x.C_old);
      MatrixOperator LambdaInvSqrt = InvertDiagonal(LambdaSqrt, InverseTol);
      xPsi.set_back(prod(xPsi.get_back(), delta_shift(LambdaSqrt, adjoint(x.QShift))));
      xPsi.set_front(prod(LambdaInvSqrt, xPsi.get_front()));
   }
   return get_spectrum(xPsi, x.QShift, MaxEigen, q, UnitOp, Tol, RightVectors, LeftVectors, ncv, true, Verbose);
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
      std::string OpL, OpR;
      double Tol = 1e-10;
      int Verbose = 0;
      int KrylovLength = 0;
      bool Symmetric = false;
      std::string Model;
      std::string OpL2 = "I";
      bool NormOnly = false;
      bool ShowMagnitude = false;
      double Spin = 0.5;
      int NMax = 3;
      std::string OpStr = "I";

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
	  "(spin, spin-su2, spin-u1, tj-u1, sf-u1, spin-z2)")
	 ("spin", prog_opt::value(&Spin), "for spin models, the value of the spin [default 0.5]")
         ("nmax", prog_opt::value(&NMax), "for Bose-Hubbard model, the max number of bosons per site")
	 ("string", prog_opt::value(&OpStr), "use this unitary operator for the transfer operator")
	 ("lcoefficients,l", prog_opt::value(&OpL), 
	  "Calculate the expansion coefficients of this operator acting on the left")
	 ("l2", prog_opt::value(&OpL2),
	  "Hack to do a 2-site operator")
	 ("rcoefficients,r", prog_opt::value(&OpR), 
	  "Calculate the expansion coefficients of this operator acting on the right")
	 ("symmetric,s", prog_opt::bool_switch(&Symmetric),
	  "Transform the state to the approximate symmetric orthonormalization constraint")
	 ("krylov,k", prog_opt::value(&KrylovLength), 
	  "Length of the Krylov sequence [default 2*num-eigenvalues]")
	 ("mag", prog_opt::bool_switch(&ShowMagnitude),
	  "Show the magnitude of the eigenvalue, instead of the real and maginary parts of the eigenvalue")
	 ("verbose,v", prog_opt_ext::accum_value(&Verbose), "increase verbosity")
	 ("norm", prog_opt::bool_switch(&NormOnly), "show the norm of the overlap, instead of the"
	  " real and imag parts")
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

      SimpleOperator MyOpL, MyOpR, MyOpL2;
      SiteBlock Site;
      if (Model == "spin")
      {
	 Site = CreateSpinSite(Spin);
      }
      else if (Model == "spin-su2")
      {
	 Site = CreateSU2SpinSite(Spin);
      }
      else if (Model == "spin-u1")
      {
	 Site = CreateU1SpinSite(Spin);
      }
      else if (Model == "spin-z2")
      {
	 Site = CreateZ2SpinSite(Spin);
      }
      else if (Model == "tj-u1")
      {
	 Site = CreateU1tJSite();
      }
      else if (Model == "sf-u1")
      {
	 Site = CreateU1SpinlessFermion();
      }
      else if (Model == "bh-u1")
      {
	 Site = CreateBoseHubbardSpinlessU1Site(NMax);
      }
      else if (Model == "hubbard-u1")
      {
         Site = CreateU1HubbardSite();
      }

      if (!OpL.empty())
	 MyOpL = Site[OpL];
      if (!OpL2.empty())
	 MyOpL2 = Site[OpL2];
      if (!OpR.empty())
	 MyOpR = Site[OpR];

      QuantumNumbers::QuantumNumberList OpLTransformsAs(1, MyOpL.TransformsAs());
      if (OpL2 != "I")
      {
	 OpLTransformsAs = transform_targets(MyOpL.TransformsAs(), MyOpL2.TransformsAs());
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      pvalue_ptr<InfiniteWavefunction> Psi1 
         = pheap::OpenPersistent(PsiStr, mp_pheap::CacheSize(), true);

      SimpleOperator UnitOp;
      if (OpStr == "I")
         UnitOp = SimpleOperator::make_identity(Psi1->Psi.get_front().LocalBasis());
      else
         UnitOp = Site[OpStr];

      // get the set of quantum numbers to show
      typedef std::set<QuantumNumbers::QuantumNumber> QSetType;
      QSetType QL;
      if (vm.count("quantumnumber") != 0)
      {
	 QL.insert(QuantumNumbers::QuantumNumber(Psi1->GetSymmetryList(), Target));
      }
      else
      {
	 QL = QuantumNumbersInBasis(Psi1->C_right.Basis1());
      }
      // show the title
      std::cout << "#sector       #n        #eigenvalue";
      if (!OpL.empty())
      {
	 if (NormOnly)
	    std::cout << "    #overlap";
	 else
	    std::cout << "    #overlap-real      #overlap-imag";
      }
      std::cout << '\n';

      for (QSetType::const_iterator qI = QL.begin(); qI != QL.end(); ++qI)
      {
	 LinearAlgebra::Vector<MatrixOperator> RightEigenvectors;
	 LinearAlgebra::Vector<MatrixOperator> LeftEigenvectors;
         LinearAlgebra::Vector<std::complex<double> > EValues =
            make_spectrum(*Psi1, *qI, MaxEigen, UnitOp, Tol, KrylovLength, Symmetric, &RightEigenvectors, 
                          &LeftEigenvectors, Verbose);

         LinearAlgebra::Vector<std::complex<double> > Overlaps(size(EValues), 0.0);
	 bool HaveOverlaps = false;

	 if (!OpL.empty() && std::find(OpLTransformsAs.begin(), OpLTransformsAs.end(), *qI)
             != OpLTransformsAs.end())
	 {
	    // calculate the MPS representation of the operator
	    
	    InfiniteWavefunction x = *Psi1;
	    MatrixOperator x_unit =  x.C_right * delta_shift(InvertDiagonal(x.C_old, InverseTol), adjoint(x.QShift));
	    LinearWavefunction xPsi = x.Psi;
	    xPsi.set_back(prod(xPsi.get_back(), x_unit));
	    if (Symmetric)
	    {
	       MatrixOperator LambdaSqrt = SqrtDiagonal(x.C_old);
	       MatrixOperator LambdaInvSqrt = InvertDiagonal(LambdaSqrt, InverseTol);
	       xPsi.set_back(prod(xPsi.get_back(), delta_shift(LambdaSqrt, adjoint(x.QShift))));
	       xPsi.set_front(prod(LambdaInvSqrt, xPsi.get_front()));
	    }
	    MatrixOperator Op = MatrixOperator::make_identity(xPsi.Basis1());
	    if (Symmetric)
	       Op = x.C_old;
	    LinearWavefunction::const_iterator I = xPsi.begin();
	    Op = adjoint(operator_prod(herm(MyOpL), herm(*I), Op, *I));
	    ++I;
	    if (I == xPsi.end() && OpL2 != "I")
	    {
	       // for the case of two-site operator with a one-site unit cell
	       I = xPsi.begin();
	    }
	    if (I != xPsi.end())
	    {
	       Op = adjoint(operator_prod(herm(MyOpL2), herm(*I), adjoint(Op), *I, *qI));
	       ++I;
	    }
	    while (I != xPsi.end())
	    {
	       Op = operator_prod(herm(UnitOp), herm(*I), Op, *I);
	       ++I;
	    }
	    //	    TRACE(Op); 
            //	    TRACE(norm_frob_sq(Op));
	    double OpNorm = 1;//norm_frob(Op);
	    //	    std::cout << "#overlaps\n";
	    for (unsigned i = 0; i < size(RightEigenvectors); ++i)
	    {
	       Overlaps[i] = inner_prod(RightEigenvectors[i], Op) / OpNorm;
	    }
	    HaveOverlaps = true;
	 }

         for (int i = 0; i < int(size(EValues)); ++i)
         {
            std::cout << std::setw(7) << boost::lexical_cast<std::string>(*qI) << "  "
                      << std::setw(6) << i << "  ";
	    if (ShowMagnitude)
               std::cout << std::setw(17) << norm_frob(EValues[i]) << "  ";
            else
               std::cout << std::setw(17) << EValues[i].real() << "  "
                         << std::setw(17) << EValues[i].imag() << "  ";

	    if (HaveOverlaps)
	    {
	       if (NormOnly)
		  std::cout << std::setw(17) << norm_frob(Overlaps[i]) << '\n';
	       else
		  std::cout << std::setw(17) << Overlaps[i].real() << "  "
			    << std::setw(17) << Overlaps[i].imag() << '\n';
	    }
	    else
	       std::cout << '\n';
         }

#if 0
	 LinearAlgebra::Matrix<std::complex<double> > Overlaps(size(Eigenvectors), size(Eigenvectors));
	 for (unsigned i = 0; i < size(Eigenvectors); ++i)
	 {
	    for (unsigned j = i; j < size(Eigenvectors); ++j)
	    {
	       Overlaps(i,j) = inner_prod(Eigenvectors[i], Eigenvectors[j]);
	       Overlaps(j,i) = conj(Overlaps(i,j));
	    }
	 }

	 TRACE(Overlaps);
	 TRACE(transform(Overlaps, LinearAlgebra::NormFrob<std::complex<double> >()));
#endif
      }

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
