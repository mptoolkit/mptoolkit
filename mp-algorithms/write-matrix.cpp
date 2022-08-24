// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/write-matrix.cpp
//
// Copyright (C) 2015-2022 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "write-matrix.h"
#include "tensor/basis.h"
#include "tensor/regularize.h"

MatrixOperator
Regularize(MatrixOperator const& M)
{
   MatrixOperator U = Regularize(M.Basis1());
   MatrixOperator V = Regularize(M.Basis2());
   return U * M * herm(V);
}

StateComponent
Regularize(StateComponent const& M)
{
   MatrixOperator U = Regularize(M.Basis1());
   MatrixOperator V = Regularize(M.Basis2());
   return prod(U, prod(M, herm(V)));
}

//
// MATLAB output format
//

void WriteMatrixFormatMATLAB(std::ostream& out, LinearAlgebra::Matrix<std::complex<double>> const& M, bool Quiet)
{
   out << "complex( [\n";
   if (!Quiet)
      "% real part\n ";
   for (int i = 0; i < M.size1(); ++i)
   {
      if (i != 0)
      {
         out << ";\n    ";
      }
      for (int j = 0; j < M.size2(); ++j)
      {
         if (j != 0)
            out << " ";
         out << M(i,j).real();
      }
   }
   out << "  ],";
   out << "[\n";
   if (!Quiet)
      out << "\n%imaginary part\n  ";
   for (int i = 0; i < M.size1(); ++i)
   {
      if (i != 0)
      {
         out << ";\n    ";
      }
      for (int j = 0; j < M.size2(); ++j)
      {
         if (j != 0)
            out << " ";
         out << M(i,j).imag();
      }
   }
   out << "  ] )";
}

void WriteRealMatrixFormatMATLAB(std::ostream& out, LinearAlgebra::Matrix<std::complex<double>> const& M, bool quiet)
{
   out << "[ ";
   for (int i = 0; i < M.size1(); ++i)
   {
      if (i != 0)
      {
         out << ";\n    ";
      }
      for (int j = 0; j < M.size2(); ++j)
      {
         if (j != 0)
            out << " ";
         out << M(i,j).real();
      }
   }
   out << "  ]";
}

void WriteMPS_MATLAB(std::ostream& out, LinearWavefunction const& Psi, MatrixOperator const& Rho, bool Quiet)
{
   LinearWavefunction::const_iterator I = Psi.begin();
   int Site = 0;
   out << "MPS = cell(1," << Psi.size() << ");\n";
   while (I != Psi.end())
   {
      if (!Quiet)
         out << "% site " << Site << "\n";
      //         out << "[\n";
      out << "MPS{1," << (Site+1) << "} = cat(3, ... \n";
      StateComponent A = Regularize(*I);
      for (int i = 0; i < A.size(); ++i)
      {
         if (!Quiet)
            out << "% site " << Site << "  basis state " << i << "\n";
         if (i != 0)
            out << ", ...\n";
         if (A.Basis1().size() != 1 || A.Basis2().size() != 1)
         {
            throw std::runtime_error("mp-matrix: error: mps has non-trivial symmetries");
         }
         WriteMatrixFormatMATLAB(out, A[i](0,0), Quiet);
         //out << "\n";
      }
      if (!Quiet)
         out << "% end of site " << Site << "\n";
      out << ");\n";

      ++I;
      ++Site;
   }
   out << "\n";
   if (!Quiet)
      out << "% end of MPS\n";

   if (!Quiet)
      out << "\n% density matrix\n";
   out << "RHO = ";
   MatrixOperator Rho1 = Regularize(Rho);
   WriteRealMatrixFormatMATLAB(out, Rho1(0,0), Quiet);
   out << ";\n";
   if (!Quiet)
      out << "% end of density matrix\n";
}

//
// Python output format
//

void WriteMatrixFormatPython(std::ostream& out, LinearAlgebra::Matrix<std::complex<double>> const& M,
                             std::string Prefix, bool Quiet)
{
   out << "np.array([[";
   for (int i = 0; i < M.size1(); ++i)
   {
      if (i != 0)
      {
         out << "],\n" << Prefix << " [";
      }
      for (int j = 0; j < M.size2(); ++j)
      {
         if (j != 0)
            out << ", ";
         out << M(i,j).real() << '+' << M(i,j).imag() << 'j';
      }
   }
   out << "]\n" << Prefix << "])";
}

void WriteRealMatrixFormatPython(std::ostream& out, LinearAlgebra::Matrix<std::complex<double>> const& M,
                                 std::string Prefix, bool quiet)
{
   out << "np.array([[";
   for (int i = 0; i < M.size1(); ++i)
   {
      if (i != 0)
      {
         out << "],\n" << Prefix << " [";
      }
      for (int j = 0; j < M.size2(); ++j)
      {
         if (j != 0)
            out << ", ";
         out << M(i,j).real();
      }
   }
   out << "]\n" << Prefix << "])";
}

void WriteMPS_Python(std::ostream& out, LinearWavefunction const& Psi, MatrixOperator const& Rho, bool Quiet)
{
   LinearWavefunction::const_iterator I = Psi.begin();
   int Site = 0;
   out <<"import numpy as np\n";
   out << "MPS = [\n";
   while (I != Psi.end())
   {
      if (!Quiet)
         out << "# site " << Site << "\n";
      out << " np.array([";
      StateComponent A = Regularize(*I);
      for (int i = 0; i < A.size(); ++i)
      {
         if (!Quiet)
            out << "  # site " << Site << "  basis state " << i << "\n";
         if (i != 0)
            out << " ,\n";
         if (A.Basis1().size() != 1 || A.Basis2().size() != 1)
         {
            throw std::runtime_error("mp-matrix: error: mps has non-trivial symmetries");
         }
         out << "  ";
         WriteMatrixFormatPython(out, A[i](0,0), "  ", Quiet);
         out << "\n";
      }
      if (!Quiet)
         out << "# end of site " << Site << "\n";
      out << "])\n";

      ++I;
      ++Site;
      if (I != Psi.end())
      {
         out << ",\n";
      }
   }
   out << "]\n";
   if (!Quiet)
      out << "# end of MPS\n";

   if (!Quiet)
      out << "\n# density matrix\n";
   out << "RHO = ";
   MatrixOperator Rho1 = Regularize(Rho);
   WriteRealMatrixFormatPython(out, Rho1(0,0), "", Quiet);
   out << "\n";
   if (!Quiet)
      out << "# end of density matrix\n";
}
