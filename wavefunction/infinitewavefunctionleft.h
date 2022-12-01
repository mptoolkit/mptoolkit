// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// wavefunction/infinitewavefunctionleft.h
//
// Copyright (C) 2015-2020 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
// InfiniteWavefunction: class to represent a linear `infinite' matrix product wavefunction
// in canonical form.
//

#if !defined(MPTOOLKIT_WAVEFUNCTION_INFINITEWAVEFUNCTIONLEFT_H)
#define MPTOOLKIT_WAVEFUNCTION_INFINITEWAVEFUNCTIONLEFT_H

#include "wavefunction/canonicalwavefunction.h"
#include "wavefunction/linearwavefunction.h"
#include "mpo/basic_finite_mpo.h"
#include "mpo/product_mpo.h"
#include <tuple>

class AttributeList;
class InfiniteWavefunctionRight;

// class to represent an infinite wavefunction in the left canonical basis
class InfiniteWavefunctionLeft : public CanonicalWavefunctionBase
{
   public:
      InfiniteWavefunctionLeft() = default;
      InfiniteWavefunctionLeft(InfiniteWavefunctionLeft const& x) = default;
      InfiniteWavefunctionLeft(InfiniteWavefunctionLeft&& x) = default;
      InfiniteWavefunctionLeft& operator=(InfiniteWavefunctionLeft const& Psi) = default;
      InfiniteWavefunctionLeft& operator=(InfiniteWavefunctionLeft&& Psi) = default;

      // named constructors

      // construction from a LinearWavefunction that is in left orthogonal form with a diagonal right Lambda matrix
      static
      InfiniteWavefunctionLeft
      ConstructFromOrthogonal(LinearWavefunction Psi,
                              QuantumNumbers::QuantumNumber const& QShift_,
                              RealDiagonalOperator const& Lambda,
                              double LogAmplitude = 0.0,
                              int Verbose = 0);

      // construct and orthogonalize from a LinearWavefunction.  The amplitude of Psi isn't used,
      // the final amplitude is set by the LogAmplitude.
      static
      InfiniteWavefunctionLeft
      Construct(LinearWavefunction Psi,
                QuantumNumbers::QuantumNumber const& QShift,
                double LogAmplitude = 0.0,
                int Verbose = 0);

      // construct and orthogonalize from a LinearWavefunction, with an approximation
      // for the right-most density matrix
      static
      InfiniteWavefunctionLeft
      Construct(LinearWavefunction Psi,
               QuantumNumbers::QuantumNumber const& QShift,
               MatrixOperator GuessRho,
               double LogAmplitude = 0.0,
               int Verbose = 0);

      // construct and orthogonalize from a LinearWavefunction, and incorporate the transfer matrix eigenvalue
      // into the LogAmplitude
      static
      InfiniteWavefunctionLeft
      ConstructPreserveAmplitude(LinearWavefunction Psi,
                                 QuantumNumbers::QuantumNumber const& QShift,
                                 double LogAmplitude,
                                 int Verbose = 0);

      // construct and orthogonalize from a LinearWavefunction, with an approximation
      // for the right-most density matrix
      static
      InfiniteWavefunctionLeft
      ConstructPreserveAmplitude(LinearWavefunction Psi,
                                 QuantumNumbers::QuantumNumber const& QShift,
                                 MatrixOperator GuessRho,
                                 double LogAmplitude,
                                 int Verbose = 0);

      // Conversion from an InfiniteWavefunctionRight
      explicit InfiniteWavefunctionLeft(InfiniteWavefunctionRight const& Psi);

      QuantumNumber qshift() const { return QShift; }

      double log_amplitude() const { return LogAmplitude; }

      // Scale the wavefunction by the complex number x.  This has the effect of adding
      // log(x).real() to the log_amplitude, and rotating the first matrix of the unit cell by
      // x / |x|
      void scale(std::complex<double> x);

      // Rotates the wavefunction to the left, by taking the left-most site
      // and moving it to the right
      void rotate_left(int Count);

      // Rotates the wavefunction to the left, by taking the right-most site
      // and moving it to the left
      void rotate_right(int Count);

      void SetDefaultAttributes(AttributeList& A) const;

      static std::string Type;

      static PStream::VersionTag VersionT;

      friend PStream::ipstream& operator>>(PStream::ipstream& in, InfiniteWavefunctionLeft& Psi);
      friend PStream::opstream& operator<<(PStream::opstream& out,
                                           InfiniteWavefunctionLeft const& Psi);
      friend void read_version(PStream::ipstream& in, InfiniteWavefunctionLeft& Psi, int Version);

      void check_structure() const;
      void debug_check_structure() const;

   private:
      explicit InfiniteWavefunctionLeft(QuantumNumber const& QShift_, double LogAmplitude_);

      // Complete the initialization, given an input Psi that satisfies the left ortho constraint, and a diagonal Lamba_r matrix
      void InitializeFromLeftOrthogonal(LinearWavefunction Psi, RealDiagonalOperator Lambda, int Verbose);

      // The quantum number shift per unit cell
      QuantumNumber QShift;

      // The wavefunction amplitude (log) per unit cell
      double LogAmplitude;

      // All functions that can modify the internal representation but preserve the canonical form
      // are friend functions.  This is so that we have a central list of such functions,
      // so can update them if the class changes.

      // spatial reflection of the wavefunction in-place.  The Basis1() / Basis2() of the reflected wavefunctions
      // are exactly the flip_conj() of the original Basis1/Basis2, with no change in gauge across the cell boundary
      // (that is, if the wavefunction is written in a reflection-symmetric basis then Psi' = Psi (up to
      // internal gauge).
      friend void inplace_reflect(InfiniteWavefunctionLeft& Psi);

      // Conjugate the wavefunction in-place
      friend void inplace_conj(InfiniteWavefunctionLeft& Psi);

      friend void inplace_qshift(InfiniteWavefunctionLeft& Psi, QuantumNumbers::QuantumNumber const& Shift);

      friend InfiniteWavefunctionLeft repeat(InfiniteWavefunctionLeft const& Psi, int Count);

      friend InfiniteWavefunctionLeft wigner_project(InfiniteWavefunctionLeft const& Psi,
                                                    SymmetryList const& FinalSL);
      friend InfiniteWavefunctionLeft ReorderSymmetry(InfiniteWavefunctionLeft const& Psi,
                                                      SymmetryList const& NewSL);
};

// Multiplication by a scalar does the same as psi.scale(x)
InfiniteWavefunctionLeft& operator*=(InfiniteWavefunctionLeft& psi, std::complex<double> x);

// Convert an infinite wavefunction to left-orthogonal form.
// This function leaves the left and right basis invariant.
// Return value is the log amplitude from the transfer matrix eigenvalue, Result' = log(sqrt(evalue))
double
left_orthogonalize(LinearWavefunction& Psi, QuantumNumbers::QuantumNumber const& QShift, double tol = 1E-14, int Verbose = 0);

// Take a wavefunction that is already in left orthogonal form, and gauge fix it so that the right
// transfer matrix eigenvector is diagonal.  Return value is the Lambda matrix on the right-hand-side.
RealDiagonalOperator
gauge_fix_left_orthogonal(LinearWavefunction& Psi, QuantumNumbers::QuantumNumber const& QShift, double tol = 1E-14, int Verbose = 0);

// This version takes a guess density operator
RealDiagonalOperator
gauge_fix_left_orthogonal(LinearWavefunction& Psi, QuantumNumbers::QuantumNumber const& QShift, MatrixOperator GuessRho, double tol = 1E-14, int Verbose = 0);

// Extend the unit cell of the wavefunction by repeating it Count number of times
InfiniteWavefunctionLeft repeat(InfiniteWavefunctionLeft const& Psi, int Count);

// returns a linear wavefunction that is in pure left-orthogonal form, and the lambda matrix.
// This is a very fast operation that only manipulates pvalue_handle objects.
// The returned Lambda matrix is in the Basis2()
std::pair<LinearWavefunction, RealDiagonalOperator>
get_left_canonical(InfiniteWavefunctionLeft const& Psi);

// N into 1 coarse-graining of a wavefunction.  The wavefunction size must be a multiple of N.
InfiniteWavefunctionLeft
coarse_grain(InfiniteWavefunctionLeft const& Psi, int N, int Verbose = 0);

// returns a linear wavefunction that is in pure right-orthogonal form.
// This preserves the basis at every bond.
// Cost is O(N) SVD's
std::tuple<RealDiagonalOperator, LinearWavefunction>
get_right_canonical(InfiniteWavefunctionLeft const& Psi);

// function to extract the local basis (as a vector of BasisList) from a wavefunction
// TODO: this implementation isn't terrible, but isn't so elegant either!
std::vector<BasisList>
inline
ExtractLocalBasis(InfiniteWavefunctionLeft const& Psi)
{
   return ExtractLocalBasis(get_left_canonical(Psi).first);
}

// When taking the inverse of the singular values, we need a cutoff value beyond which we ignore the
// small elements.  This is set to the environment variable MP_INVERSE_TOL, or if that variable
// is not defined, it is set to InverseTolDefault.
double const InverseTolDefault = 1E-7;
extern double const InverseTol;

// This version allows the wavefunctions and operator to have different sizes.
// The overlap is returned as a quantity per length, which is the lowest
// common multiple of x.size(), y.size(), StringOp.size()
// The length is returned as the second component of the tuple

// This version calculates n eigenvalues
std::tuple<std::vector<std::complex<double>>, int>
overlap(InfiniteWavefunctionLeft const& x, ProductMPO const& StringOp,
	       InfiniteWavefunctionLeft const& y, int n,
	       QuantumNumbers::QuantumNumber const& Sector, bool UseAmplitude = true, double Tol = 1E-12, int Verbose = 0);

inline
std::tuple<std::complex<double>, int>
overlap(InfiniteWavefunctionLeft const& x, ProductMPO const& StringOp,
	       InfiniteWavefunctionLeft const& y,
	       QuantumNumbers::QuantumNumber const& Sector, bool UseAmplitude, double Tol, int Verbose)
{
	auto r = overlap(x, StringOp, y, 1, Sector, UseAmplitude, Tol, Verbose);
   return std::make_tuple(std::get<0>(r)[0], std::get<1>(r));  // Could be improved with C++17
}

// inject_left for a BasicFiniteMPO.  This can have support on multiple wavefunction unit cells
MatrixOperator
inject_left(MatrixOperator const& m,
            InfiniteWavefunctionLeft const& Psi1,
            BasicFiniteMPO const& Op,
            InfiniteWavefunctionLeft const& Psi2);

// Reflect a wavefunction in place
void inplace_reflect(InfiniteWavefunctionLeft& Psi);

// Conjugate a wavefunction in place
void inplace_conj(InfiniteWavefunctionLeft& Psi);

// Spatial reflection of a wavefunction
InfiniteWavefunctionRight reflect(InfiniteWavefunctionLeft const& Psi);

// version of reflect where we apply a local operator also
//InfiniteWavefunctionRight reflect(InfiniteWavefunction const& Psi,
// std::vector<SimpleOperator> const& Op);

// Calculates a normalized expectation value over the wavefunction,
// i.e. <Psi|Op|Psi> / <Psi|Psi>
// Op.size() must be a multiple of Psi.size()
std::complex<double>
expectation(InfiniteWavefunctionLeft const& Psi, BasicFiniteMPO const& Op);

inline
void
InfiniteWavefunctionLeft::debug_check_structure() const
{
#if !defined(NDEBUG)
   this->check_structure();
#endif
}

#endif
