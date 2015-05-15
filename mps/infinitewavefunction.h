// -*- C++ -*- $Id$
//
// InfiniteWavefunction: class to represent a linear `infinite' matrix product wavefunction.
//
// here we mean `infinite' in the iterative sense; a sequence of matrices that can be used to
// generate a translationally invariant state (possibly with a non-trivial unit cell)
// but it is not necessarily at the fixed point.
//
// We store:
// The Center matrix for some arbitrary size lattice, C_0
// the wavefunction for the size plus one unit cell, consisting
// of the right Center matrix C_R and the left-orthonormalized A-matrices B_L.
// We also keep the quantum number that determines the density (ie, quantum number per
// unit cell).  This logically fits in at the left hand of the unit cell.
//
// By convention, we take C_old to be always diagonal.  In practice, 
// this means that C_0.Basis1() == C_0.Basis2().
//
// The right basis for C_0 is the same as the right basis of C_R.  The left basis
// of C_0 is the same as the left basis of the left-most A-matrix, but shifted by
// adjoint(Shift).  This means that (C_0)^{-1} can be applied to the left-most A matrix.
//
// If we rotate the wavefunction to the left and obtain right-orthonormalized matrices B_R
// and the Center matrix C_L (which must be shifted in quantum number)
// then the iDMRG wavefunction transformation for
// the 2*N site lattice is C_R C_0^{-1} C_L.
// At the fixed point, C_0 = C_R = C_L.  I'm not sure if care needs to be taken to get
// the right signs of the density matrix eigenvectors to achieve this condition.
//
// The wavefunction part satisfies DeltaShift(Psi.C_right.Basis2(), QShift) == Psi.Basis1()
// C_old is defined on Psi.Basis1()
// C_right is defined on Psi.Basis2()
// C_old must be diagonal and have C_old.Basis1() == C_old.Basis2()
// C_right may have a different Basis1 to Basis2, and need not be diagonal.

#if !defined(INFINITEWAVEFUNCTION_H_SDJKLCHUIORYH8945Y89Y98)
#define INFINITEWAVEFUNCTION_H_SDJKLCHUIORYH8945Y89Y98

#include "linearwavefunction.h"
#include "mpo/finite_mpo.h"
#include "mpo/product_mpo.h"

class InfiniteWavefunction
{
   public:
      typedef LinearWavefunction::iterator       iterator;
      typedef LinearWavefunction::const_iterator const_iterator;

      InfiniteWavefunction() {}

      QuantumNumbers::SymmetryList GetSymmetryList() const { return C_old.GetSymmetryList(); }

      // size of the unit cell
      int size() const { return Psi.size(); }

      // the target quantum number, per unit cell
      QuantumNumbers::QuantumNumber const& shift() const { return QShift; }

      // returns the left-most basis.
      VectorBasis Basis1() const { return Psi.Basis1(); }
      // returns the right-most basis.
      VectorBasis Basis2() const { return Psi.Basis2(); }

      iterator begin() { return Psi.begin(); }
      iterator end() { return Psi.end(); }

      const_iterator begin() const { return Psi.begin(); }
      const_iterator end() const { return Psi.end(); }

      MatrixOperator C_old;  // By convention, this is kept in a diagonal basis

      QuantumNumbers::QuantumNumber QShift;
      LinearWavefunction Psi;  // left-orthogonalized
      MatrixOperator C_right;

      mutable AttributeList Attr;

   friend PStream::opstream& operator<<(PStream::opstream& out, InfiniteWavefunction const& psi);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, InfiniteWavefunction& psi);
};

// When taking the inverse of the singular values, we need a cutoff value beyond which we ignore the
// small elements.  This is set to the environment variable MP_INVERSE_TOL, or if that variable
// is not defined, it is set to InverseTolDefault.
double const InverseTolDefault = 1E-7;
extern double const InverseTol;

// rotate sites in the unit cell, by taking the left-most site and putting 
// it on the right hand side, repeat Count times
InfiniteWavefunction rotate_left(InfiniteWavefunction const& Psi, int Count);

// Join two wavefunctions together.  Precondition: Psi1.Basis2() == Psi2.Basis1()
InfiniteWavefunction join(InfiniteWavefunction const& Psi1, InfiniteWavefunction const& Psi2);

// Join two wavefunctions together, with a delta shift of Psi1.  Equivalent to
// join(delta_shift(Psi1, Psi2.shift()), Psi2)
// Precondition: Psi1.Basis2() == Psi2.Basis2()
InfiniteWavefunction join_shift(InfiniteWavefunction const& Psi1, InfiniteWavefunction const& Psi2);

// Extend the unit cell of the wavefunction by repeating it Count number of times
InfiniteWavefunction repeat(InfiniteWavefunction const& Psi, int Count);

// returns a linear wavefunction that is in pure left-orthogonal form.
// Note that this is the opposite of the 'normal form' used for finite-size MPS.
// If Psi is not orthonormalized, then the result is unspecified.
LinearWavefunction get_orthogonal_wavefunction(InfiniteWavefunction const& Psi);

// returns a linear wavefunction that is right-orthogonalized
LinearWavefunction get_right_orthogonal_wavefunction(InfiniteWavefunction const& Psi);

// returns the fidelity of the density matrices of C_old and C_right.
// If the fidelity is 1.0, the state is orthonormalized.
double orthogonality_fidelity(InfiniteWavefunction const& x);

// explicitly orthogonalizes the iMPS
void orthogonalize(InfiniteWavefunction& x);

// calculates the overlap of two iMPS, per unit cell.
// The eigenvector can be in any allowable symmetry sector.
std::complex<double> overlap(InfiniteWavefunction const& x,  InfiniteWavefunction const& y,
                             QuantumNumbers::QuantumNumber const& Sector, 
			     int Iter = 20, double Tol = 1E-12, int Verbose = 0);

// This version allows a string operator also
std::complex<double> overlap(InfiniteWavefunction const& x, FiniteMPO const& StringOp,
                             InfiniteWavefunction const& y,
                             QuantumNumbers::QuantumNumber const& Sector, 
			     int Iter = 20, double Tol = 1E-12, int Verbose = 0);

std::complex<double> overlap(InfiniteWavefunction const& x, ProductMPO const& StringOp,
                             InfiniteWavefunction const& y,
                             QuantumNumbers::QuantumNumber const& Sector, 
			     int Iter = 20, double Tol = 1E-12, int Verbose = 0);

// Spatial reflection of a wavefunction
InfiniteWavefunction reflect(InfiniteWavefunction const& Psi);

// version of reflect where we apply a local operator also
InfiniteWavefunction reflect(InfiniteWavefunction const& Psi, std::vector<SimpleOperator> const& Op);


#endif
