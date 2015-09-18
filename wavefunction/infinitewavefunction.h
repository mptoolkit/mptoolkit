// -*- C++ -*-
//
// InfiniteWavefunction: class to represent a linear `infinite' matrix product wavefunction
// in canonical form.
//

#if !defined(MPTOOLKIT_WAVEFUNCTION_INFINITEWAVEFUNCTION_H)
#define MPTOOLKIT_WAVEFUNCTION_INFINITEWAVEFUNCTION_H

#include "wavefunction/canonicalwavefunction.h"
#include "wavefunction/linearwavefunction.h"
#include "mpo/finite_mpo.h"
#include "mpo/product_mpo.h"
#include <boost/tuple/tuple.hpp>

// class to represent an infinite wavefunction in the left canonical basis
class InfiniteWavefunctionLeft : public CanonicalWavefunctionBase
{
   public:
      InfiniteWavefunctionLeft() {}

      // construction from a LinearWavefunction (in left-canonical form with lambda matrix on the right)
      InfiniteWavefunctionLeft(LinearWavefunction const& Psi, MatrixOperator const& Lambda, 
			       QuantumNumbers::QuantumNumber const& QShift_);

      // constructs and canonicalizes the wavefunction
      InfiniteWavefunctionLeft(LinearWavefunction const& Psi, QuantumNumbers::QuantumNumber const& QShift_);

      InfiniteWavefunctionLeft(InfiniteWavefunctionLeft const& Psi) : CanonicalWavefunctionBase(Psi), QShift(Psi.QShift) {}

      InfiniteWavefunctionLeft& operator=(InfiniteWavefunctionLeft const& Psi)
      { CanonicalWavefunctionBase::operator=(Psi); QShift = Psi.QShift; return *this; }

      RealDiagonalOperator lambda_r() const { return this->lambda(this->size()-1); }

      QuantumNumber qshift() const { return QShift; }

      // Rotates the wavefunction to the left, by taking the left-most site and moving it to the right
      void rotate_left(int Count);

      // Rotates the wavefunction to the left, by taking the right-most site and moving it to the left
      void rotate_right(int Count);

      // returns the orthogonality fidelity.  Normally this should be epsilon
      double orthogonality_fidelity() const;

      void check_structure() const;
      void debug_check_structure() const;

      static PStream::VersionTag VersionT;

      friend PStream::ipstream& operator>>(PStream::ipstream& in, InfiniteWavefunctionLeft& Psi);
      friend PStream::opstream& operator<<(PStream::opstream& out, InfiniteWavefunctionLeft const& Psi);
      friend void read_version(PStream::ipstream& in, InfiniteWavefunctionLeft& Psi, int Version);

   private:
      void Initialize(LinearWavefunction const& Psi, MatrixOperator const& Lambda);

      QuantumNumber QShift;

      friend void inplace_reflect(InfiniteWavefunctionLeft& Psi);
      friend void inplace_conj(InfiniteWavefunctionLeft& Psi);
      friend InfiniteWavefunctionLeft repeat(InfiniteWavefunctionLeft const& Psi, int Count);
};

// Convert a, infinite wavefunction to left-canonical form, 
// and returns the Lambda matrix on the right-hand-side.
RealDiagonalOperator
left_canonicalize(LinearWavefunction& Psi, QuantumNumbers::QuantumNumber const& QShift);

// Extend the unit cell of the wavefunction by repeating it Count number of times
InfiniteWavefunctionLeft repeat(InfiniteWavefunctionLeft const& Psi, int Count);

// returns a linear wavefunction that is in pure left-orthogonal form.
// This is a very fast operation that only manipulates pvalue_handle objects.
std::pair<LinearWavefunction, RealDiagonalOperator>
get_left_canonical(InfiniteWavefunctionLeft const& Psi);

// returns a linear wavefunction that is in pure right-orthogonal form.
// The asymmetry in the return type between get_left_canonical and get_right_canonical
// is because for a left-canonical wavefunction, the Lambda matrix is already in diagonal form,
// whereas converting a left-canonical wavefunction into a right-canonical wavefunction give
// an additional unitary matrix, that we cannot wrap around to the other end of the wavefunction
// as it would change the basis there.
// This function does an SVD on each MPS
boost::tuple<MatrixOperator, RealDiagonalOperator, LinearWavefunction>
get_right_canonical(InfiniteWavefunctionLeft const& Psi);

// function to extract the local basis (as a vector of BasisList) from a wavefunction
// TODO: this implementation isn't terrible, but isn't so elegant either!
std::vector<BasisList>
inline
ExtractLocalBasis(InfiniteWavefunctionLeft const& Psi)
{
   return ExtractLocalBasis(get_left_canonical(Psi).first);
}

#if 0
class InfiniteRightWavefunctionRight : public CanonicalWavefunctionBase
{
   public:
      // construction and orthogonalization from a LinearWavefunction
      explicit InfiniteRightCanonicalWavefunction(LinearWavefunction const& Psi);

      InfiniteRightCanonicalWavefunction(InfiniteRightCanonicalWavefunction const& Psi);

      InfiniteRightCanonicalWavefunction& operator=(InfiniteRightCanonicalWavefunction const& Psi);

      void check_structure() const;
      void debug_check_structure() const;

      friend PStream::ipstream& operator>>(PStream::ipstream& in, RightCanonicalWavefunction& Psi);
      friend PStream::opstream& operator<<(PStream::opstream& out, RightCanonicalWavefunction const& Psi);
};
#endif


// When taking the inverse of the singular values, we need a cutoff value beyond which we ignore the
// small elements.  This is set to the environment variable MP_INVERSE_TOL, or if that variable
// is not defined, it is set to InverseTolDefault.
double const InverseTolDefault = 1E-7;
extern double const InverseTol;

// calculates the overlap of two iMPS, per unit cell.
// The eigenvector can be in any allowable symmetry sector.
std::complex<double> overlap(InfiniteWavefunctionLeft const& x,  InfiniteWavefunctionLeft const& y,
                             QuantumNumbers::QuantumNumber const& Sector, 
			     int Iter = 20, double Tol = 1E-12, int Verbose = 0);

// This version allows a string operator also
std::complex<double> overlap(InfiniteWavefunctionLeft const& x, FiniteMPO const& StringOp,
                             InfiniteWavefunctionLeft const& y,
                             QuantumNumbers::QuantumNumber const& Sector, 
			     int Iter = 20, double Tol = 1E-12, int Verbose = 0);

std::complex<double> overlap(InfiniteWavefunctionLeft const& x, ProductMPO const& StringOp,
                             InfiniteWavefunctionLeft const& y,
                             QuantumNumbers::QuantumNumber const& Sector, 
			     int Iter = 20, double Tol = 1E-12, int Verbose = 0);

// Reflect a wavefunction in place
void inplace_reflect(InfiniteWavefunctionLeft& Psi);

// Conjugate a wavefunction in place
void inplace_conj(InfiniteWavefunctionLeft& Psi);

// Spatial reflection of a wavefunction
//InfiniteWavefunctionRight reflect(InfiniteWavefunction const& Psi);

// version of reflect where we apply a local operator also
//InfiniteWavefunctionRight reflect(InfiniteWavefunction const& Psi, std::vector<SimpleOperator> const& Op);


#endif
