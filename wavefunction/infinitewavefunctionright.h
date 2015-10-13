// -*- C++ -*-
//
// InfiniteWavefunctionRight: class to represent a linear `infinite' matrix product wavefunction
// in right canonical form.
//

#if !defined(MPTOOLKIT_WAVEFUNCTION_INFINITEWAVEFUNCTIONRIGHT_H)
#define MPTOOLKIT_WAVEFUNCTION_INFINITEWAVEFUNCTIONRIGHT_H

#include "wavefunction/canonicalwavefunction.h"
#include "wavefunction/linearwavefunction.h"
#include "mpo/finite_mpo.h"
#include "mpo/product_mpo.h"
#include <boost/tuple/tuple.hpp>

// class to represent an infinite wavefunction in the left canonical basis
class InfiniteWavefunctionRight : public CanonicalWavefunctionBase
{
   public:
      InfiniteWavefunctionRight() {}

      // construction from a LinearWavefunction (in right-canonical form with 
      // lambda matrix on the left)
      InfiniteWavefunctionRight(MatrixOperator const& Lambda, LinearWavefunction const& Psi, 
				QuantumNumbers::QuantumNumber const& QShift_);

      // constructs and canonicalizes the wavefunction
      InfiniteWavefunctionRight(LinearWavefunction const& Psi, 
				QuantumNumbers::QuantumNumber const& QShift_);

      InfiniteWavefunctionRight(InfiniteWavefunctionRight const& Psi) 
	 : CanonicalWavefunctionBase(Psi), QShift(Psi.QShift) {}

      InfiniteWavefunctionRight& operator=(InfiniteWavefunctionRight const& Psi)
      { CanonicalWavefunctionBase::operator=(Psi); QShift = Psi.QShift; return *this; }

      QuantumNumber qshift() const { return QShift; }

      // Rotates the wavefunction to the left, by taking the left-most site and moving 
      // it to the right
      void rotate_left(int Count);

      // Rotates the wavefunction to the left, by taking the right-most site and moving 
      // it to the left
      void rotate_right(int Count);

      // returns the orthogonality fidelity.  Normally this should be epsilon
      double orthogonality_fidelity() const;

      static PStream::VersionTag VersionT;

      friend PStream::ipstream& operator>>(PStream::ipstream& in, 
					   InfiniteWavefunctionRight& Psi);
      friend PStream::opstream& operator<<(PStream::opstream& out, 
					   InfiniteWavefunctionRight const& Psi);
      friend void read_version(PStream::ipstream& in, InfiniteWavefunctionRight& Psi, int Version);

      void check_structure();
      void debug_check_structure();

   private:
      void Initialize(MatrixOperator const& Lambda, LinearWavefunction const& Psi);

      QuantumNumber QShift;

      // All functions that can modify the internal representation but preserve the canonical form
      // are friend functions.  This is so that we have a central list of such functions,
      // so can update them if the class changes.
      friend void inplace_reflect(InfiniteWavefunctionRight& Psi);
      friend void inplace_conj(InfiniteWavefunctionRight& Psi);
      friend InfiniteWavefunctionRight repeat(InfiniteWavefunctionRight const& Psi, int Count);
};

class InfiniteWavefunctionLeft;

// Convert a, infinite wavefunction to right-canonical form, 
// and returns the Lambda matrix on the right-hand-side.
RealDiagonalOperator
right_canonicalize(LinearWavefunction& Psi, QuantumNumbers::QuantumNumber const& QShift);

// Extend the unit cell of the wavefunction by repeating it Count number of times
InfiniteWavefunctionRight repeat(InfiniteWavefunctionRight const& Psi, int Count);

// returns a linear wavefunction that is in pure right-orthogonal form.
// This is a very fast operation that only manipulates pvalue_handle objects.
std::pair<LinearWavefunction, RealDiagonalOperator>
get_right_canonical(InfiniteWavefunctionRight const& Psi);

// returns a linear wavefunction that is in pure left-orthogonal form.
// The asymmetry in the return type between get_right_canonical and get_left_canonical
// is because for a right-canonical wavefunction, the Lambda matrix is already in diagonal form,
// whereas converting a right-canonical wavefunction into a left-canonical wavefunction gives
// an additional unitary matrix, that we cannot wrap around to the other end of the wavefunction
// as it would change the basis there.
// This function does an SVD on each MPS.
// Often the caller may want to construct
// U*D*herm(U), and U*Psi, as being the right canonical wavefunction in the same basis as Psi.
boost::tuple<MatrixOperator, RealDiagonalOperator, LinearWavefunction>
get_left_canonical(InfiniteWavefunctionRight const& Psi);

// function to extract the local basis (as a vector of BasisList) from a wavefunction
// TODO: this implementation isn't terrible, but isn't so elegant either!
std::vector<BasisList>
inline
ExtractLocalBasis(InfiniteWavefunctionRight const& Psi)
{
   return ExtractLocalBasis(get_right_canonical(Psi).first);
}


// Reflect a wavefunction in place
void inplace_reflect(InfiniteWavefunctionRight& Psi);

// Conjugate a wavefunction in place
void inplace_conj(InfiniteWavefunctionRight& Psi);

// Spatial reflection of a wavefunction
InfiniteWavefunctionLeft reflect(InfiniteWavefunctionRight const& Psi);

// version of reflect where we apply a local operator also
//InfiniteWavefunctionRight reflect(InfiniteWavefunctionLeft const& Psi, 
// std::vector<SimpleOperator> const& Op);

inline
void
InfiniteWavefunctionRight::debug_check_structure()
{
#if !defined(NDEBUG)
   this->check_structure();
#endif
}

#endif
