// -*- C++ -*-
//
// InfiniteWavefunction: class to represent a linear `infinite' matrix product wavefunction
// in canonical form.
//

#if !defined(MPTOOLKIT_WAVEFUNCTION_INFINITEWAVEFUNCTIONLEFT_H)
#define MPTOOLKIT_WAVEFUNCTION_INFINITEWAVEFUNCTIONLEFT_H

#include "wavefunction/canonicalwavefunction.h"
#include "wavefunction/linearwavefunction.h"
#include "mpo/finite_mpo.h"
#include "mpo/product_mpo.h"
#include <boost/tuple/tuple.hpp>

class AttributeList;

// class to represent an infinite wavefunction in the left canonical basis
class InfiniteWavefunctionLeft : public CanonicalWavefunctionBase
{
   public:
      InfiniteWavefunctionLeft() {}

      // named constructors

      // construction from a LinearWavefunction (in left-canonical form with lambda 
      // matrix on the right)
      static
      InfiniteWavefunctionLeft ConstructFromOrthogonal(LinearWavefunction const& Psi, 
						       MatrixOperator const& Lambda, 
						       QuantumNumbers::QuantumNumber const& QShift_,
						       int Verbose = 0);

      // construct and orthogonalize from a LinearWavefunction
      static
      InfiniteWavefunctionLeft Construct(LinearWavefunction const& Psi, 
					 QuantumNumbers::QuantumNumber const& QShift, 
					 int Verbose = 0);

      // construct and orthogonalize from a LinearWavefunction, with an approximation
      // for the right-most density matrix
      static
      InfiniteWavefunctionLeft Construct(LinearWavefunction const& Psi, 
					 MatrixOperator const& GuessRho, 
					 QuantumNumbers::QuantumNumber const& QShift, 
					 int Verbose = 0);

      InfiniteWavefunctionLeft(InfiniteWavefunctionLeft const& Psi) 
	 : CanonicalWavefunctionBase(Psi), QShift(Psi.QShift) {}

      InfiniteWavefunctionLeft& operator=(InfiniteWavefunctionLeft const& Psi)
      { CanonicalWavefunctionBase::operator=(Psi); QShift = Psi.QShift; return *this; }

      QuantumNumber qshift() const { return QShift; }

      // Rotates the wavefunction to the left, by taking the left-most site 
      // and moving it to the right
      void rotate_left(int Count);

      // Rotates the wavefunction to the left, by taking the right-most site 
      // and moving it to the left
      void rotate_right(int Count);

      // returns the orthogonality fidelity.  Normally this should be epsilon
      double orthogonality_fidelity() const;

      void SetDefaultAttributes(AttributeList& A) const;

      static PStream::VersionTag VersionT;

      friend PStream::ipstream& operator>>(PStream::ipstream& in, InfiniteWavefunctionLeft& Psi);
      friend PStream::opstream& operator<<(PStream::opstream& out, 
					   InfiniteWavefunctionLeft const& Psi);
      friend void read_version(PStream::ipstream& in, InfiniteWavefunctionLeft& Psi, int Version);

      void check_structure() const;
      void debug_check_structure() const;

   private:
      explicit InfiniteWavefunctionLeft(QuantumNumber const& QShift_);

      void Initialize(LinearWavefunction const& Psi, MatrixOperator const& Lambda, int Verbose);

      QuantumNumber QShift;

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

      friend InfiniteWavefunctionLeft repeat(InfiniteWavefunctionLeft const& Psi, int Count);

      friend InfiniteWavefunctionLeft wigner_project(InfiniteWavefunctionLeft const& Psi,
						    SymmetryList const& FinalSL);
      friend InfiniteWavefunctionLeft ReorderSymmetry(InfiniteWavefunctionLeft const& Psi, 
						      SymmetryList const& NewSL);
};

class InfiniteWavefunctionRight;

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
// This function does an SVD on each MPS.
// Often the caller may want to construct
// U*D*herm(U), and U*Psi, as being the right canonical wavefunction in the same basis as Psi.
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

// When taking the inverse of the singular values, we need a cutoff value beyond which we ignore the
// small elements.  This is set to the environment variable MP_INVERSE_TOL, or if that variable
// is not defined, it is set to InverseTolDefault.
double const InverseTolDefault = 1E-7;
extern double const InverseTol;

// calculates the overlap of two iMPS, per unit cell.
// The eigenvector can be in any allowable symmetry sector.
// x and y must have the same size
std::pair<std::complex<double>, StateComponent>
overlap(InfiniteWavefunctionLeft const& x,  InfiniteWavefunctionLeft const& y,
	QuantumNumbers::QuantumNumber const& Sector, 
	int Iter = 20, double Tol = 1E-12, int Verbose = 0);

// This version allows a string operator also.
// This version is deprecated.
std::complex<double> overlap(InfiniteWavefunctionLeft const& x, FiniteMPO const& StringOp,
                             InfiniteWavefunctionLeft const& y,
                             QuantumNumbers::QuantumNumber const& Sector, 
			     int Iter = 20, double Tol = 1E-12, int Verbose = 0);

// This version allows the wavefunctions and operator to have different sizes.
// The overlap is returned as a quantity per length, which is the lowest
// common multiple of x.size(), y.size(), StringOp.size()
std::tuple<std::complex<double>, int, StateComponent>
overlap(InfiniteWavefunctionLeft const& x, ProductMPO const& StringOp,
	InfiniteWavefunctionLeft const& y,
	QuantumNumbers::QuantumNumber const& Sector, 
	int Iter = 20, double Tol = 1E-12, int Verbose = 0);

// Reflect a wavefunction in place
void inplace_reflect(InfiniteWavefunctionLeft& Psi);

// Conjugate a wavefunction in place
void inplace_conj(InfiniteWavefunctionLeft& Psi);

// Spatial reflection of a wavefunction
InfiniteWavefunctionRight reflect(InfiniteWavefunctionLeft const& Psi);

// version of reflect where we apply a local operator also
//InfiniteWavefunctionRight reflect(InfiniteWavefunction const& Psi, 
// std::vector<SimpleOperator> const& Op);

inline
void
InfiniteWavefunctionLeft::debug_check_structure() const
{
#if !defined(NDEBUG)
   this->check_structure();
#endif
}

#endif
