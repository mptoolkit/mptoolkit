// -*- C++ -*-

// IBC wavefunction.
// We store the Left semi-infinite strip (in left-orthogonal form)
// the Right semi-infinite strip (in right-orthogonal form)
// and the Window
// Number of sites of the Left unit cell that have been incorporated
// into the Window.  If this is zero, then the Window is made up
// of complete unit cells.  If it becomes equal to Left.size(), then
// we want to shift quantum numbers of Left and reset WindowLeftSites to zero.
//
// Note: we store Left and Right as InfiniteWavefunction's, but the
// Lambda matrices don't have physical signifance (except a long way away from
// the window), since the MPS doesn't have a canonical form.
// We can use the Lambda matrices only to get fixed points in the semi-infinite
// strip.
//
// The window is stored as a WavefunctionSectionLeft.  This is derived from CanonicalWavefunctionBase,
// and therefore is in the form
//         |                     |
// -lambda-A-lambda- ... -lambda-A-lambda-
// with diagonal Lambda matrices on both boundaries.  Beyond these boundaries, there is necessarily
// a unitary matrix to connect this to the basis of the semi-infinite strips.

#if !defined(MPTOOLKIT_WAVFUNCTION_IBC_H)
#define MPTOOLKIT_WAVFUNCTION_IBC_H

#include "infinitewavefunctionleft.h"
#include "infinitewavefunctionright.h"
#include "canonicalwavefunction.h"

// Note that Basis1(), Basis2() refer to the basis of the edge lambda matrices, NOT the LeftU.Basis1()
// and RightU.Basis2() !
class WavefunctionSectionLeft : public CanonicalWavefunctionBase
{
   public:
      WavefunctionSectionLeft();

      WavefunctionSectionLeft(WavefunctionSectionLeft const& Psi) = default;

      explicit WavefunctionSectionLeft(InfiniteWavefunctionLeft const& Psi);

      // Constructs a WavefunctionSectionLeft from a LinearWavefunction that is
      // in left-orthogonal form, with the Lambda matrix at the right-hand edge
      // (this isn't a proper WavefunctionSectionLeft form until the lambda matrix
      // is re-incorporated into the MPS)
      // We overload on rvalues here.
      // We could also use a CentreWavefunction with no loss of efficiency if we wanted
      static
      WavefunctionSectionLeft ConstructFromLeftOrthogonal(LinearWavefunction&& Psi, 
							  MatrixOperator const& Lambda,
							  int Verbose = 0);

      static
      WavefunctionSectionLeft ConstructFromLeftOrthogonal(LinearWavefunction const& Psi, 
							  MatrixOperator const& Lambda,
							  int Verbose = 0);
      void check_structure() const;
      void debug_check_structure() const;

      MatrixOperator const& LeftU() const { return LeftU_; }
      MatrixOperator const& RightU() const { return RightU_; }

      static PStream::VersionTag VersionT;

   private:
      MatrixOperator LeftU_, RightU_;

      friend void inplace_reflect(WavefunctionSectionLeft& Psi);
      friend void inplace_conj(WavefunctionSectionLeft& Psi);

      friend WavefunctionSectionLeft wigner_project(WavefunctionSectionLeft const& Psi,
						    SymmetryList const& FinalSL);
      friend WavefunctionSectionLeft ReorderSymmetry(WavefunctionSectionLeft const& Psi, 
						     SymmetryList const& NewSL);

      friend PStream::opstream& operator<<(PStream::opstream& out, WavefunctionSectionLeft const& Psi);
      friend PStream::ipstream& operator>>(PStream::ipstream& in, WavefunctionSectionLeft& Psi);
};

std::pair<LinearWavefunction, MatrixOperator>
get_left_canonical(WavefunctionSectionLeft const& Psi);

void inplace_reflect(WavefunctionSectionLeft& Psi);

void inplace_conj(WavefunctionSectionLeft& Psi);

class IBCWavefunction
{
   public:
      IBCWavefunction();

      IBCWavefunction(IBCWavefunction const& Psi) = default;

      IBCWavefunction(InfiniteWavefunctionLeft const& Left_,
		      WavefunctionSectionLeft const& Window_,
		      InfiniteWavefunctionRight const& Right_,
		      int Offset = 0);

      IBCWavefunction(InfiniteWavefunctionLeft const& Left_,
		      WavefunctionSectionLeft const& Window_,
		      InfiniteWavefunctionRight const& Right_,
		      int Offset,
		      int WindowLeft,
		      int WindowRight);

      SymmetryList GetSymmetryList() const { return Window.GetSymmetryList(); }

      int window_size() const { return Window.size(); }

      int window_offset() const { return WindowOffset; }

      void SetDefaultAttributes(AttributeList& A) const;

      static PStream::VersionTag VersionT;

      void check_structure() const;
      void debug_check_structure() const;

      // private:

      // Number of sites of the Left unit cell that have been incorporated into
      // the Window (from 0 .. Left.size()-1)
      int WindowLeftSites;

      // Number of sites of the Right unit cell that have been incorporated into
      // the Window (from 0 .. Right.size()-1)
      int WindowRightSites;

      int WindowOffset;  // site index of the first site of the Window

      InfiniteWavefunctionLeft Left;
      WavefunctionSectionLeft Window;
      InfiniteWavefunctionRight Right;

      friend void inplace_reflect(IBCWavefunction& Psi);
      friend void inplace_conj(IBCWavefunction& Psi);

      friend IBCWavefunction wigner_project(IBCWavefunction const& Psi,
					    SymmetryList const& FinalSL);
      friend IBCWavefunction ReorderSymmetry(IBCWavefunction const& Psi, 
					     SymmetryList const& NewSL);

      friend PStream::opstream& operator<<(PStream::opstream& out, IBCWavefunction const& Psi);
      friend PStream::ipstream& operator>>(PStream::ipstream& in, IBCWavefunction& Psi);

};

// Reflect a wavefunction in place
void inplace_reflect(IBCWavefunction& Psi);

// Conjugate a wavefunction in place
void inplace_conj(IBCWavefunction& Psi);

// Spatial reflection of a wavefunction
IBCWavefunction reflect(IBCWavefunction const& Psi);

inline
void
WavefunctionSectionLeft::debug_check_structure() const
{
#if !defined(DNDEBUG)
   this->check_structure();
#endif
}

#endif
