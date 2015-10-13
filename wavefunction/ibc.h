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

#if !defined(MPTOOLKIT_WAVFUNCTION_IBC_H)
#define MPTOOLKIT_WAVFUNCTION_IBC_H

#include "infinitewavefunctionleft.h"
#include "infinitewavefunctionright.h"
#include "canonicalwavefunction.h"

class WavefunctionSectionLeft : public CanonicalWavefunctionBase
{
};

class IBCWavefunction
{
   public:
      IBCWavefunction(InfiniteWavefunctionLeft const& Left_,
		      WavefunctionSectionLeft const& Window_,
		      InfiniteWavefunctionRight const& Right_,
		      int Offset);

      // Number of sites of the Left unit cell that have been incorporated into
      // the Window (from 0 .. Left.size()-1)
      int WindowLeftSites;

      // Number of sites of the Right unit cell that have been incorporated into
      // the Window (from 0 .. Right.size()-1)
      int WindowSitesRight;

      int WindowOffset;  // site index of the first site of the Window

      InfiniteWavefunctionLeft Left;
      WavefunctionSectionLeft Window;
      InfiniteWavefunctionRight Right;
};

