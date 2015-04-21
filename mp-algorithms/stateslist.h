// -*- C++ -*- $Id$
/*
  stateslist.h

  utility class for parsing the NumStates line.

  Created 2001-04-26 Ian McCulloch
  Heavily modified 2015-04-21

  a StatesList is a sequence of directives, separated by whitespace or commas.
  A single directive is of the following forms:

  N..MxS             meaning: perform S sweeps, with the number of states increasing uniformly from N to M
  N+IxS              meaning: perform S sweeps, starting with N states and increasing by I each sweel
  NxS                meaning: perform S sweeps with N states in each sweep
  N                  meaning: perform a single sweep with N states.

  If just a single sweep is specified, then additionally flags can be specified, in any order,
  w : wait for convergence (keep sweeping until the convergence criteria is satisfied)
  t : test for convergence (do the convergence test and report only, don't wait if the criteria is not satisfied)
  s : save the wavefunction at the end of the sweep (if combined with w, only save once converged)
  v : calculate the variance at the end of this sweep
*/

#if !defined(MPTOOLKIT_MP_ALGORITHMS_TATESLIST_H)
#define MPTOOLKIT_MP_ALGORITHMS_STATESLIST_H

#include "common/trace.h"
#include "pstream/pstream.h"
#include "mps/density.h"  // for definition of StatesInfo
#include <string>
#include <vector>
#include <iostream>

struct StateParams
{
   // The default constructor initializes the members to the default values
   StateParams();

   int NumStates;
   bool Wait;
   bool Test;
   bool Save;
   bool Variance;

   // Named constructor to continue a sequence of StateParams.
   // Currently this just copies the NumStates and sets everything else false
   static StateParams ContinueFrom(StateParams const& Info);
};

PStream::opstream& operator<<(PStream::opstream& out, StateParams const& Info);
PStream::ipstream& operator>>(PStream::ipstream& in, StateParams& Info);

std::ostream& operator<<(std::ostream& out, StateParams const& Info);

class StatesList
{
   public:
      explicit StatesList(char const* str);
      
      explicit StatesList(std::string const& str);

      StatesList(int NumSweeps, int NumStates);

      void Append(char const* str);

      void Append(std::string const& s);

      // repeat the previous sweep n times (including the exising; add n-1 new sweeps)
      void Repeat(int n);

      int size() const { return Info.size(); }

      StateParams const& operator[](int i) const { return Info[i]; }

   private:
      void AppendToken(char const* str);

      std::vector<StateParams> Info;

   friend PStream::opstream& operator<<(PStream::opstream& out, StatesList const& s);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, StatesList& s);
};

std::ostream& operator<<(std::ostream& out, StatesList const& States);

#endif
