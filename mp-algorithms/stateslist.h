// -*- C++ -*- $Id$
/*
  readstates.h

  utility class for parsing the NumStates line.

  Created 2001-04-26 Ian McCulloch

  2002-10-27: Added the WeightsList for controlling density matrix weights for multiple wavefunctions.

  2002-12-11: Extended the syntax to allow 'xN' and '+N'.  Usage is, eg 10x5 (equivalent to 10 10 10 10 10),
              or 10x5+5 (equivalent to 10 15 20 25 30).  'x' and '+' fields can occur in any order.
              Also added 'w' as an alternative to '*' to wait for convergence.  '*' is a bit confusing
              when combined with '+'.
	      Also added the Append() method to increase the StatesList on the fly.

  2007-04-30: Added TruncationError and Broadening fields to the states list.  This is
              done with, for example, 100r1e-5, meaning keep max 100 states with truncation error
              1e-5.  If the number of states is zero, then it is set to DefaultMaxStates.
              Broadening works a similar way, with a suffix b1e-5.
*/

#if !defined(READSTATES_H_FDSKJFDS89U8UFRHUFRHUFR895T7Y547THFDS8)
#define READSTATES_H_FDSKJFDS89U8UFRHUFRHUFR895T7Y547THFDS8

#include "common/trace.h"
#include "pstream/pstream.h"
#include "matrixproduct/density.h"  // for definition of StatesInfo
#include <string>
#include <vector>
#include <iostream>

class StatesList
{
   public:
      StatesList() {}

      ~StatesList();

      StatesList(char const* str);
      
      StatesList(std::string const& str);

      void Append(char const* str);

      // overrides the current settings and sets CalculatePhysics() to false for all sweeps
      void OverrideNoPhysics();

      // overrides the current settings and sets CalculateCorrelations() to false for all sweeps
      void OverrideNoCorrelations();

      int NumSweeps() const { return Info.size(); }
      int NumStates(int s) const { RANGE_CHECK(s, 0, int(Info.size()-1)); return Info[s].NumStates; }
      bool CalculatePhysics(int s) const { RANGE_CHECK(s, 0, int(Info.size())-1); return Info[s].Physics; }
      bool CalculateCorrelations(int s) const { RANGE_CHECK(s, 0, int(Info.size())-1); return Info[s].Corrs; }

      double TruncationError(int s) const { RANGE_CHECK(s, 0, int(Info.size())-1); return Info[s].Trunc; }
      double Broadening(int s) const { RANGE_CHECK(s, 0, int(Info.size())-1); return Info[s].Broad; }

      double MixFactor(int s) const { RANGE_CHECK(s, 0, int(Info.size())-1); return Info[s].MixFactor; }

      bool WaitConverge(int s) const { RANGE_CHECK(s, 0, int(Info.size())-1); return Info[s].Wait; }
      bool SaveState(int s) const { RANGE_CHECK(s, 0, int(Info.size())-1); return Info[s].Save; }
      bool TestConverge(int s) const { RANGE_CHECK(s, 0, int(Info.size())-1); return Info[s].Test; }

      StatesInfo GetStatesInfo(StatesInfo const& Default, int s) const;

      // returns true if any physical operators are wanted, on any sweep
      bool AnyPhysics() const;

      // returns true if any correlations are wanted, on any sweep
      bool AnyCorrelations() const;

      void ShowOptions(std::ostream& out) const;

   private:
      struct InfoType
      {
	 int NumStates;
	 bool Physics;
	 bool Corrs;
	 bool Wait;
         bool Test;
	 bool Save;
	 double Trunc;
	 double Broad;
         double MixFactor;
      };

      std::vector<InfoType> Info;

   friend PStream::opstream& operator<<(PStream::opstream& out, InfoType const& i);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, InfoType& i);

   friend PStream::opstream& operator<<(PStream::opstream& out, StatesList const& s);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, StatesList& s);
};

class WeightsList
{
   public:
      WeightsList() {}

      WeightsList(char const* str);

      int NumWeights() const { return Weights.size(); }

      double Weight(int n) const { return Weights[n]; }

      void ShowOptions(std::ostream& out) const;

  private:
      std::vector<double> Weights;

   friend PStream::opstream& operator<<(PStream::opstream& out, WeightsList const& s);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, WeightsList& s);
};

#endif
