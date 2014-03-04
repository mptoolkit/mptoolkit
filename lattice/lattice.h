// -*- C++ -*- $Id$
//
// A lattice is a UnitCell plus a collection
// of operators.  Local operators are obtained
// as-needed from the UnitCell (as a FiniteMPO).
// Other operators can be defined, either a
// FiniteMPO, TriangularMPO or PeriodicMPO.

#if !defined(LATTICE_H_SDJKHFWUY89AS)
#define LATTICE_H_SDJKHFWUY89AS

#include "lattice_operator.h"

class Lattice
{
   public:

   private:
      std::map<std::string, LatticeOperator> DefinedOperators;

};


#endif
