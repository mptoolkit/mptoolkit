// -*- C++ -*- $Id$

// An InfiniteLattice is defined over a UnitCell, and also defines operators
// with infinite support (ie Triangular mpo's).
//
// We really should also define also a ProductMPO.  This would represent the
// infinite product of operators defined on the unit cell, and also have a QShift, with
// the requirement that Basis1() == delta_shift(Basis2(), QShift).
//
// The unit cell of the operators is allowed to a a multiple of the lattice UnitCell.

#if !defined(MPTOOLKIT_LATTICE_INFINITELATTICE_H)
#define MPTOOLKIT_LATTICE_INFINITELATTICE_H

#include "unitcell.h"
#include <vector>
#include "pheap/pvalueptr.h"

// InfiniteMPO is currently a typedef for TriangularMPO,
// but will 
typedef TriangularMPO InfiniteMPO;

class InfiniteLattice
{
   public:
      typedef InfiniteMPO value_type;       // will eventually become some kind of variant
      typedef std::map<std::string, value_type> operator_map_type;
      typedef operator_map_type::const_iterator const_iterator;

      InfiniteLattice(UnitCell const& uc);

      UnitCell& GetUnitCell() const { return UnitCell_; }

      // returns true if the given InfiniteMPO exists
      bool operator_exists(std::string const& s) const;

      // Lookup the specified InfiniteMPO
      InfiniteMPO const& Operator(std::string const& Op) const;

      // iterators over the infinite support operators
      const_iterator begin() const { return OperatorMap_.begin(); }
      const_iterator end() const { return OperatorMap_.end(); }

      // Add the operator to the lattice
      void add(std::string const& Op, InfiniteMPO const& Op);

   private:
      UnitCell UnitCell_;
      opertor_map_type OperatorMap_;

   friend PStream::opstream& operator<<(PStream::opstream& out, InfiniteLattice const& L);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, InfiniteLattice& L);
};

#endif
