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
#include "mpo/triangular_mpo.h"
#include <vector>
#include "pheap/pvalueptr.h"

// InfiniteMPO is currently a typedef for TriangularMPO,
// but will eventually be some kind of variant for a TriangularMPO and ProductMPO
typedef TriangularMPO InfiniteMPO;

class InfiniteLattice
{
   public:
      typedef std::map<std::string, TriangularMPO> triangular_map_type;
      typedef triangular_map_type::const_iterator  const_triangular_iterator;

      InfiniteLattice();

      explicit InfiniteLattice(UnitCell const& uc);

      UnitCell const& GetUnitCell() const { return UnitCell_; }

      // returns true if the given InfiniteMPO exists
      bool triangular_operator_exists(std::string const& s) const;

      // Lookup the specified InfiniteMPO
      TriangularMPO const& TriangularOperator(std::string const& Op) const;

      TriangularMPO TriangularOperatorFunction(std::string const& Op,
					       std::vector<std::complex<double> > const& Params) const;	 

      TriangularMPO& operator[](std::string const& Op);
      TriangularMPO const& operator[](std::string const& Op) const;

      // iterators over the infinite support operators
      const_triangular_iterator begin_triangular() const { return OperatorMap_.begin(); }
      const_triangular_iterator end_triangular() const { return OperatorMap_.end(); }

      // Add the operator to the lattice
      void add_triangular(std::string const& Name, InfiniteMPO const& Op);

   private:
      UnitCell UnitCell_;
      triangular_map_type OperatorMap_;

   friend PStream::opstream& operator<<(PStream::opstream& out, InfiniteLattice const& L);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, InfiniteLattice& L);
};

// Constucts a TriangularMPO from the summation over unit cell translations of a finite MPO.
// The Op must have a size() that is a multiple of UnitCellSize, and the local basis of the
/// unit cell must be consistent
TriangularMPO sum_unit(UnitCell const& Cell, FiniteMPO const& Op, LatticeCommute Com);

TriangularMPO sum_unit(UnitCellMPO const& Op);

// Constructs a zero triangular MPO
TriangularMPO make_zero(UnitCell const& Cell);

#endif
