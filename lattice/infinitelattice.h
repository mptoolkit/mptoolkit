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
#include "mpo/product_mpo.h"
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

      typedef std::map<std::string, TriangularMPO> product_map_type;
      typedef triangular_map_type::const_iterator  const_product_iterator;

      InfiniteLattice();

      explicit InfiniteLattice(UnitCell const& uc);

      UnitCell const& GetUnitCell() const { return UnitCell_; }

      // TriangularMPO functions

      // returns true if the given InfiniteMPO exists
      bool triangular_operator_exists(std::string const& s) const;

      // Lookup the named opreator
      TriangularMPO& Triangular(std::string const& Op);
      TriangularMPO const& Triangular(std::string const& Op) const;

      TriangularMPO& operator[](std::string const& Op) { return this->Triangular(Op); }
      TriangularMPO const& operator[](std::string const& Op) const { return this->Triangular(Op); }

      // invoke the given function
      TriangularMPO TriangularOperatorFunction(std::string const& Op,
					       std::vector<std::complex<double> > const& Params) const;	 

      // iterators over the infinite support operators
      const_triangular_iterator begin_triangular() const { return OperatorMap_.begin(); }
      const_triangular_iterator end_triangular() const { return OperatorMap_.end(); }

      // Add the operator to the lattice
      void add_triangular(std::string const& Name, InfiniteMPO const& Op);

      // ProductMPO functions

      bool product_operator_exists(std::string const& s) const;

      ProductMPO const& Product(std::string const& Op) const;
      ProductMPO& Product(std::string const& Op);

      ProductMPO ProductOperatorFunction(std::string const& Op,
					 std::vector<std::complex<double> > const& Params) const;	 

   private:
      UnitCell UnitCell_;
      triangular_map_type OperatorMap_;

   friend PStream::opstream& operator<<(PStream::opstream& out, InfiniteLattice const& L);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, InfiniteLattice& L);
};

// Constucts a TriangularMPO from the summation over unit cell translations of a finite MPO.
// The Op must have a size() that is a multiple of SiteListTypeSize, which must itself be an
// integer multiple of SiteListType.size().
TriangularMPO sum_unit(SiteListType const& SiteList, FiniteMPO const& JW, FiniteMPO const& Op, int UnitCellSize);

TriangularMPO sum_unit(SiteListType const& SiteList, FiniteMPO const& Op, LatticeCommute Com, int UnitCellSize);

TriangularMPO sum_unit(SiteListType const& SiteList, FiniteMPO const& Op, LatticeCommute Com);

TriangularMPO sum_unit(UnitCellMPO const& Op, int UnitCellSize);

TriangularMPO sum_unit(UnitCellMPO const& Op);

// Variant of sum_unit where we add the kink operator (generally will be unitary) to the left hand side
TriangularMPO sum_kink(SiteListType const& SiteList, FiniteMPO const& Kink,
		       FiniteMPO const& Op, LatticeCommute Com, int UnitCellSize);

TriangularMPO sum_kink(UnitCellMPO const& Kink, UnitCellMPO const& Op, int UnitCellSize);

TriangularMPO sum_kink(UnitCellMPO const& Kink, UnitCellMPO const& Op);

TriangularMPO sum_k(SiteListType const& SiteList, std::complex<double> const& k,
		       FiniteMPO const& Op, LatticeCommute Com, int UnitCellSize);

TriangularMPO sum_k(std::complex<double> const& k, UnitCellMPO const& Op, int UnitCellSize);

TriangularMPO sum_k(std::complex<double> const& k, UnitCellMPO const& Op);


// Constructs a zero triangular MPO
TriangularMPO make_zero(SiteListType const& SiteList);

#endif
