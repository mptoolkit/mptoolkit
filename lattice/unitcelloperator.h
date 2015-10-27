// -*- C++ -*- $Id$
//
// A proxy class to represent a UnitCellMPO

#if !defined(MPTOOLKIT_LATTICE_UNITCELLOPERATOR_H)
#define MPTOOLKIT_LATTICE_UNITCELLOPERATOR_H

#include "lattice/unitcell.h"

// Local operators, using the notation Op(cell)(index),
// are always r-values, since we cannot redefine them.  We
// must return a local operator as a UnitCellMPO by value, since it is constructed
// as needed.
//
// Finite operators, using the notation Op(cell), are also r-values.
//
// We can use a UnitCellOperator as an L-value only in the form of Op = UnitCellMPO

class UnitCellOperator;

class UnitCellOperatorAtCell
{
   public:
      UnitCellOperatorAtCell() = delete;
      UnitCellOperatorAtCell& operator=(UnitCellOperatorAtCell const&) = delete;

      UnitCellOperatorAtCell(UnitCell const& Cell_, std::string const& Name_, int n_);

      // returns the local operator at the given site
      UnitCellMPO operator[](int i) const;

      // conversion to UnitCellMPO, only possible if the unit cell is 1 site
      operator UnitCellMPO() const;

   private:
      UnitCellOperatorAtCell(UnitCellOperatorAtCell const&) = default;
      UnitCellOperatorAtCell(UnitCellOperatorAtCell&&) = default;

      UnitCell const* Cell;
      std::string Name;
      int n;

      friend class UnitCellOperator;  // because UnitCellOperator needs the copy constructor
};

class UnitCellOperator
{
   public:
      UnitCellOperator() = delete;
      UnitCellOperator(UnitCellOperator const&) = delete;
      UnitCellOperator& operator=(UnitCellOperator const&) = delete;

      UnitCellOperator(UnitCell& Cell_, std::string const& Name_);

      UnitCellOperatorAtCell operator()(int n) const;

      // Returns a local operator at the given site number in the 0'th cell
      UnitCellMPO operator[](int i) const;

      // conversion to a UnitCellMPO, as an operator spanning the 0'th cell
      // We can't return a const& in the const version here, because if
      // the unit cell is 1-site we do an implicit conversion to a local operator,
      // and therefore the UnitCellMPO is a temporary.
      operator UnitCellMPO() const;
      operator UnitCellMPO&();

      UnitCellOperator& operator=(UnitCellMPO const& Op);

      // we could define this - conversion from complex to an operator
      //      UnitCellOperator& operator=(std::complex<double> c);

   private:
      UnitCell* Cell;
      std::string Name;
};

#include "unitcelloperator.cc"

#endif

