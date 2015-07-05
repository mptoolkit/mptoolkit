// -*- C++ -*- $Id$
//
// A proxy class to represent a UnitCellMPO

#if !defined(MPTOOLKIT_LATTICE_UNITCELLOPERATOR_H)
#define MPTOOLKIT_LATTICE_UNITCELLOPERATOR_H

#include "lattice/unitcell.h"

class UnitCellOperatorAtCell
{
   public:
      UnitCellOperatorAtCell();  // not defined

      UnitCellOperatorAtCell(UnitCell const& Cell_, std::string const& Name_, int n_);

      // returns the local operator at the given site
      UnitCellMPO operator[](int i) const;

      // conversion to UnitCellMPO, only possible if the unit cell is 1 site
      operator UnitCellMPO() const;

   private:
      UnitCell const* Cell;
      std::string Name;
      int n;
};

class UnitCellOperator
{
   public:
      UnitCellOperator();  // not defined

      UnitCellOperator(UnitCell& Cell_, std::string const& Name_);

      UnitCellOperatorAtCell operator()(int n);

      // Returns a local operator at the given site number in the 0'th cell
      UnitCellMPO operator[](int i) const;

      // conversion to a UnitCellMPO, as an operator spanning the 0'th cell
      operator UnitCellMPO&();

      operator UnitCellMPO() const;

      UnitCellOperator& operator=(UnitCellMPO const& Op);

      // we could define this - conversion from complex to an operator
      //      UnitCellOperator& operator=(std::complex<double> c);

   private:
      UnitCell* Cell;
      std::string Name;
};

#include "unitcelloperator.cc"

#endif

