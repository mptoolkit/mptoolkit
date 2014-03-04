// -*- C++ -*- $Id$
//
// A UnitCellOperator is an operator that exists within a UnitCell, ie it has finite support
// and is represented by a FiniteMPO and a JordanWignerString.  The JordanWignerString is
// stored as a string, so that we can commute properly with arbitrary other unit cells.
// The JordanWignerString should itself be a UnitCellOperator that has a trivial ("I") 
// `nested` J-W string.

// ** This is currently obsolete, the operators are directly part of the UnitCell.
// We probably want to turn this class into some kind of proxy, similar to the old OperatorAtSite.

#if !defined(MPTOOLKIT_LATTICE_UNITCELLOPERATOR_H)
#define MPTOOLKIT_LATTICE_UNITCELLOPERATOR_H

#include "unitcell.h"
#include "mpo/finite_mpo.h"
#include "mpo/triangular_mpo.h"

class UnitCellOperator
{
   public:

      // returns the number of sites that this operator is defined over.
      // This is always a multiple of the unit cell size.
      int size() const;

      // returns the Jordan-Wigner string associated with this operator
      std::string JordanWignerString() const;

      // Returns a finite version of the operator.  If it is
      // a triangular operator, then close off the boundaries
      // to make it finite on the unit cell.
      // PRECONDITION: operator is finite or triangular
      FiniteMPO AsFiniteMPO() const;

      // returns the UnitCell of this operator
      UnitCell const& GetUnitCell() const { return *pUnitCell; }

      // shorthand for GetUnitCell().size()
      int UnitCellSize() const { return pUnitCell->size(); }

   private:
      pvalue_ptr<UnitCell> pUnitCell;
      std::string JWString;
      boost::variant<FiniteMPO, TriangularMPO> Operator;
};

UnitCellOperator operator+(UnitCellOperator const& A, UnitCellOperator const& B);

UnitCellOperator operator*(UnitCellOperator const& A, UnitCellOperator const& B);

UnitCellOperator operator*(double x, UnitCellOperator const& Op);

UnitCellOperator operator*(std::complex<double> x, UnitCellOperator const& Op);


#endif

