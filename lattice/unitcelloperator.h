// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// lattice/unitcelloperator.h
//
// Copyright (C) 2014-2016 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Reseach publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER
//
// A proxy class to represent a UnitCellMPO

#if !defined(MPTOOLKIT_LATTICE_UNITCELLOPERATOR_H)
#define MPTOOLKIT_LATTICE_UNITCELLOPERATOR_H

#include "lattice/unitcell.h"

// Local operators, using the notation Op(cell)[index],
// are always r-values, since we cannot redefine them.  We
// must return a local operator as a UnitCellMPO by value, since it is constructed
// as needed.
//
// Finite operators, using the notation Op(cell), can be l-values or r-values.
// The l-value form can be used to define the operator.
//
// We can use a UnitCellOperator as an L-value only in the form of Op = UnitCellMPO

class UnitCellOperator;

class UnitCellOperatorAtCell
{
   public:
      UnitCellOperatorAtCell() = delete;
      UnitCellOperatorAtCell& operator=(UnitCellOperatorAtCell const&) = delete;

      UnitCellOperatorAtCell(UnitCell& Cell_, std::string const& Name_, int n_);

      // use to define the MPO
      UnitCellOperatorAtCell& operator=(UnitCellMPO const& Op);

      // returns the local operator at the given site
      UnitCellMPO operator[](int i) const;

      // conversion to UnitCellMPO, only possible if the unit cell is 1 site
      operator UnitCellMPO() const;

   private:
      UnitCellOperatorAtCell(UnitCellOperatorAtCell const&) = default;
      UnitCellOperatorAtCell(UnitCellOperatorAtCell&&) = default;

      UnitCell* Cell;
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

      // Set (or reset) the description of the operator
      void set_description(std::string s) const;

      // we could define this - conversion from complex to an operator
      //      UnitCellOperator& operator=(std::complex<double> c);

   private:
      UnitCell* Cell;
      std::string Name;
};

#include "unitcelloperator.cc"

#endif
