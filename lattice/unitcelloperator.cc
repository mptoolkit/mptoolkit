// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// lattice/unitcelloperator.cc
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

// UnitCellOperatorAtCell

inline
UnitCellOperatorAtCell::UnitCellOperatorAtCell(UnitCell& Cell_, std::string const& Name_, int n_)
   : Cell(&Cell_), Name(Name_), n(n_)
{
}

inline
UnitCellOperatorAtCell&
UnitCellOperatorAtCell::operator=(UnitCellMPO const& Op)
{
   Cell->assign_operator(Name, Op, n);
   return *this;
}

inline
UnitCellMPO
UnitCellOperatorAtCell::operator[](int i) const
{
   return Cell->local_operator(Name, n, i);
}

inline
UnitCellOperatorAtCell::operator UnitCellMPO() const
{
   return (*Cell)(Name, n);
}

// UnitCellOperator

inline
UnitCellOperator::UnitCellOperator(UnitCell& Cell_, std::string const& Name_)
   : Cell(&Cell_), Name(Name_)
{
}

inline
UnitCellOperatorAtCell
UnitCellOperator::operator()(int n) const
{
   return UnitCellOperatorAtCell(*Cell, Name, n);
}

inline
UnitCellMPO
UnitCellOperator::operator[](int i) const
{
   return Cell->local_operator(Name, i);
}

inline
UnitCellOperator::operator UnitCellMPO&()
{
   return (*Cell)[Name];
}

inline
UnitCellOperator::operator UnitCellMPO() const
{
   return (*Cell)[Name];
}

inline
UnitCellOperator&
UnitCellOperator::operator=(UnitCellMPO const& Op)
{
   Cell->assign_operator(Name, Op);
   return *this;
}

inline
void
UnitCellOperator::set_description(std::string s) const
{
   (*Cell)[Name].set_description(std::move(s));
}
