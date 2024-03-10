// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// lattice/finitelattice.cc
//
// Copyright (C) 2013-2016 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

template <typename Visitor>
typename Visitor::result_type
UnitCell::apply_visitor(Visitor const& v) const
{
   return Data_->apply_visitor(v);
}

template <typename Visitor>
typename Visitor::result_type
UnitCell::apply_visitor(Visitor const& v)
{
   return Data_.mutate()->apply_visitor(v);
}
