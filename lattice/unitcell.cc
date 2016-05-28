// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// lattice/unitcell.cc
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#if 0
template <typename Visitor>
typename Visitor::result_type
UnitCell::apply_for_each_site(Visitor const& v) const
{
   return Data_.apply_visitor(v);
}

template <typename Visitor>
typename Visitor::result_type
UnitCell::apply_for_each_site(Visitor const& v)
{
   return Data_.apply_visitor(v);
}
#endif
