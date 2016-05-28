// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mpo/generic_mpo.cc
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

inline
void
GenericMPO::debug_check_structure() const
{
#if !defined(NDEBUG)
   this->check_structure();
#endif
}
