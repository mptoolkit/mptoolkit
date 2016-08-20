// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// misc/oogl/oogl.cc
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

namespace oogl
{

template <class VertexListType>
BasicPolyhedra<VertexListType>::BasicPolyhedra()
{
}

template <class VertexListType>
BasicPolyhedra<VertexListType>::BasicPolyhedra(VertexListType const& VList_) : VList(VList_)
{
}

template <class VertexListType>
OoglPrimitive* BasicPolyhedra<VertexListType>::clone() const
{
   return new BasicPolyhedra<VertexListType>(*this);
}

template <class VertexListType>
void BasicPolyhedra<VertexListType>::write(std::ostream& out) const
{
   out << VertexListType::prefix() << "OFF\n"
       << num_vertices() << ' '
       << num_faces() << ' '
       << num_edges() << '\n'
       << vertex_list();
   std::copy(begin(), end(), std::ostream_iterator<Face>(out, "\n"));
}

template <class VertexListType>
void
BasicPolyhedra<VertexListType>::set_vertex_list(VertexListType const& VList_)
{
   VList = VList_;
}

template <class VertexListType>
inline
BasicPolyhedra<VertexListType>&
BasicPolyhedra<VertexListType>::operator*=(Transform const& t)
{
   VList *= t;
   return *this;
}

} // namespace oogl
