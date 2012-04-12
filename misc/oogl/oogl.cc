// -*- C++ -*- $Id$

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
