// -*- C++ -*- $Id$

template <typename Visitor>
typename Visitor::result_type
Lattice::apply_visitor(Visitor const& v) const
{
   return Data_.apply_visitor(v);
}

template <typename Visitor>
typename Visitor::result_type
Lattice::apply_visitor(Visitor const& v)
{
   return Data_.apply_visitor(v);
}
