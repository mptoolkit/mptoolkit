// -*- C++ -*- $Id$

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
