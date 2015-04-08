// -*- C++ -*- $Id$

inline
void
GenericMPO::debug_check_structure() const
{
#if !defined(NDEBUG)
   this->check_structure();
#endif
}
