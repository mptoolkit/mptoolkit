// -*- C++ -*- $Id$

      const_index1_type index1(Index const& i) const;
 { return LinearAlgebra::index1(this->as_derived(), i); }

      const_index2_type index2(Index const& i) const;
{ return LinearAlgebra::index2(this->as_derived(), i); }

