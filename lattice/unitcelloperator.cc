// -*- C++ -*- $Id$

inline
UnitCellOperatorAtCell::UnitCellOperatorAtCell(UnitCell const& Cell_, std::string const& Name_, int n_)
   : Cell(&Cell_), Name(Name_), n(n_)
{
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


inline
UnitCellOperator::UnitCellOperator(UnitCell& Cell_, std::string const& Name_)
   : Cell(&Cell_), Name(Name_)
{
}

inline
UnitCellOperatorAtCell
UnitCellOperator::operator()(int n)
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
   (*Cell)[Name] = Op;
   return *this;
}
