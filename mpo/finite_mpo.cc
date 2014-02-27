// -*- C++ -*- $Id$

inline
std::ostream& operator<<(std::ostream& out, FiniteMPO const& x)
{
   for (FiniteMPO::const_iterator I = x.begin(); I != x.end(); ++I)
   {
      out << *I << '\n';
   }
   return out;
}

inline
FiniteMPO::FiniteMPO(GenericMPO const& Other)
   : Data(Other)
{
   //   CHECK(Data.back().Basis2().is_identity())("Finite operator: right basis must be scalar");
   CHECK(Data.back().Basis2().is_regular())("Finite operator: left basis must be regular");
}
