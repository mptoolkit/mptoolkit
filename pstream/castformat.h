// -*- C++ -*- $Id$

#include <boost/cast.hpp>

template <typename DestType, typename SrcType>
inline
DestType format_cast(SrcType f)
{
   return Private::FormatCastHelper<DestType, 
     CurrentFormat::IsTrivialConversion<SrcType, DestType::value>::apply(f);
}

namespace Private
{

template <typename DestType, bool isTrivial>
struct FormatCastHelper
{
   static DestType apply(SrcType f) { return boost::numeric_cast<DestType>(f); }
};

template <typename DestType>
struct FormatCastHelper<DestType, true>
{
   static DestType apply(DestType f) { return static_cast<DestType>(f); }
};

} // namespace Private

