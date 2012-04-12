// -*- C++ -*- $Id$

//
// Attribute
//

template <typename T>
Attribute::Attribute(T const& x)
   : value_(boost::lexical_cast<std::string>(x))
{
}

template <typename T>
Attribute& Attribute::operator=(T const& x)
{
   value_ = boost::lexical_cast<std::string>(x);
   return *this;
}

template <typename T>
Attribute& Attribute::operator+=(T const& x)
{
   value_ += boost::lexical_cast<std::string>(x);
   return *this;
}

#if 0
template <typename T>
Attribute::operator T() const
{
   return boost::lexical_cast<T>(value_);
}
#endif

template <typename T>
T Attribute::as() const
{
   return boost::lexical_cast<T>(value_);
}

template <typename T>
T Attribute::get_or_default(T const& Def) const
{
   if (value_.empty())
      return Def;
   return boost::lexical_cast<T>(value_);
}
