// -*- C++ -*-

// inlines

inline
bool operator==(signed_half_int h1, signed_half_int h2)
{
   return h1.twice() == h2.twice();
}

inline
bool operator!=(signed_half_int h1, signed_half_int h2)
{
   return h1.twice() != h2.twice();
}

inline
bool operator<(signed_half_int h1, signed_half_int h2)
{
   return h1.twice() < h2.twice();
}

inline
bool operator>(signed_half_int h1, signed_half_int h2)
{
   return h1.twice() > h2.twice();
}

inline
bool operator<=(signed_half_int h1, signed_half_int h2)
{
   return h1.twice() <= h2.twice();
}

inline
bool operator>=(signed_half_int h1, signed_half_int h2)
{
   return h1.twice() >= h2.twice();
}

inline
half_int abs(signed_half_int h)
{
   return half_int(h.twice_magnitude(), signed_half_int::twice_tag());
}

inline
int sign(signed_half_int h)
{
   return h.sign();
}

#if defined(USE_PSTREAM)
inline
PStream::opstream& operator<<(PStream::opstream& out, const signed_half_int& H)
{
   return out << H.N2;
}

inline
PStream::ipstream& operator>>(PStream::ipstream& in, signed_half_int& H)
{
   return in >> H.N2;
}
#endif // #defined (USE_PSTREAM)
