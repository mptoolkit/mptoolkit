// ENDHEADER

#if !defined(MPTOOLKIT_QUANTUMNUMBERS_QREP_H)
#define MPTOOLKIT_QUANTUMNUMBERS_QREP_H

#include "quantumnumbers/detail/traits.h"

namespace QuantumNumbers
{

class symmetry_not_integral : public std::runtime_error
{
   public:
      symmetry_not_integral(std::string const& Name)
         : std::runtime_error("symmetry is not integral: " + Name) {}
};

template <typename Symmetry>
inline
typename std::enable_if<detail::has_func_size<Symmetry>::value, int>::type
size(Symmetry const& s)
{
   return s.size();
}

template <typename Symmetry>
inline
typename std::enable_if<detail::has_func_enumerate<Symmetry>::value,
                        typename Symmetry::value_type>::type
enumerate(Symmetry const& s, int n)
{
   return s.enumerate(n);
}

template <typename Symmetry>
inline
typename std::enable_if<detail::has_func_name<Symmetry>::value, std::string>::type
name(Symmetry const& s)
{
   return s.name();
}

template <typename Symmetry>
inline
typename std::enable_if<!detail::has_func_name<Symmetry>::value
                    && detail::has_func_static_name<Symmetry>::value, std::string>::type
name(Symmetry const& s)
{
   return Symmetry::static_name();
}

template <typename Symmetry>
inline
typename std::enable_if<detail::has_func_qdim<Symmetry>::value, real>::type
qdim(Symmetry const& s, typename Symmetry::value_type n)
{
   static_assert(std::is_same<decltype(s.qdim(n)), real>::value);
   return s.qdim(n);
}

template <typename Symmetry>
inline
typename std::enable_if<!detail::has_func_qdim<Symmetry>::value
    && detail::has_func_degree<Symmetry>::value, real>::type
qdim(Symmetry const& s, typename Symmetry::value_type n)
{
   return s.degree(n);
}

template <typename Symmetry>
inline
typename std::enable_if<detail::has_func_degree<Symmetry>::value, int>::type
degree(Symmetry const& s, typename Symmetry::value_type n)
{
   return s.degree(n);
}

template <typename Symmetry>
inline
typename std::enable_if<detail::has_func_qdim<Symmetry>::value
                    && !detail::has_func_degree<Symmetry>::value, int>::type
degree(Symmetry const& s, typename Symmetry::value_type n)
{
   throw symmetry_not_integral(name(s));
}

template <typename Symmetry>
inline
typename Symmetry::value_type
adjoint(Symmetry const& s, typename Symmetry::value_type n)
{
   return s.adjoint(n);
}

template <typename Symmetry>
inline
bool
is_scalar(Symmetry const& s, typename Symmetry::value_type n)
{
   return s.is_scalar(n);
}

template <typename Symmetry>
inline
typename std::enable_if<detail::has_func_fusion_22<Symmetry>::value,
                        real>::type
fusion_31(Symmetry const& s, typename Symmetry::value_type a,
          typename Symmetry::value_type b,
          typename Symmetry::value_type c,
          typename Symmetry::value_type d,
          typename Symmetry::value_type ab,
          typename Symmetry::value_type bc)
{
   static_assert(std::is_same<decltype(s.fusion_31(a,b,c,d,ab,bc)), real>::value);
   return s.fusion_31(a,b,c,d,ab,cd);
}

template <typename Symmetry>
inline
typename std::enable_if<!detail::has_func_fusion_22<Symmetry>::value
&& detail::has_func_coupling_6j<Symmetry>::value,
                        real>::type
fusion_31(Symmetry const& s,
          typename Symmetry::value_type a,
          typename Symmetry::value_type b,
          typename Symmetry::value_type c,
          typename Symmetry::value_type d,
          typename Symmetry::value_type ab,
          typename Symmetry::value_type cd)
{
   static_assert(std::is_same<decltype(s.coupling_6j(a,b,c,d,e,f)), real>::value);
   // fusion in terms of coupling_6j
}

template <typename Symmetry>
inline
typename std::enable_if<detail::has_func_fusion_22<Symmetry>::value,
                        real>::type
fusion_22(Symmetry const& s,
          typename Symmetry::value_type a,
          typename Symmetry::value_type b,
          typename Symmetry::value_type c,
          typename Symmetry::value_type d,
          typename Symmetry::value_type e,
          typename Symmetry::value_type f)
{
   static_assert(std::is_same<decltype(s.fusion_22(a,b,c,d,e,f)), real>::value);
   return s.fusion_22(a,b,c,d,e,f);
}

template <typename Symmetry>
inline
typename std::enable_if<!detail::has_func_fusion_22<Symmetry>::value
&& detail::has_func_coupling_6j<Symmetry>::value,
                        real>::type
fusion_22(Symmetry const& s,
          typename Symmetry::value_type a,
          typename Symmetry::value_type b,
          typename Symmetry::value_type c,
          typename Symmetry::value_type d,
          typename Symmetry::value_type e,
          typename Symmetry::value_type f)
{
   static_assert(std::is_same<decltype(s.coupling_6j(a,b,c,d,e,f)), real>::value);
   // fusion in terms of coupling_6j
}

//
// coupling_9j
//

template <typename Symmetry>
inline
typename std::enable_if<detail::has_func_coupling_9j<Symmetry>::value,
                        real>::type
coupling_9j(Symmetry const& s,
            typename Symmetry::value_type a,
            typename Symmetry::value_type b,
            typename Symmetry::value_type c,
            typename Symmetry::value_type d,
            typename Symmetry::value_type e,
            typename Symmetry::value_type f,
            typename Symmetry::value_type g,
            typename Symmetry::value_type h,
            typename Symmetry::value_type i)
{
   static_assert(std::is_same<decltype(s.coupling_6j(a,b,c,d,e,f)), real>::value);
   return s.coupling_9j(a,b,c,d,e,f,g,h,i);
}


template <typename Symmetry>
inline
typename std::enable_if<!detail::has_func_coupling_9j<Symmetry>::value
&& detail::has_func_coupling_6j<Symmetry>::value,
                        real>::type
coupling_9j(Symmetry const& s,
            typename Symmetry::value_type a,
            typename Symmetry::value_type b,
            typename Symmetry::value_type c,
            typename Symmetry::value_type d,
            typename Symmetry::value_type e,
            typename Symmetry::value_type f,
            typename Symmetry::value_type g,
            typename Symmetry::value_type h,
            typename Symmetry::value_type i)
{
   static_assert(std::is_same<decltype(s.coupling_6j(a,b,c,d,e,f)), real>::value);
   // 9j in terms of coupling_6j
}


}; // namespace QuantumNumbers


#endif
