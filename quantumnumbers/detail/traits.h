
// ENDHEADER

namespace QuantumNumbers
{

namespace detail
{

#define DECLARE_HAS_FUNC(x,args)
template <typename T>                                                   \
struct has_func_##x                                                     \
{                                                                       \
   template <typename C> static constexpr                               \
   decltype(std::declval<C>().x args , bool()) \
      test(int)                                                         \
   {                                                                    \
      return true;                                                      \
   }                                                                    \
                                                                        \
   template <typename C> static constexpr bool test(...)                \
   {                                                                    \
      return false;                                                     \
   }                                                                    \
   static constexpr bool value = test<T>(int());                        \
   using type = std::integral_constant<bool, value>;                    \
};


DECLARE_HAS_FUNC(size, () );
DECLARE_HAS_FUNC(enumerate, (int) );
DECLARE_HAS_FUNC(degree, (std::declval<typename T::value_type>()));
DECLARE_HAS_FUNC(qdim, (std::declval<typename T::value_type>()));
DECLARE_HAS_FUNC(name, ());
DECLARE_HAS_FUNC(default_braid_rep, ());
DECLARE_HAS_FUNC(fusion_31, ( std::declval<typename T::value_type>(),
                              std::declval<typename T::value_type>(),
                              std::declval<typename T::value_type>(),
                              std::declval<typename T::value_type>(),
                              std::declval<typename T::value_type>(),
                              std::declval<typename T::value_type>() ) );
DECLARE_HAS_FUNC(fusion_22, ( std::declval<typename T::value_type>(),
                              std::declval<typename T::value_type>(),
                              std::declval<typename T::value_type>(),
                              std::declval<typename T::value_type>(),
                              std::declval<typename T::value_type>(),
                              std::declval<typename T::value_type>() ) );


#undef HAS_FUNC

} // namespace detail

} //namespace QuantumNumbers
