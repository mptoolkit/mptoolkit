
#include <type_traits>
#include <iostream>

template <typename T>
inline
std::enable_if_t<std::is_move_constructible<T>::value && !std::is_lvalue_reference<T>::value, T>
copy(T&& x) // x is a universal reference here, which isn't really what we want
{
   std::cout << "universal version\n";
   return std::move(x);
}

template <typename T>
struct foo {};

template <typename T>
inline
foo<T> copy(foo<T> const& f)
{
   std::cout << "partially specialized\n";
   return f;
}

inline
foo<int> copy(foo<int> const& f)
{
   std::cout << "explicitly specialized\n";
   return f;
}


foo<int> test_foo() { return foo<int>(); }

int main()
{
   copy(test_foo());

   foo<int> f;
   copy(f);

   copy(std::move(f));
}
