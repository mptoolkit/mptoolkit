
#include <iostream>
#include <iomanip>

using namespace std;

template <typename T>
struct base {};

template <typename T>
struct derived : public base<T> {};

typedef char yes[2];
typedef char no[1];

no& foo(...);

template <typename T>
yes& foo(base<T> const& b);

int main()
{
  int i;
  base<long> b;
  derived<double> d;
  cout << sizeof(foo(i)) << '\n' << sizeof(foo(b)) << '\n' << sizeof(foo(d)) << '\n';
}
