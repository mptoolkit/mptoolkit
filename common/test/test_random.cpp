
#include "common/randutil.h"
#include <iostream>

int main()
{
   // the 10000th consecutive invocation of a default-contructed std::mt19937 is required to produce the value 4123659995.

   for (int i = 0; i <  9999; ++i)
   {
      randutil::u_rand();
   }

   unsigned x = randutil::u_rand();
   std::cout << x << '\n';
}
