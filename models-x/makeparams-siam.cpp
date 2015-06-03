// -*- C++ -*- $Id$

#include <iostream>
#include <cmath>
#include "common/math_const.h"
#include <boost/lexical_cast.hpp>

using math_const::pi;

using namespace std;

double xi(int n, double Lambda)
{
   double Result = (1.0 - pow(Lambda, -n-1))
      / sqrt(1.0 - pow(Lambda, -2*n-1))
      / sqrt(1.0 - pow(Lambda, -2*n-3));
   return Result;
}

int main(int argc, char** argv)
{
   if (argc != 6)
   {
      cerr << "usage: makeparams-siam <Bandwidth> <Gamma> <epsilonD> <Lambda> <N>\n";
      return 1;
   }
   double Bandwidth = boost::lexical_cast<double>(argv[1]);
   double Gamma = boost::lexical_cast<double>(argv[2]);
   double EpsilonD = boost::lexical_cast<double>(argv[3]);
   double Lambda = boost::lexical_cast<double>(argv[4]);
   int N = boost::lexical_cast<int>(argv[5]);

   std::cout.precision(20);
   std::cout << (EpsilonD) << '\n';

   double Epsilon = 0;
   double Hop0 = std::sqrt(2 * Gamma * Bandwidth / pi);
   std::cout << Hop0 << ' ' << Epsilon << '\n';

   for (int i = 0; i < N-1; ++i)
   {
      double Hop = (Bandwidth / 2.0) * (1.0 + (1.0 / Lambda)) * pow(Lambda, -i/2.0) * xi(i, Lambda);
      std::cout << Hop << ' ' << Epsilon << '\n';
   }
}
