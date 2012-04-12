// -*- C++ -*- $Id$

#include <iostream>
#include <cmath>
#include "common/math_const.h"
#include <boost/lexical_cast.hpp>

using namespace std;

int main(int argc, char** argv)
{
   if (argc != 5)
   {
      cerr << "usage: makeparams-siam <Bandwidth> <Gamma> <epsilonD> <N>\n";
      return 1;
   }

   double Bandwidth = boost::lexical_cast<double>(argv[1]);
   double Gamma = boost::lexical_cast<double>(argv[2]);
   double EpsilonD = boost::lexical_cast<double>(argv[3]);
   int N = boost::lexical_cast<int>(argv[4]);
   double Hop = Bandwidth / 2.0;
   double Epsilon = 0;

   std::cout.precision(20);
   std::cout << (EpsilonD) << '\n';
   std::cout << Gamma << ' ' << Epsilon << '\n';
   for (int i = 1; i < N-1; ++i)
   {
      std::cout << Hop << ' ' << Epsilon << '\n';
   }
}
