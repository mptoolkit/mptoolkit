
#include <boost/lexical_cast.hpp>
#include <iostream>

using namespace std;

int main(int argc, char** argv)
{
   if (argc != 5)
   {
      std::cerr << "usage: make-tj-op <L> <HoppingRange> <SpinRange> <CoulombRange>\n";
      exit(1);
   }

   int L = boost::lexical_cast<int>(argv[1]);
   int HopR = boost::lexical_cast<int>(argv[2]);
   int SpinR = boost::lexical_cast<int>(argv[3]);
   int CoulombR = boost::lexical_cast<int>(argv[4]);

   for (int i = 1; i <= L; ++i)
   {
      for (int j = i+1; j <= L; ++j)
      {
         if (j-i <= HopR)
            std::cout <<  " PairHopg("<<i<<","<<j<<")";
         if (j-i <= SpinR)
            std::cout <<  " SpinHop("<<i<<","<<j<<")";
         if (j-i <= CoulombR)
            std::cout <<  " VgHop("<<i<<","<<j<<")";
      }
   }
}
