
#include <map>
#include <vector>
#include <iostream>
#include <tuple>
#include <complex>
#include <cmath>

using complex = std::complex<double>;

std::vector<int> QNums = {0,1};

complex rtt = std::exp(-complex(0,1)*M_PI*4.0/5.0);

// The known R-matrix elements R^{ab}_c
std::map<std::tuple<int,int,int>, complex> R = {{{0,0,0},1},
						{{0,1,1},1},
						{{1,0,1},1},
						{{1,1,0},rtt}};

const complex p = (1.0 + std::sqrt(5.0)) / 2.0;
const complex p2 = std::sqrt(1.0/p);

// The F-symbols
std::map<std::tuple<int,int,int,int,int,int>, complex>
F = {{{0,0,0,0,0,0},1},

     // 1-particle
     {{0,0,1,1,1,0},1},
     {{0,1,0,1,1,1},1},
     {{1,0,0,1,0,1},1},

     // 2-particle
     {{0,1,1,1,0,1},1},
     {{1,0,1,0,1,1},1},
     {{1,1,0,0,1,0},1},

     // 3-particle, off-diagonal
     {{1,1,1,1,0,0},1.0/p},
     {{1,1,1,1,0,1},p2},
     {{1,1,1,1,1,0},p2},
     {{1,1,1,1,1,1},-1.0/p}};



int main()
{
   // loop through pairs of known R-matrix elements and
   // use this as the right-hand-side of the hexagon equation
   // to find any missing R-matrix elements.
   for (auto i : R)
   {
      int x1 = std::get<0>(i.first);
      int x3 = std::get<1>(i.first);
      int c = std::get<2>(i.first);
      for (auto j : R)
      {
	 if (std::get<0>(j.first) != x1)
	    continue;

	 int x2 = std::get<1>(j.first);
	 int a = std::get<2>(j.first);

	 // find all possible F-symbols (F^4_{213})^c_a
	 for (auto f : F)
	 {
	    if (std::get<0>(f.first) != x2)
	       continue;
	    if (std::get<1>(f.first) != x1)
	       continue;
	    if (std::get<2>(f.first) != x3)
	       continue;
	    if (std::get<4>(f.first) != c)
	       continue;
	    if (std::get<5>(f.first) != a)
	       continue;

	    int x4 = std::get<3>(f.first);

	    complex RHFactor =  F[{x2,x1,x3,x4,c,a}] * i.second * j.second;

	    std::cout << "we have a candidate for R^a_{bc} with a="
		      << x4 << ", b = " << x1
		      << " with RHS = " << RHFactor << '\n';
	    //std::cout << "x1=" << x1 << ", x2=" << x2 << ", x3=" << x3
	    //      << ", x4=" << x4 << ", a=" << a << ", c=" << c << '\n';
	    // now iterate over possible values of b
	    for (int b : QNums)
	    {
	       if (F.find({x2,x3,x1,x4,c,b}) != F.end()
		   && F.find({x1,x2,x3,x4,b,a}) != F.end())
	       {
		  std::cout << "c=" << b << '\n';
		  complex LHFactor = F[{x2,x3,x1,x4,c,b}] * 
		     F[{x1,x2,x3,x4,b,a}];
		  if (R.find({x1,b,x4}) != R.end())
		  {
		     if (std::abs(R[{x1,b,x4}] - (RHFactor / LHFactor))
			 > 1e-10)
		     {
			std::cout << "***Mismatch!!!***\n";
			std::cout << "We expect "
				  << "R^{" << x4 << "}_{" << x1 << ","
				  << b << "} = " << (RHFactor / LHFactor)
				  << '\n';
			std::cout << "but the value is " << R[{x1,b,x4}] << '\n';
		     }
		     else
		     {
			std::cout << "consistent.\n";
		     }
		  }
		  else
		  {
		  std::cout << "c = " << b << " left-hand-factor "
			    << LHFactor << '\n';
		  std::cout << "if there is no multiplicity, then "
			    << "R^{" << x4 << "}_{" << x1 << ","
			    << b << "} = " << (RHFactor / LHFactor)
			    << '\n';
		  }
	       }
	    }
	    std::cout << '\n';
	 }
      }
   }
}
