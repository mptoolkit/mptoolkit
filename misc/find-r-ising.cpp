
#include <map>
#include <vector>
#include <iostream>
#include <tuple>
#include <complex>

using complex = std::complex<double>;

std::vector<int> QNums = {0,1,2};

// The known R-matrix elements R^{ab}_c
std::map<std::tuple<int,int,int>, complex> R = {{{0,0,0},1},
						{{0,1,1},1},
						{{1,0,1},1},
						{{0,2,1},1},
						{{2,0,1},1}};

// The F-symbols
std::map<std::tuple<int,int,int,int,int,int>, complex>
F = {{{0,0,0,0,0,0},1},
     {{0,0,1,1,1,0},1},
     {{0,0,2,2,2,0},1},
     {{0,1,0,1,1,1},1},
     {{0,2,0,2,2,2},1},
     {{1,0,0,1,0,1},1},
     {{2,0,0,2,0,2},1},
     {{1,0,1,0,1,1},1},
     {{1,1,0,0,1,0},1},
     {{1,1,1,1,0,0},1}};



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
		  complex LHFactor = F[{x2,x3,x1,x4,c,b}] * 
		     F[{x1,x2,x3,x4,b,a}];
		  std::cout << "c = " << b << " left-hand-factor "
			    << LHFactor << '\n';
		  std::cout << "if there is no multiplicity, then "
			    << "R^{" << x4 << "}_{" << x1 << ","
			    << b << "} = " << (RHFactor / LHFactor)
			    << '\n';
		  if (R.find({x1,b,x4}) != R.end())
		  {
		     if (std::abs(R[{x1,b,x4}] - (RHFactor / LHFactor))
			 > 1e-10)
		     {
			std::cout << "***Mismatch!!!***\n";
		     }
		  }
	       }
	    }
	    std::cout << '\n';
	 }
      }
   }
}
