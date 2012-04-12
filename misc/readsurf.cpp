
#include <iostream>
#include <map>
#include <complex>
#include <cmath>
#include <boost/lexical_cast.hpp>
#include "linearalgebra/matrix.h"
#include "misc/oogl/mesh.h"

using namespace std;
using namespace LinearAlgebra;

typedef map<int,double> row_type;
typedef map<int, row_type> mat_type;

int main(int argc, char** argv)
{
   int x,y;
   complex<double> c;

   mat_type m;

   if (argc != 3)
   {
      cerr << "usage: readsurf <lattice-size> <theta-delta>\n";
      exit(1);
   }
   
   double Sz = boost::lexical_cast<double>(argv[1]);
   double Th = boost::lexical_cast<double>(argv[2]);

   cin >> x >> y >> c;
   //   double z = exp(log(abs(c.real()))/Sz);
   double z = log(abs(c.real()))/Sz;
   m[x][y] = z;
   m[y][x] = z;
   m[x][x] = 0;
   m[y][y] = 0;

   int umax = std::max(x,y), umin = std::min(x,y);
   double zmin = z, zmax = z;

   while (cin >> x >> y >> c)
   {
      umin = std::min(umin, x);
      umax = std::max(umax, x);
      umin = std::min(umin, y);
      umax = std::max(umax, y);
      //      double z = exp(log(abs(c.real()))/Sz);
      double z = log(abs(c.real()))/Sz;
      zmax = std::max(zmax, z);
      zmin = std::min(zmin, z);
      m[x][y] = z;
      m[y][x] = z;
      m[x][x] = 0;
      m[y][y] = 0;
   }
   //zmax = std::max(zmax, 0.0);

   double ZSize = double(umax - umin + 1);

   mat_type MDeriv;
   double vMin = 0, vMax = 0;
   std::map<int, double> vLastRow;
   double vLastCol;
   for (row_type::const_iterator j = m.begin()->second.begin(); j != m.begin()->second.end(); ++j)
   {
      vLastRow[j->first] = j->second;
   }   
   for (mat_type::const_iterator i = m.begin(); i != m.end(); ++i)
   {
      vLastCol = i->second.begin()->second;
      for (row_type::const_iterator j = i->second.begin(); j != i->second.end(); ++j)
      {
	 double v = std::sqrt((vLastRow[j->first] - j->second) * (vLastRow[j->first] - j->second)
			      + (vLastCol - j->second) * (vLastCol - j->second));
	 MDeriv[i->first][j->first] = v;
	 //	 TRACE(i->first)(j->first)(v);
	 if (vMin == 0)
	    vMin = v;
	 else
	    vMin = std::min(v, vMin);
	 if (vMax == 0)
	    vMax = v;
	 else
	    vMax = std::max(v, vMax);
	 vLastCol = j->second;
	 vLastRow[j->first] = j->second;
      }
   }

   Matrix<double> M(umax-umin+1, umax-umin+1, 1);
   Matrix<oogl::Color> C(umax-umin+1, umax-umin+1);
   for (mat_type::const_iterator i = m.begin(); i != m.end(); ++i)
   {
      for (row_type::const_iterator j = i->second.begin(); j != i->second.end(); ++j)
      {
         M(i->first-umin, j->first-umin) = j->second * ZSize;
	 //	 double Color = double(j->second - zmin) / (zmax - zmin);
	 double Color =(MDeriv[i->first][j->first] - vMin) / (vMax - vMin);
         C(i->first-umin, j->first-umin) = oogl::Color::HSV(359.0 * Color, 1, 1);
      }
   }

   oogl::OoglObject MyMesh(oogl::ColorZMesh(M, C));
   std::cout << MyMesh << std::endl;
}
