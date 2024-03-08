
#include <iostream>
#include <cmath>
#include <map>
#include <set>
#include <fstream>
#include <deque>
#include "common/trace.h"

#include "matrix.h"

int constexpr MaxRate = 4;
int constexpr MaxCumulants = 4;

double& re(std::complex<double>& c)
{
    return reinterpret_cast<double (&)[2]>(c)[0];
}

double& im(std::complex<double>& c)
{
    return reinterpret_cast<double (&)[2]>(c)[1];
}

const double& re(const std::complex<double>& c)
{
    return reinterpret_cast<const double (&)[2]>(c)[0];
}

const double& im(const std::complex<double>& c)
{
    return reinterpret_cast<const double (&)[2]>(c)[1];
}

struct GridPoint
{
   double b;
   double t;
   int nrate;
   std::array<std::complex<double>, MaxRate> rate;
   int nc;
   std::array<std::complex<double>, MaxCumulants> c;

   GridPoint() : b(0.0), t(0.0), nrate(0), nc(0) {}

   GridPoint(double b_, double t_, std::complex<double> r) : b(b_), t(t_), nrate(1), rate{{r}}, nc(0) {}

   template <typename Container>
   GridPoint(double b_, double t_, std::complex<double> r, Container c_) : b(b_), t(t_), nrate(1), rate{{r}}, nc(c_.size()) { std::copy(c_.begin(), c_.end(), c.begin()); }

   GridPoint(double b_, double t_, std::complex<double> r, std::complex<double> c1_) : b(b_), t(t_), nrate(1), rate{{r}}, nc(1), c{{c1_}} {}
};

double normalize_angle(double a)
{
   a = fmod(a, 2*M_PI);
   if (a < -M_PI)
      a += M_PI*2;
   if (a > M_PI)
      a -= M_PI*2;
   return a;
}

std::ostream& operator<<(std::ostream& out, GridPoint const& p)
{
   out << p.b << ' ' << p.t << " 0 ";
   for (int i = 0; i < p.nrate; ++i)
   {
      out << p.rate[i].real() << ' ' << (normalize_angle(p.rate[i].imag())*180.0/M_PI) << ' ';
   }
   for (int i = 0; i < p.nc; ++i)
   {
      out << p.c[i].real() << ' ' << p.c[i].imag() << ' ';
   }
   return out;
}

GridPoint Parse(int NumRates, int NumCumulants, std::string const& s)
{
   GridPoint r;
   r.nrate = NumRates;
   r.nc = NumCumulants;
   std::istringstream ss(s);
   double junk;
   ss >> r.b >> r.t >> junk;
   for (int i = 0; i < NumRates; ++i)
   {
      ss >> re(r.rate[i]) >> im(r.rate[i]);
      im(r.rate[i]) = normalize_angle(im(r.rate[i]) * M_PI / 180.0);  // convert to radians
   }
   for (int i = 0; i < NumCumulants; ++i)
      ss >> re(r.c[i]) >> im(r.c[i]);
   if (!ss)
      throw std::runtime_error("bad grid point");

   return r;
}

class MapValues
{
   public:
      MapValues() {}

      MapValues(std::set<double> bmap, std::set<double> tmap);

      bool bexists(double b) const { return b_map.count(b); }
      bool texists(double t) const { return t_map.count(t); }

      int bsize() const { return b_map.size(); }
      int tsize() const { return t_map.size(); }

      double b(int i) const { return bi_map[i]; }
      double t(int i) const { return ti_map[i]; }

      // returns the largest index i such that this->b(i) <= b
      int bm(double b) const;

      // returns the largest index j such that this->t(j) <= t
      int tm(double t) const;

   private:
      std::map<double, int> b_map, t_map;
      std::vector<double> bi_map, ti_map;
};

int MapValues::bm(double b) const
{
   auto i = b_map.find(b);
   if (i != b_map.end())
      return i->second;

   auto j = std::upper_bound(bi_map.begin(), bi_map.end(), b);
   assert(j != bi_map.end());
   return (j-bi_map.begin())-1;
}

int MapValues::tm(double t) const
{
   auto i = t_map.find(t);
   if (i != t_map.end())
      return i->second;

   auto j = std::upper_bound(ti_map.begin(), ti_map.end(), t);
   assert(j != ti_map.end());
   return (j-ti_map.begin())-1;
}

MapValues::MapValues(std::set<double> bs, std::set<double> ts)
{
   int i = 0;
   for (double  b : bs)
   {
      bi_map.push_back(b);
      b_map.insert({b, i++});
   }
   i = 0;
   for (double  t : ts)
   {
      ti_map.push_back(t);
      t_map.insert({t, i++});
   }
}

MapValues MapFromFile(int NumRates, int NumCumulants, std::string f)
{
   std::ifstream in(f);
   std::set<double> bs, ts;
   while (in)
   {
      std::string str;
      std::getline(in, str);
      if (!str.empty() || !(str[0] == '#'))
      {
         try
         {
            GridPoint p = Parse(NumRates, NumCumulants, str);
            bs.insert(p.b);
            ts.insert(p.t);
         }
         catch (...)
         {
            // ignore errors
         }
      }
   }
   return MapValues(bs, ts);
}

// Function to make a fine-grained map that increases the number of grid points by bfactor or tfactor
MapValues MakeInterpolationMap(MapValues const& m, int bfactor, int tfactor)
{
   int bsize = m.bsize();
   int tsize = m.tsize();

   std::set<double> bs;
   bs.insert(m.b(0));
   for (int bi = 0; bi < bsize-1; ++bi)
   {
      double b1 = m.b(bi);
      double b2 = m.b(bi+1);
      for (int j = 1; j < bfactor; ++j)
      {
         bs.insert(b1 + j * ((b2-b1)/bfactor));
      }
      bs.insert(b2);
   }

   std::set<double> ts;
   ts.insert(m.t(0));
   for (int ti = 0; ti < tsize-1; ++ti)
   {
      double t1 = m.t(ti);
      double t2 = m.t(ti+1);
      for (int j = 1; j < tfactor; ++j)
      {
         ts.insert(t1 + j * ((t2-t1)/tfactor));
      }
      ts.insert(t2);
   }

   return MapValues(bs, ts);
}

struct GridType
{
   GridType(MapValues MyMap) : Map(MyMap), Valid(Map.bsize(), Map.tsize(),false), Grid(Map.bsize(), Map.tsize()) {}

   void LoadFromFile(int NumRates, int NumCumulants, std::string f);

   void Interpolate(GridType const& g);

   int bsize() const { return Map.bsize(); }
   int tsize() const { return Map.tsize(); }

   int bm(double b) const { return Map.bm(b); }
   int tm(double t) const { return Map.tm(t); }

   GridPoint operator()(double b, double t) const { return Grid(Map.bm(b), Map.tm(t)); }

   bool is_valid(double b, double t) const { return Valid(Map.bm(b), Map.tm(t)); }

   // Find a set of up to 4 grid points that are nerest to (b,t)
   std::vector<GridPoint> FindGridPointsNearest(double b, double t) const;

   MapValues Map;
   matrix<int> Valid;
   matrix<GridPoint> Grid;
};

void GridType::LoadFromFile(int NumRates, int NumCumulants, std::string f)
{
   std::ifstream in(f);
   while (in)
   {
      std::string str;
      std::getline(in, str);
      if (!str.empty() || !(str[0] == '#'))
      {
         try
         {
            GridPoint p = Parse(NumRates, NumCumulants, str);
            int bi = Map.bm(p.b);
            int ti = Map.tm(p.t);
            Valid(bi,ti) = true;
            Grid(bi,ti) = p;
         }
         catch (...)
         {
            // ignore errors, they will be marked as Valid=false
         }
      }
   }
}

// Find up to 4 points nearest to (b,t) in the grid
std::vector<GridPoint> GridType::FindGridPointsNearest(double b, double t) const
{
   std::vector<GridPoint> Result;
   int bi = Map.bm(b);
   int ti = Map.tm(t);

   //TRACE(b)(t)(bi)(ti)(Map.b(bi))(Map.t(ti))(Valid(bi,ti));

   // probe to the bottom left
   while (bi >= 0 && ti >= 0 && !Valid(bi, ti))
   {
      if ((ti+bi)%2 == 0 && ti > 0)
         --ti;
      else
         --bi;
   }
   //TRACE(bi)(ti);
   if (bi >= 0 && ti >= 0)
   {
      Result.push_back(Grid(bi,ti));
      // early return if (b,t) is actually on a grid point
      if (Map.b(bi) == b && Map.t(ti) == t)
         return Result;
   }

   // probe to the top left
   bi = Map.bm(b);
   ti = Map.tm(t)+1;
   while (ti < Map.tsize() && bi >= 0 && !Valid(bi, ti))
   {
      if ((ti+bi)%2 == 0 && ti < Map.tsize()-1)
         ++ti;
      else
         --bi;
   }
   if (ti < Map.tsize() && bi >= 0)
      Result.push_back(Grid(bi,ti));

   // probe to the bottom right
   bi = Map.bm(b)+1;
   ti = Map.tm(t);
   while (bi < Map.bsize() && ti >= 0 && !Valid(bi,ti))
   {
      if ((ti+bi)%2 == 0 && ti > 0)
         --ti;
      else
         ++bi;
   }
   if (bi < Map.bsize() && ti >= 0)
      Result.push_back(Grid(bi,ti));

   // probe to the top right
   bi = Map.bm(b)+1;
   ti = Map.tm(t)+1;
   while (bi < Map.bsize() && ti < Map.tsize() && !Valid(bi,ti))
   {
      if ((ti+bi)%2 == 0 && ti < Map.tsize()-1)
         ++ti;
      else
         ++bi;
   }
   if (bi < Map.bsize() && ti < Map.tsize())
      Result.push_back(Grid(bi,ti));
   return Result;
}

GridPoint InterpolateFrom(GridPoint x, double b, double t)
{
   assert(x.nc >= 1);
   std::complex<double> dx = {b-x.b, t-x.t};
   std::complex<double> dr = x.c[x.nc-1];
   std::deque<std::complex<double>> c;
   for (int i = x.nc-1; i > 0; --i)
   {
      dr = x.c[i-1] - (dx / (i+1.0)) * dr;
   }
   dr = x.rate[0] - dx*dr;
   c.push_back(x.c[0] - dx*x.c[1]);

   dr = x.rate[0] - dx*x.c[0] + 0.5*dx*dx*x.c[1];
   c.push_back({0.0,0.0});

   return GridPoint(b, t, dr, c);
}

void GridType::Interpolate(GridType const& g)
{
   double b_last = g.Map.b(0);
   double t_last = g.Map.t(0);

   for (int bi = 0; bi < Map.bsize(); ++bi)
   {
      double b = Map.b(bi);
      // if we are on a grid point then update b_last
      if (g.Map.bexists(b))
         b_last = b;

      for (int ti = 0; ti < Map.tsize(); ++ti)
      {
         double t = Map.t(ti);
         // if we are on a grid point then update t_last, but only if it is a valid grid point
         if (g.Map.texists(t) && g.is_valid(b_last, t))
         {
            t_last = t;
         }

         // If we're exactly on a grid point, then don't do an interpolation
         if (b_last == b && t_last == t)
         {
            Grid(bi, ti) = g(b_last, t_last);
            Valid(bi,ti) = true;
            continue;
         }

         auto GridPoints = g.FindGridPointsNearest(b, t);

         for (GridPoint& gg : GridPoints)
         {
            gg = InterpolateFrom(gg, b, t);
         }
         // the interpolation to keep is the one that has the largest (negative) real rate
         int keep = 0;
         for (int i = 1; i < GridPoints.size(); ++i)
         {
            if (GridPoints[i].rate[0].real() > GridPoints[keep].rate[0].real())
               keep = i;
         }

         Grid(bi, ti) = GridPoints[keep];
         Valid(bi, ti) = true;
      }
   }
}

std::ostream& operator<<(std::ostream& out, GridType const& g)
{
   for (int bi = 0; bi < g.Map.bsize(); ++bi)
   {
      for (int ti = 0; ti < g.Map.tsize(); ++ti)
      {
         if (g.Valid(bi, ti))
            out << g.Grid(bi, ti) << '\n';
      }
      out << '\n';
   }
   return out;
}

int main(int argc, char** argv)
{
   if (argc != 3)
   {
      std::cerr << "expected: interpolate <moments-file> <factor>\n";
      return 1;
   }

   int NumRates = 1;
   int NumCumulants = 2;

   std::string File = argv[1];
   int Factor = std::stoi(argv[2]);
   MapValues Map = MapFromFile(NumRates, NumCumulants, File);
   std::cerr << "size is " << Map.bsize() << " by " << Map.tsize() << '\n';

   GridType Grid(Map);
   Grid.LoadFromFile(NumRates, NumCumulants, File);

   MapValues IMap = MakeInterpolationMap(Map, Factor, Factor);
   std::cerr << "new size is " << IMap.bsize() << " by " << IMap.tsize() << '\n';

   GridType IGrid(IMap);
   IGrid.Interpolate(Grid);

   std::cout << IGrid << '\n';
}
