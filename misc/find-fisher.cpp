// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// misc/find-fisher.cpp
//
// Copyright (C) 2022 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#include "common/polynomial.h"
#include "misc/simple-matrix.h"
#include <complex>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <set>

using Poly = Polynomial<std::complex<double>>;

double timestep;
double betastep;

// The matrix of data starts from (t0,b0) at the top left (0,0) entry and
// advances in time by t0 for each row, and advances by imaginary time b0 for each
// column.
matrix<Poly> Rate;

std::vector<double> tGrid, bGrid;
int tGridSize, bGridSize;

int RowOf(double t)
{
   auto I = std::lower_bound(tGrid.begin(), tGrid.end(), t);
   return I-tGrid.begin();
}

int ColOf(double b)
{
   auto I = std::lower_bound(bGrid.begin(), bGrid.end(), b);
   return I-bGrid.begin();
}

double normalize_angle(double x)
{
   int sign = 1;
   if (x < 0)
   {
      sign = -1;
      x = -x;
   }
   x = std::fmod(x, M_PI*2);
   if (x > M_PI)
      x = x - 2*M_PI;
   if (sign == -1)
      x = -x;
   return x;
}

// Information structure about a Fisher line that passes through a link in the lattice.
// We assume that information about the grid location is already available; we only
// store information about the line itself here.
struct FisherInfo
{
   static int NumTags;

   enum StatusType { Empty, Zero, NoZero };
   StatusType Status;

   // -z at the location of the zero
   std::complex<double> zm;
   double Real;
   double ImagDiff;
   int tag; // tag is initialized to zero and can be used by callers, eg to
   // identify specific fisher lines
   double DerivMod;

   FisherInfo() : Status(Empty), tag(-1) {}
   FisherInfo(std::complex<double> zm_, double Real_, double ImagDiff_, double DerivMod_)
      : Status(Zero), zm(zm_), Real(Real_), ImagDiff(ImagDiff_), DerivMod(DerivMod_), tag(NumTags++) {}
   explicit FisherInfo(bool b) : Status(NoZero), tag(-1) {}

   bool HasZero() const { if (Status == Empty) {std::cerr << "fatal: Fisher info is not intialized.\n"; std::abort(); } return Status == Zero; }
};

std::ostream& operator<<(std::ostream& out, FisherInfo const& f)
{
   if (f.Status == FisherInfo::Empty)
   {
      out << "empty!";
      return out;
   }
   if (f.Status == FisherInfo::NoZero)
   {
      out << "no zero";
      return out;
   }
   out << f.zm << ' ' << f.Real << ' ' << f.ImagDiff << ' ' << f.tag;
   return out;
}

int FisherInfo::NumTags = 0;

// the FisherInfo lives on the bonds of the lattice
matrix<FisherInfo> HorizBonds;
matrix<FisherInfo> VertBonds;

std::map<int,int> ClusterSize; // size of the cluster of Fisher lines as a function of tag
// this is so we can remove single-site clusters, which are likely to be outliers
int const MinCluster = 1; // set to 0 to include all clusters

// detects a DQPT on the line between (time,beta) points (T0,B0) and (T1,B1)
FisherInfo
ProbeLine(int T0, int B0, int T1, int B1)
{
   // -z for the 0 and 1 points
   std::complex<double> z0m = {-bGrid[B0], -tGrid[T0]};
   std::complex<double> z1m = {-bGrid[B1], -tGrid[T1]};

   auto zlm = z1m-z0m;

   Poly f0 = Rate(T0, B0);
   Poly f1 = Rate(T1, B1);

   // For a crossover we need to have a crossing in (f0-f1).
   // This could be an accident though, if there is no DQPT but the polynomials
   // happen to cross at some accuracy epsilon.  So the seond criteria is that
   // there should be a finite discontinuity in imaginary part.
   double const ImagEps = 0.06;
   // The rate function along the line from 0 to zl as x goes from 0 to 1 is
   // f0(x*zl) - f1((x-1.0)*zl)
   auto ratediff = [f0,f1,zlm](double x) { return f0(x*zlm) - f1((x-1.0)*zlm); };
   double r0 = ratediff(0.0).real();
   double r1 = ratediff(1.0).real();
   double const eps = 1E-12;

   #if 0
   std::cout << "f0 is " << f0 << " f1 is " << f1 << '\n';
   std::cout << "z0 is " << z0 << " z1 is " << z1 << '\n';
   std::cout << "f0(z0) is " << f0(0.0) << " f1(z0) is " << f1(-zl) << '\n';
   std::cout << "f1(z1) is " << f1(0.0) << " f0(z1) is " << f0(zl) << '\n';
   #endif

   if (r0 * r1 >= 0.0)
      return FisherInfo(false);

   // we have a candidate.  Do a binary search for the location
   double x0 = 0.0;
   double x1 = 1.0;
   double x = 0.5*(x0+x1);
   while (x1-x0 > eps)
   {
      if (r0 * ratediff(x).real() > 0)
         x0 = x;
      else
         x1 = x;
      x = 0.5*(x0+x1);
   }
   if (std::abs(normalize_angle(ratediff(x).imag())) < ImagEps)
      return FisherInfo(false);

   #if 0
   std::cout << "found a crossing at " << -(z0m + x*zlm) << " wih x = " << x << '\n';
   std::cout << "T0=" << T0 << " B0=" << B0 << " T1=" << T1 << " B1=" << B1 << '\n';
   std::cout << "f0 is " << f0(x*zlm) << '\n';
   std::cout << "f1 is " << f1((x-1.0)*zlm) << '\n';
   std::cout << "imag difference is " << std::remainder(ratediff(x).imag(),2*M_PI) << '\n';
   #endif

   #if 0
   std::cout << "T0=" << T0 << " B0=" << B0 << " T1=" << T1 << " B1=" << B1 << ' ';
   std::cout << "f0 is " << f0(x*zlm) << ' ';
   std::cout << "f1 is " << f1((x-1.0)*zlm) << '\n';
   #endif

   return FisherInfo(z0m + x*zlm, f0(x*zlm).real(), normalize_angle(ratediff(x).imag()),
      std::abs(f0.derivative(x*zlm)));
}

// This function walks the HorizBonds and VertBonds arrays of fisher discoontinuities
// and tags them as distinct fisher lines, by setting the tag element of the FisherInfo
// to a distinct number for each line.  Fisher lines that intersect are regarded as distinct,
// and each line that intersects at a plaquette is given a different tag.
void IdentifyFisherLines()
{
   // On entry, each identified discontinuity has a different tag.  Merge the tags
   // using the Hoshen–Kopelman algorithm.
   // The general strategy is to iterate through each plaquette in the lattice, and
   // if there are exactly two tags on the edges of the plaquette then merge them.
   std::vector<int> CanonicalTag(FisherInfo::NumTags);
   int n = 0;
   for (auto& t : CanonicalTag)
   {
      t = n++;
   }

   auto FindTag = [&CanonicalTag] (int t) { int tnew = t; while (CanonicalTag[tnew] != tnew) { tnew = CanonicalTag[tnew];} while (CanonicalTag[t] != t) { int z = CanonicalTag[t]; CanonicalTag[t] = tnew; t = z; } return tnew; };

   auto MergeTag = [&CanonicalTag,FindTag] (int t1, int t2) { CanonicalTag[FindTag(t1)] = FindTag(t2); };

   // iterate over all plaquettes
   for (int i = 0; i < int(tGrid.size()-1); ++i)
   {
      for (int j = 0; j < int(bGrid.size()-1); ++j)
      {
         FisherInfo& Left = VertBonds(i,j);
         FisherInfo& Right = VertBonds(i,j+1);
         FisherInfo& Down = HorizBonds(i,j);
         FisherInfo& Up = HorizBonds(i+1,j);

         int NumLines = 0;
         if (Left.tag >= 0)
            ++NumLines;
         if (Right.tag >= 0)
            ++NumLines;
         if (Up.tag >= 0)
            ++NumLines;
         if (Down.tag >= 0)
            ++NumLines;

         if (NumLines == 2)
         {
            int tag1 = std::max(std::max(Left.tag, Right.tag), std::max(Up.tag, Down.tag));
            if (Left.tag >= 0 && Left.tag != tag1)
               MergeTag(Left.tag, tag1);
            else if (Right.tag >= 0 && Right.tag != tag1)
               MergeTag(Right.tag, tag1);
            else if (Up.tag >= 0 && Up.tag != tag1)
               MergeTag(Up.tag, tag1);
            else if (Down.tag >= 0 && Down.tag != tag1)
               MergeTag(Down.tag, tag1);
         }
      }
   }

   // Now clean up the labels
   for (int i = 0; i < VertBonds.rows(); ++i)
   {
      for (int j = 0; j < VertBonds.cols(); ++j)
      {
         if (VertBonds(i,j).HasZero())
         {
            VertBonds(i,j).tag = FindTag(VertBonds(i,j).tag);
            ++ClusterSize[VertBonds(i,j).tag];
            //std::cout << "found tag " << VertBonds(i,j).tag << '\n';
         }
      }
   }
   for (int i = 0; i < HorizBonds.rows(); ++i)
   {
      for (int j = 0; j < HorizBonds.cols(); ++j)
      {
         if (HorizBonds(i,j).HasZero())
         {
            HorizBonds(i,j).tag = FindTag(HorizBonds(i,j).tag);
            ++ClusterSize[HorizBonds(i,j).tag];
            //std::cout << "found tag " << HorizBonds(i,j).tag << '\n';
         }
      }
   }
}

enum class Direction { Left, Right, Up, Down };

std::ostream& operator<<(std::ostream& out, Direction d)
{
   if (d == Direction::Left)
      out << "left";
   if (d == Direction::Right)
      out << "right";
   if (d == Direction::Up)
      out << "up";
   if (d == Direction::Down)
      out << "down";
   return out;
}

// Here T,B are interpreted as a plaquette.
// Attempt to move to the plaquette in the direction d,
// and update d to point to the new direction of the Fisher line at that plaquette.
// If the move isn't possible because we are already at the end of a Fisher line
// (because we are at the end of the grid, or the next plaquette doesn't have an
// outgoing Fisher line), then return false.
bool MoveFisherPlaquette(int& T, int& B, Direction& d, int tag)
{
   int Tnew = T;
   int Bnew = B;
   if (d == Direction::Up)
   {
      assert(HorizBonds(T+1,B).tag == tag);
      if (T == tGridSize-2) return false;
      ++Tnew;
   }
   if (d == Direction::Down)
   {
      assert(HorizBonds(T,B).tag == tag);
      if (T == 0) return false;
      --Tnew;
   }
   if (d == Direction::Left)
   {
      assert(VertBonds(T,B).tag == tag);
      if (B == 0) return false;
      --Bnew;
   }
   if (d == Direction::Right)
   {
      assert(VertBonds(T,B+1).tag == tag);
      if (B == bGridSize-2) return false;
      ++Bnew;
   }

   // These are flipped from what we expect, since it is the direction we are
   // coming from, not moving to
   if (d != Direction::Right && VertBonds(Tnew,Bnew).tag == tag)
      d = Direction::Left;
   else if (d != Direction::Left && VertBonds(Tnew,Bnew+1).tag == tag)
      d = Direction::Right;
   else if (d != Direction::Up && HorizBonds(Tnew,Bnew).tag == tag)
      d = Direction::Down;
   else if (d != Direction::Down && HorizBonds(Tnew+1,Bnew).tag == tag)
      d = Direction::Up;
   else
      return false;
   T = Tnew;
   B = Bnew;
   return true;
}

// We identify the end point of a Fisher line by specifing plaquette
// and a direction which specifies the bond.
// We can be 'one past the end' of the grid of plaquettes, as long as we
// are then referring to the left or down bonds as required.
// Once we get to the end point, we 'turn around' by moving to the next
// plaquette and flipping the direction.  This means that we might end up
// on a plaquette that is one past the end, on any side of the grid.
void FindFisherEndPoint(int& T, int& B, Direction& d, int tag)
{
   // can we move plaquette?  If not, we're already at an end point
   while (MoveFisherPlaquette(T,B,d,tag))
   {
      // do nothing
   }
   // Now flip direction
   if (d == Direction::Left)
   {
      --B;
      d = Direction::Right;
   }
   else if (d == Direction::Right)
   {
      ++B;
      d = Direction::Left;
   }
   else if (d == Direction::Up)
   {
      ++T;
      d = Direction::Down;
   }
   else if (d == Direction::Down)
   {
      --T;
      d = Direction::Up;
   }
}

// Walk a Fisher line that starts at plaquette (T,B) direction d and output the coordinates
void ShowWalkFisherLine(int T, int B, Direction d, int tag)
{
   // (T,B) is a plaquette location
   bool OK = true;
   while (OK)
   {
      FisherInfo f;
      if (d == Direction::Left)
         f = VertBonds(T,B);
      else if (d == Direction::Right)
         f = VertBonds(T,B+1);
      else if (d == Direction::Up)
         f = HorizBonds(T+1,B);
      else if (d == Direction::Down)
         f = HorizBonds(T,B);

      if (f.tag < 0)
      {
         std::cerr << "error! tag is not defined at time " << tGrid[T] << " beta " << bGrid[B]
            << " direction " << d << '\n';
         std::abort();
      }

      std::cout << T << ' ' << B << ' ' << d << ' ';
      if (d == Direction::Up || d == Direction::Left)
      {
         std::cout << -f.zm.real() << ' ' << -f.zm.imag() << ' ' << f.Real << ' ' << -f.ImagDiff
            << ' ' << f.DerivMod << ' ' << tag << '\n';
      }
      else
      {
         std::cout << -f.zm.real() << ' ' << -f.zm.imag() << ' ' << f.Real << ' ' << f.ImagDiff
            << ' ' << f.DerivMod << ' ' << tag << '\n';
      }

      #if 0
      if (std::abs(f.zm.real()) < 1E-8)
      {
         std::cout << "error at " << T << ' ' << B << ' ' << d << << tGrid[T] << ' ' << bGrid[B] << '\n';
         std::cout << f << '\n';
      }
      #endif
      OK = MoveFisherPlaquette(T, B, d, tag);
   }
}

void ShowAllFisherLines()
{
   // go through each Fisher line and identify the start point, and display it
   std::set<int> ShownTags;

   for (int i = 0; i < VertBonds.rows(); ++i)
   {
      for (int j = 0; j < VertBonds.cols(); ++j)
      {
         int tag = VertBonds(i,j).tag;
         if (tag >= 0 && ClusterSize[tag] > MinCluster)
         {
            if (ShownTags.find(tag) == ShownTags.end())
            {
               ShownTags.insert(tag);
               int T = i;
               int B = j;
               Direction d = Direction::Left;
               //std::cerr << "Searching for Fisher line starting " << T << ' ' << B << ' ' << d << ' ' << tag << '\n';
               FindFisherEndPoint(T,B, d, VertBonds(i,j).tag);
               //std::cerr << "Fisher end point is " << T << ' ' << B << ' ' << d << '\n';
               ShowWalkFisherLine(T,B, d, VertBonds(i,j).tag);
               std::cout << '\n' << '\n';
            }
         }
      }
   }
   for (int i = 0; i < HorizBonds.rows(); ++i)
   {
      for (int j = 0; j < HorizBonds.cols(); ++j)
      {
         int tag = HorizBonds(i,j).tag;
         if (tag >= 0 && ClusterSize[tag] > MinCluster)
         {
            if (ShownTags.find(tag) == ShownTags.end())
            {
               ShownTags.insert(tag);
               int T = i;
               int B = j;
               Direction d = Direction::Down;
               //std::cerr << "Searching for Fisher line starting " << T << ' ' << B << ' ' << d << ' ' << tag << '\n';
               FindFisherEndPoint(T,B, d, HorizBonds(i,j).tag);
               //std::cerr << "Fisher end point is " << T << ' ' << B << ' ' << d << '\n';
               ShowWalkFisherLine(T,B, d, HorizBonds(i,j).tag);
               std::cout << '\n' << '\n';
            }
         }
      }
   }
}

int main()
{
   std::map<std::tuple<double,double>, Poly> AllData;
   // Read data files from standard input
   std::string line;
   while (std::getline(std::cin, line))
   {
      if (line.empty() || line[0] == '#')
         continue;

      std::istringstream In(line);
      double t, b;
      double r0r, r0c, r1r, r1c, r2r, r2c;
      In >> t >> b >> r0r >> r0c >> r1r >> r1c >> r2r >> r2c;

      double r, c;
      Poly p;
      p[0] = {r0r, r0c};
      int deg = 1;
      while (In >> r >> c)
      {
         p[deg] = {r,c};
         ++deg;
      }
      AllData.insert({{t,b}, p});
   }
   if (AllData.empty())
   {
      std::cerr << "No data!\n";
      return 1;
   }
   // find the dimensions of the grid
   std::set<double> T, B;
   for (auto const& x : AllData)
   {
      T.insert(std::get<0>(x.first));
      B.insert(std::get<1>(x.first));
   }
   tGrid = std::vector<double>(T.begin(), T.end());
   tGridSize = tGrid.size();
   bGrid = std::vector<double>(B.begin(), B.end());
   bGridSize = bGrid.size();
   std::cerr << "Time grid is size " << tGrid.size() << " starting " << tGrid[0]
      << " finishing " << tGrid.back() << '\n';
   std::cerr << "Beta grid is size " << bGrid.size() << " starting " << bGrid[0]
      << " finishing " << bGrid.back() << '\n';

   Rate = matrix<Poly>(tGrid.size(), bGrid.size());
   for (auto const& x : AllData)
   {
      int r = RowOf(std::get<0>(x.first));
      int c = ColOf(std::get<1>(x.first));
      Rate(r,c) = x.second;
   }

   // Identify the fisher lines as crossings in horizontal and vertical bonds
   // Vertical bonds capture a horizontal crossing Fisher line at some fixed
   // value of B.
   // Horizontal lines capture a vertical crossing Fisher line at some fixed
   // value of T.
   // The number of vertical bonds is tGrid.size()-1 by bGrid.size()
   // the horizontal bonds arae tGrid.size() by bGrid.size()-1
   VertBonds = matrix<FisherInfo>(tGrid.size()-1, bGrid.size());
   for (int i = 0; i < tGrid.size()-1; ++i)
   {
      for (int j = 0; j < bGrid.size(); ++j)
      {
         VertBonds(i,j) = ProbeLine(i,j, i+1,j);
      }
   }

   HorizBonds = matrix<FisherInfo>(tGrid.size(), bGrid.size()-1);
   for (int i = 0; i < tGrid.size(); ++i)
   {
      for (int j = 0; j < bGrid.size()-1; ++j)
      {
         HorizBonds(i,j) = ProbeLine(i,j, i,j+1);
      }
   }

   // Now identify continuous fisher lines.
   IdentifyFisherLines();
   ShowAllFisherLines();
}
