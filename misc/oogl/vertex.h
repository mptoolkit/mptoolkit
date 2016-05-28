// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// misc/oogl/vertex.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Reseach publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER
//
// Classes Vertex, Color, Transform, VertexList, ColorVertexList, Face
// for OOGL
//

#if !defined(VERTEX_H_CHDUY8975849YHP8EW9HP)
#define VERTEX_H_CHDUY8975849YHP8EW9HP

#include <iostream>
#include <vector>
#include <utility>

namespace oogl
{

class Transform;

struct Vertex
{
   double x, y, z;

   Vertex() {}
   Vertex(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

   Vertex& operator+=(Vertex const& v);
   Vertex& operator-=(Vertex const& v);

   Vertex& operator*=(double a);

   Vertex& operator*=(Transform const& t);

   // named constructors for various coordinate systems
   static Vertex cylindrical(double r, double theta, double z);
   static Vertex spherical(double r, double theta, double phi);
   static Vertex cartesian(double x, double y, double z);
};

Vertex operator-(Vertex const& v);

Vertex operator+(Vertex const& v1, Vertex const& v2);
Vertex operator-(Vertex const& v1, Vertex const& v2);

Vertex operator*(double a, Vertex const& v);
Vertex operator*(Vertex const& v, double a);

bool operator==(Vertex const& v1, Vertex const& v2);
bool operator!=(Vertex const& v1, Vertex const& v2);

std::ostream& operator<<(std::ostream& out, Vertex const& v);

// Color is a representation of an RGBA color.  
// There is a named constructor to convert from HSV.
struct Color
{
   double R, G, B, A;

   Color() {}
  
   Color(double R_, double G_, double B_, double A_ = 1.0)
     : R(R_), G(G_), B(B_), A(A_) {}

   // Named ctor to convert from HSV coordinates.
   // H = [0,360], S = [0,1], V = [0,1]
   // S is saturation, 1 = pure (no white), 0 = white
   // V is brightness, 1 = brightest, 0 = black
   static Color HSV(double H, double S, double V, double A = 1.0);

   // Named ctor for RDB data, not needed as this is the default, but for symmetry
   static Color RGB(double R, double G, double B, double A = 1.0)
   { return Color(R,G,B,A); }
};

// The predefined colors are in 
namespace Colors
{

   extern Color const black;
   extern Color const white;
   extern Color const red;
   extern Color const green;
   extern Color const blue;
   extern Color const cyan;
   extern Color const yellow;
   extern Color const magenta;

} // namespace Colors


std::ostream& operator<<(std::ostream& out, Color const& c);


// For some reason, a transform is not an oogl-object.
struct Transform
{
   public:
      Transform() {}

     double operator()(int i, int j) const { return TransformMat[i][j]; }
     double& operator()(int i, int j) { return TransformMat[i][j]; }

     Transform& operator*=(Transform const& t);

     static Transform identity;
     static Transform rot_x(double theta);
     static Transform rot_y(double theta);
     static Transform rot_z(double theta);

     static Transform trans_x(double x);
     static Transform trans_y(double y);
     static Transform trans_z(double z);

     static Transform scale(double s);

     static Transform scale_x(double s);
     static Transform scale_y(double s);
     static Transform scale_z(double s);

     // returns a transformation that shifts by the amount v
     static Transform translate(Vertex const& v);

  private:
     double TransformMat[4][4];

     struct ident_tag {};
     Transform(ident_tag);
};

Vertex operator*(Vertex const& v, Transform const& t);

Transform operator*(Transform const& t1, Transform const& t2);

std::ostream& operator<<(std::ostream& out, Transform const& t);

// VertexList is not an oogl-object
class VertexList
{
   public:
      typedef Vertex                              value_type;
      typedef value_type&                             reference;
      typedef value_type*                             pointer;
      typedef std::vector<value_type>::iterator       iterator;
      typedef std::vector<value_type>::const_iterator const_iterator;

      VertexList() {}

      VertexList(Vertex const& v) { VList.push_back(v); }
      VertexList(Vertex const& v1, Vertex const& v2) 
         { VList.push_back(v1); VList.push_back(v2); }

      int size() const { return VList.size(); }

      void push_back(Vertex const& v) { VList.push_back(v); }
      void pop_back() { VList.pop_back(); }

      iterator begin() { return VList.begin(); }
      iterator end() { return VList.end(); }

      const_iterator begin() const { return VList.begin(); }
      const_iterator end() const { return VList.end(); }

      static char const* prefix();

      VertexList& operator*=(Transform const& t)
      {
	 for (unsigned i = 0; i < VList.size(); ++i)
	 {
	    VList[i] *= t;
	 }
	 return *this;
      }

   private:
      std::vector<value_type> VList;
};

std::ostream& operator<<(std::ostream& out, VertexList const& VList);

// ColorVertexList adds an RGBA Color to the VertexList.  Again it is not an oogl-object
class ColorVertexList
{
   public:
      typedef std::pair<Vertex, Color>                value_type;
      typedef value_type&                             reference;
      typedef value_type*                             pointer;
      typedef std::vector<value_type>::iterator       iterator;
      typedef std::vector<value_type>::const_iterator const_iterator;

      ColorVertexList() {}

      int size() const { return VList.size(); }

      void push_back(Vertex const& v, Color const& c) { VList.push_back(value_type(v, c)); }

      iterator begin() { return VList.begin(); }
      iterator end() { return VList.end(); }

      const_iterator begin() const { return VList.begin(); }
      const_iterator end() const { return VList.end(); }

      static char const* prefix();

      ColorVertexList& operator*=(Transform const& t)
      {
	 for (unsigned i = 0; i < VList.size(); ++i)
	 {
	    VList[i].first *= t;
	 }
	 return *this;
      }

   private:
      std::vector<value_type> VList;
};

std::ostream& operator<<(std::ostream& out, std::pair<Vertex, Color> const& VC);

std::ostream& operator<<(std::ostream& out, ColorVertexList const& VList);

class Face : public std::vector<int>
{
   public:
      Face() {}
      Face(size_t Size) : std::vector<int>(Size) {}
      template <class Iter>
	 Face(Iter first, Iter last) : std::vector<int>(first, last) {}
};

std::ostream& operator<<(std::ostream& out, Face const& f);

} // namespace oogl

#endif
