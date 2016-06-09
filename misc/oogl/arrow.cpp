// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// misc/oogl/arrow.cpp
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

#include "arrow.h"
#include <math.h>
#include "math_const.h"

using namespace oogl;

ColorPolyhedra CreateColorArrow(double BackLength, double FrontLength, 
                                double InnerRadius, double OuterRadius, int NumFaces,
                                Color BackColor, Color MidColor, Color PointColor)
{
   // make the vertex list
   ColorVertexList Vertices;
   
   double const BackPlane = -BackLength;
   double const BasePoint = 0;
   double const Point = FrontLength;
   double const r = InnerRadius;
   double const rOuter = OuterRadius;

   double x,y;

   // base of the body of the arrow
   for (int i = 0; i < NumFaces; ++i)
   {
      x = cos(i * math_const::pi * 2 / NumFaces);
      y = sin(i * math_const::pi * 2 / NumFaces);
      Vertices.push_back(Vertex(x*r, y * r, BackPlane), BackColor);
      
   }

   // top of the body of the arrow
   for (int i = 0; i < NumFaces; ++i)
   {
      x = cos(i * math_const::pi * 2 / NumFaces);
      y = sin(i * math_const::pi * 2 / NumFaces);
      Vertices.push_back(Vertex(x * r, y * r, BasePoint), MidColor);
   }
 
   // base of the point of the arrow
   for (int i = 0; i < NumFaces; ++i)
   {
      x = cos(i * math_const::pi * 2 / NumFaces);
      y = sin(i * math_const::pi * 2 / NumFaces);
      Vertices.push_back(Vertex(x * rOuter, y * rOuter, BasePoint), MidColor);
   }

   // top of the arrow
   Vertices.push_back(Vertex(0,0,Point), PointColor);

   // make the actual arrow
   ColorPolyhedra MyArrow(Vertices);

   // make the faces

   // Base of the arrow
   std::vector<int> BaseArrowFace;
   for (int i = 0; i < NumFaces; ++i)
   {
      BaseArrowFace.push_back(i);
   }
   MyArrow.append_face(BaseArrowFace.begin(), BaseArrowFace.end());

   // faces of the body of the arrow
   for (int i = 0; i < NumFaces; ++i)
   {
      MyArrow.append_face(i, (i+1)%NumFaces, (i+1)%NumFaces + NumFaces, i+NumFaces);
   }

   // Base of the point
   std::vector<int> BaseFace;
   for (int i = 0; i < NumFaces; ++i)
   {
      BaseFace.push_back(NumFaces*2 + i);
   }
   MyArrow.append_face(BaseFace.begin(), BaseFace.end());
  
   // faces of the point
   for (int i = 0; i < NumFaces; ++i)
   {
      MyArrow.append_face(NumFaces*2+i, NumFaces*2 + (i+1)%NumFaces, NumFaces*3);
   }

   return MyArrow;
}

ColorPolyhedra Create2DArrow(double BackLength, double FrontLength, 
                             double MinWidth, double MaxWidth,
                             Color BackColor, Color MidColor, Color PointColor)
{
   // make the vertex list
   ColorVertexList Vertices;
   
   // bottom of arrow
   Vertices.push_back(Vertex(-MinWidth/2, -BackLength, 0), BackColor);
   Vertices.push_back(Vertex(MinWidth/2, -BackLength, 0), BackColor);

   // Central portion
   Vertices.push_back(Vertex(MinWidth/2, 0, 0), MidColor);
   Vertices.push_back(Vertex(-MinWidth/2, 0, 0), MidColor);

   // Top portion
   Vertices.push_back(Vertex(-MaxWidth/2, 0, 0), MidColor);
   Vertices.push_back(Vertex(MaxWidth/2, 0, 0), MidColor);
   Vertices.push_back(Vertex(0, FrontLength, 0), PointColor);

   ColorPolyhedra MyArrow(Vertices);

   MyArrow.append_face(0, 1, 2, 3);
   MyArrow.append_face(4, 5, 6);

   return MyArrow;
}

oogl::ColorPolyhedra CreateBox(double XLength, double YLength, oogl::Color Color)
{
   ColorVertexList Vertices;

   Vertices.push_back(Vertex(-XLength/2, -YLength/2, 0), Color);
   Vertices.push_back(Vertex(-XLength/2,  YLength/2, 0), Color);
   Vertices.push_back(Vertex( XLength/2,  YLength/2, 0), Color);
   Vertices.push_back(Vertex( XLength/2, -YLength/2, 0), Color);

   ColorPolyhedra MyPoly(Vertices);
   MyPoly.append_face(0,1,2,3);

   return MyPoly;
}
