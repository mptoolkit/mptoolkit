// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// misc/oogl/oogl.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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

#include "oogl.h"
#include <math.h>
#include "common/sortsearch.h"
#include <algorithm>
#include <set>
#include <sstream>
#include <ostream>
#include <iterator>

namespace oogl
{

Vertex& Vertex::operator*=(double a)
{
   x *= a;
   y *= a;
   z *= a;
   return *this;
}

Vertex& Vertex::operator+=(Vertex const& v)
{
   x += v.x;
   y += v.y;
   z += v.z;
   return *this;
}

Vertex& Vertex::operator-=(Vertex const& v)
{
   x -= v.x;
   y -= v.y;
   z -= v.z;
   return *this;
}

Vertex& Vertex::operator*=(Transform const& t)
{
   double xNew = x * t(0,0) + y * t(1,0) + z * t(2,0) + t(3,0);
   double yNew = x * t(0,1) + y * t(1,1) + z * t(2,1) + t(3,1);
   z = x * t(0,2) + y * t(1,2) + z * t(2,2) + t(3,2);
   x = xNew;
   y = yNew;
   return *this;
}

Vertex cylindrical(double r, double theta, double z)
{
   return Vertex(r * cos(theta), r * sin(theta), z);
}

Vertex spherical(double r, double theta, double phi)
{
   return Vertex(r * cos(theta) * sin(phi), r * sin(theta) * sin(phi), r * cos(phi));
}

Vertex cartesian(double x, double y, double z)
{
   return Vertex(x,y,z);
}

Vertex operator-(Vertex const& v)
{
   return Vertex(-v.x, -v.y, -v.z);
}

Vertex operator+(Vertex const& v1, Vertex const& v2)
{
   return Vertex(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

Vertex operator-(Vertex const& v1, Vertex const& v2)
{
   return Vertex(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

Vertex operator*(double a, Vertex const& v)
{
   return Vertex(a * v.x, a * v.y, a * v.z);
}

Vertex operator*(Vertex const& v, double a)
{
   return Vertex(v.x * a, v.y * a, v.z * a);
}

bool operator==(Vertex const& v1, Vertex const& v2)
{
   return v1.x == v2.x && v1.y == v2.y && v1.z == v2.z;
}

bool operator!=(Vertex const& v1, Vertex const& v2)
{
   return v1.x != v2.x || v1.y != v2.y || v1.z != v2.z;
}

std::ostream& operator<<(std::ostream& out, Vertex const& v)
{
   return out << v.x << ' ' << v.y << ' ' << v.z;
}

//
// Color
//

Color Color::HSV(double H, double S, double V, double A)
{
   H /= 60;                        // sector 0 to 5
   int i = int(floor(H));
   double f = H - i;               // fractional part of H
   double p = V * (1 - S);
   double q = V * (1 - S * f);
   double t = V * (1 - S * (1 - f));
   switch (i)
   {
      case 0  : return Color(V,t,p,A);
      case 1  : return Color(q,V,p,A);
      case 2  : return Color(p,V,t,A);
      case 3  : return Color(p,q,V,A);
      case 4  : return Color(t,p,V,A);
      default : return Color(V,p,q,A);
   }
}

std::ostream& operator<<(std::ostream& out, Color const& c)
{
   // force printing the decimal place, to avoid confusion with the
   // integer (0..255) format
   std::ostringstream Str;
   Str.setf(std::ios_base::fixed, std::ios_base::floatfield);
   Str.precision(6);
   Str << c.R << ' ' << c.G << ' ' << c.B << ' ' << c.A;
   return out << Str.str();
}

//
// Transform
//

Transform Transform::identity = Transform(ident_tag());

Transform::Transform(ident_tag)
{
   for (int i = 0; i < 4; ++i)
   {
      for (int j = 0; j < 4; ++j)
      {
         TransformMat[i][j] = 0;
      }
      TransformMat[i][i] = 1.0;
   }
}

Transform& Transform::operator*=(Transform const& t)
{
   *this = *this * t;
   return *this;
}

Transform Transform::rot_x(double theta)
{
   Transform Ret = Transform::identity;
   Ret(1,1) = Ret(2,2) = cos(theta);
   Ret(2,1) = -sin(theta);
   Ret(1,2) = sin(theta);
   return Ret;
}

Transform Transform::rot_y(double theta)
{
   Transform Ret = Transform::identity;
   Ret(0,0) = Ret(2,2) = cos(theta);
   Ret(0,2) = -sin(theta);
   Ret(2,0) = sin(theta);
   return Ret;
}

Transform Transform::rot_z(double theta)
{
   Transform Ret = Transform::identity;
   Ret(0,0) = Ret(1,1) = cos(theta);
   Ret(1,0) = -sin(theta);
   Ret(0,1) = sin(theta);
   return Ret;
}

Transform Transform::trans_x(double x)
{
   Transform Ret = Transform::identity;
   Ret(3,0) = x;
   return Ret;
}

Transform Transform::trans_y(double y)
{
   Transform Ret = Transform::identity;
   Ret(3,1) = y;
   return Ret;
}

Transform Transform::trans_z(double z)
{
   Transform Ret = Transform::identity;
   Ret(3,2) = z;
   return Ret;
}

Transform Transform::scale(double s)
{
   Transform Ret = Transform::identity;
   Ret(0,0) = s;
   Ret(1,1) = s;
   Ret(2,2) = s;
   return Ret;
}

Transform Transform::scale_x(double s)
{
   Transform Ret = Transform::identity;
   Ret(0,0) = s;
   return Ret;
}

Transform Transform::scale_y(double s)
{
   Transform Ret = Transform::identity;
   Ret(1,1) = s;
   return Ret;
}

Transform Transform::scale_z(double s)
{
   Transform Ret = Transform::identity;
   Ret(2,2) = s;
   return Ret;
}

Transform Transform::translate(Vertex const& v)
{
   Transform Ret = Transform::identity;
   Ret(3,0) = v.x;
   Ret(3,1) = v.y;
   Ret(3,2) = v.z;
   return Ret;
}

Transform operator*(Transform const& t1, Transform const& t2)
{
   Transform Ret;
   for (int i = 0; i < 4; ++i)
   {
      for (int j = 0; j < 4; ++j)
      {
         Ret(i,j) = t1(i,0) * t2(0,j) + t1(i,1) * t2(1,j) + t1(i,2) * t2(2,j) + t1(i,3) * t2(3,j);
      }
   }
   return Ret;
}

Vertex operator*(Vertex const& v, Transform const& t)
{
   return Vertex(v.x * t(0,0) + v.y * t(1,0) + v.z * t(2,0) + t(3,0),
                 v.x * t(0,1) + v.y * t(1,1) + v.z * t(2,1) + t(3,1),
                 v.x * t(0,2) + v.y * t(1,2) + v.z * t(2,2) + t(3,2));
}

std::ostream& operator<<(std::ostream& out, Transform const& t)
{
   for (int i = 0; i < 4; ++i)
   {
      for (int j = 0; j < 4; ++j)
      {
         out << t(i,j) << ' ';
      }
      out << '\n';
   }
   return out;
}

//
// TransformList
//

OoglPrimitive* TransformList::clone() const
{
   return new TransformList(*this);
}

void TransformList::write(std::ostream& out) const
{
   out << "TLIST\n";
   std::copy(begin(), end(), std::ostream_iterator<Transform>(out, "\n"));
}

//
// VertexList
//

std::ostream& operator<<(std::ostream& out, VertexList const& VList)
{
   std::copy(VList.begin(), VList.end(), std::ostream_iterator<VertexList::value_type>(out, "\n"));
   return out;
}

char const* VertexList::prefix()
{
   return "";
}

//
// ColorVertexList
//

std::ostream& operator<<(std::ostream& out, std::pair<Vertex, Color> const& VC)
{
   return out << VC.first << ' ' << VC.second;
}

std::ostream& operator<<(std::ostream& out, ColorVertexList const& VList)
{
   std::copy(VList.begin(), VList.end(), std::ostream_iterator<ColorVertexList::value_type>(out, "\n"));
   return out;
}

char const* ColorVertexList::prefix()
{
   return "C";
}

//
// Face
//

std::ostream& operator<<(std::ostream& out, Face const& f)
{
   out << f.size() << ' ';
   std::copy(f.begin(), f.end(), std::ostream_iterator<int>(out, " "));
   return out;
}

//
// PolyhedraBase
//

PolyhedraBase::PolyhedraBase()
{
}


void
PolyhedraBase::append_face(int c1, int c2, int c3)
{
   Face f(3);
   f[0] = c1;
   f[1] = c2;
   f[2] = c3;
   Faces.push_back(f);
}

void
PolyhedraBase::append_face(int c1, int c2, int c3, int c4)
{
   Face f(4);
   f[0] = c1;
   f[1] = c2;
   f[2] = c3;
   f[3] = c4;
   Faces.push_back(f);
}

void
PolyhedraBase::append_face(int c1, int c2, int c3, int c4, int c5)
{
   Face f(5);
   f[0] = c1;
   f[1] = c2;
   f[2] = c3;
   f[3] = c4;
   f[4] = c5;
   Faces.push_back(f);
}

void
PolyhedraBase::append_face(int c1, int c2, int c3, int c4, int c5, int c6)
{
   Face f(6);
   f[0] = c1;
   f[1] = c2;
   f[2] = c3;
   f[3] = c4;
   f[4] = c5;
   f[5] = c6;
   Faces.push_back(f);
}

void
PolyhedraBase::append_face(int c1, int c2, int c3, int c4, int c5, int c6, int c7)
{
   Face f(7);
   f[0] = c1;
   f[1] = c2;
   f[2] = c3;
   f[3] = c4;
   f[4] = c5;
   f[5] = c6;
   f[6] = c7;
   Faces.push_back(f);
}

void
PolyhedraBase::append_face(int c1, int c2, int c3, int c4, int c5, int c6, int c7, int c8)
{
   Face f(8);
   f[0] = c1;
   f[1] = c2;
   f[2] = c3;
   f[3] = c4;
   f[4] = c5;
   f[5] = c6;
   f[6] = c7;
   f[7] = c8;
   Faces.push_back(f);
}

int
PolyhedraBase::num_edges() const
{
   std::set<std::pair<int, int> > EdgeList;
   for (const_iterator I = begin(); I != end(); ++I)
   {
      for (unsigned i = 0; i < I->size() - 1; ++i)
      {
         int v1 = (*I)[i];
         int v2 = (*I)[i+1];
         sort2(v1,v2);
         EdgeList.insert(std::make_pair(v1,v2));
      }
      if (I->size() > 2)
      {
         int v1 = (*I)[0];
         int v2 = (*I)[I->size()-1];
         sort2(v1,v2);
         EdgeList.insert(std::make_pair(v1,v2));
      }
   }
   return EdgeList.size();
}

//
// VectorList
//

VectorList::VectorList()
{
}

OoglPrimitive* VectorList::clone() const
{
  return new VectorList(*this);
}

void VectorList::write(std::ostream& out) const
{
   out << "VECT\n";
   out << NVertices.size() << ' ' << Vertices.size() << ' ' << Colors.size() << '\n';
   std::copy(NVertices.begin(), NVertices.end(), std::ostream_iterator<int>(out, " "));
   out << '\n';
   std::copy(NColors.begin(), NColors.end(), std::ostream_iterator<int>(out, " "));
   out << '\n';
   std::copy(Vertices.begin(), Vertices.end(), std::ostream_iterator<Vertex>(out, "\n"));
   std::copy(Colors.begin(), Colors.end(), std::ostream_iterator<Color>(out, " "));
}

void VectorList::append_polyline(VertexList const& VList)
{
   NVertices.push_back(VList.size());
   NColors.push_back(0);

   std::copy(VList.begin(), VList.end(), std::back_inserter(Vertices));
}

void VectorList::append_polyline_closed(VertexList const& VList)
{
   NVertices.push_back(-VList.size());
   NColors.push_back(0);

   std::copy(VList.begin(), VList.end(), std::back_inserter(Vertices));
}

void VectorList::append_polyline(VertexList const& VList, Color const& c)
{
   NVertices.push_back(VList.size());
   std::copy(VList.begin(), VList.end(), std::back_inserter(Vertices));

   NColors.push_back(1);
   Colors.push_back(c);
}

void VectorList::append_polyline_closed(VertexList const& VList, Color const& c)
{
   NVertices.push_back(-VList.size());
   NColors.push_back(1);

   std::copy(VList.begin(), VList.end(), std::back_inserter(Vertices));
   Colors.push_back(c);
}

void VectorList::append_polyline(ColorVertexList const& VList)
{
   NVertices.push_back(VList.size());
   NColors.push_back(VList.size());

   for (ColorVertexList::const_iterator I = VList.begin(); I != VList.end(); ++I)
   {
      Vertices.push_back(I->first);
      Colors.push_back(I->second);
   }
}

void VectorList::append_polyline_closed(ColorVertexList const& VList)
{
   NVertices.push_back(-VList.size());
   NColors.push_back(VList.size());

   for (ColorVertexList::const_iterator I = VList.begin(); I != VList.end(); ++I)
   {
      Vertices.push_back(I->first);
      Colors.push_back(I->second);
   }
}

//
// OoglObject
//

OoglObject::OoglObject()
  : SymbolName(), Appear(), Data(NULL)
{
}

OoglObject::OoglObject(std::string const& SymbolName_, Appearance const& Appear_, OoglPrimitive const& Data_)
  : SymbolName(SymbolName_), Appear(Appear_), Data(Data_.clone())
{
}

OoglObject::OoglObject(Appearance const& Appear_, OoglPrimitive const& Data_)
  : SymbolName(), Appear(Appear_), Data(Data_.clone())
{
}

OoglObject::OoglObject(std::string const& SymbolName_, OoglPrimitive const& Data_)
  : SymbolName(SymbolName_), Appear(), Data(Data_.clone())
{
}

OoglObject::OoglObject(OoglPrimitive const& Data_)
  : SymbolName(), Appear(), Data(Data_.clone())
{
}

OoglObject::OoglObject(OoglObject const& Obj)
  : SymbolName(Obj.SymbolName), Appear(Obj.Appear), Data(Obj.Data->clone())
{
}

OoglObject& OoglObject::operator=(OoglObject const& Obj)
{
   SymbolName = Obj.SymbolName;
   Appear = Obj.Appear;
   OoglPrimitive* Temp = Obj.Data->clone();
   delete Data;
   Data = Temp;
   return *this;
}

OoglObject::~OoglObject()
{
   delete Data;
}

void OoglObject::set_primitive(OoglPrimitive const& Data_)
{
   OoglPrimitive* Temp = Data_.clone();
   delete Data;
   Data = Temp;
}

std::ostream& operator<<(std::ostream& out, OoglObject const& Obj)
{
   out << "{ ";
   if (!Obj.symbol_name().empty())
     out << "define " << Obj.symbol_name() << ' ';

   if (!Obj.appearance().empty())
     out << "appearance {\n" << Obj.appearance() << "}\n";

   out << Obj.primitive() << " }";
   return out;
}

//
// List
//

List::List(List const& L)
{
   for (const_iterator I = L.begin(); I != L.end(); ++I)
   {
      Data.push_back((*I));
   }
}

List& List::operator=(List const& L)
{
   // this is not exception safe
   if (&L == this) return *this;

   Data.clear();
   for (const_iterator I = L.begin(); I != L.end(); ++I)
   {
      Data.push_back((*I));
   }

   return *this;
}

OoglPrimitive* List::clone() const
{
   return new List(*this);
}

void List::write(std::ostream& out) const
{
   out << "LIST\n";

   for (const_iterator I = begin(); I != end(); ++I)
   {
      out << *I << '\n';
   }
}

//
// External
//

External::External(std::string const& FileName_)
  : FileName(FileName_)
{
}

OoglPrimitive* External::clone() const
{
   return new External(FileName);
}

void External::write(std::ostream& out) const
{
   out << " < " << FileName;
}

//
// DefinedSymbol
//

DefinedSymbol::DefinedSymbol(std::string const& SymbolName_)
  : SymbolName(SymbolName_)
{
}

OoglPrimitive* DefinedSymbol::clone() const
{
   return new DefinedSymbol(SymbolName);
}

void DefinedSymbol::write(std::ostream& out) const
{
   out << " : " << SymbolName;
}

//
// OjbectTransform
//

char ObjectTransform::location_string[5][7] = {"global", "camera", "ndc", "screen", "local"};


ObjectTransform::ObjectTransform()
  : Location(local), Origin(local), OriginPoint(0,0,0),
  Geom(), TList()
{
}

ObjectTransform::ObjectTransform(location_type Locat, location_type Orig, Vertex const& OrigP,
                OoglObject const& G, TransformList const& TL)
  : Location(Locat), Origin(Orig), OriginPoint(OrigP),
  Geom(G), TList(TL)
{
}

ObjectTransform::ObjectTransform(location_type Locat, OoglObject const& G, TransformList const& TL)
  : Location(Locat), Origin(local), OriginPoint(0,0,0),
  Geom(G), TList(TL)
{
}

ObjectTransform::ObjectTransform(location_type Orig, Vertex const& OrigP,
                OoglObject const& G, TransformList const& TL)
  : Location(local), Origin(Orig), OriginPoint(OrigP),
  Geom(G), TList(TL)
{
}

ObjectTransform::ObjectTransform(OoglObject const& G, TransformList const& TL)
  : Location(local), Origin(local), OriginPoint(0,0,0),
  Geom(G), TList(TL)
{
}

OoglPrimitive* ObjectTransform::clone() const
{
   return new ObjectTransform(Location, Origin, OriginPoint, Geom, TList);
}

void ObjectTransform::write(std::ostream& out) const
{
   out << "INST\n";

   if (Location != local)
     out << "location " << location_string[Location] << '\n';

   if (Origin != local && OriginPoint != Vertex(0,0,0))
     out << "origin " << location_string[Origin] << ' ' << OriginPoint << '\n';

   out << "geom " << Geom << '\n';

   if (TList.size() > 0)
     out << "transforms { " << TList << "}\n";
}

//
// pre-defined colors
//

namespace Colors
{

   Color const black = Color(0,0,0);
   Color const white = Color(1,1,1);
   Color const red = Color(1,0,0);
   Color const green = Color(0,1,0);
   Color const blue = Color(0,0,1);
   Color const cyan = Color(0,1,1);
   Color const yellow = Color(1,1,0);
   Color const magenta = Color(1,0,1);

} // namespace Colors

//
// Appearance
//

bool operator==(Appearance const& a1, Appearance const& a2)
{
   return a1.ShowFace == a2.ShowFace
      && a1.ShowEdge == a2.ShowEdge
      && a1.ShowVECT == a2.ShowVECT
      ;
}

bool Appearance::MaterialProperties::empty() const
{
   return AmbientReflection.is_empty() && AmbientColor.is_empty()
      && DiffuseReflection.is_empty() && DiffuseColor.is_empty()
      && SpecularReflection.is_empty() && SpecularColor.is_empty()
      && Shininess.is_empty() && BackDiffuseColor.is_empty()
      && Alpha.is_empty() && EdgeColor.is_empty() && NormalColor.is_empty();
}

bool Appearance::empty() const
{
   return ShowFace.is_empty() && ShowEdge.is_empty() && ShowVECT.is_empty()
      && ShowTransparent.is_empty() && ShowNormal.is_empty()
      && NormalLength.is_empty() && Evert.is_empty() && ShowTextures.is_empty()
      && BackCull.is_empty() && HandleConcave.is_empty() && ShadeLines.is_empty()
      && Shading.is_empty() && LineWidth.is_empty() && BezierDivide.is_empty()
      && Material.empty();
}

void ShowProperty(std::ostream& out, std::string Name, Property<bool> const& p)
{
   if (p.is_active())
   {
      out << (p.value() ? '+' : '-') << Name << '\n';
   }
}

void ShowProperty(std::ostream& out, std::string Name, Property<int> const& p)
{
   if (p.is_active())
   {
      out << Name << ' ' << p.value() << '\n';
   }
}

void ShowProperty(std::ostream& out, std::string Name, Property<double> const& p)
{
   if (p.is_active())
   {
      out << Name << ' ' << p.value() << '\n';
   }
}

void ShowProperty(std::ostream& out, std::string Name, Property<std::pair<int, int> > const& p)
{
   if (p.is_active())
   {
      out << Name << ' ' << p.value().first << ' ' << p.value().second << '\n';
   }
}

void ShowProperty(std::ostream& out, std::string Name, Property<Color> const& p)
{
   if (p.is_active())
   {
      out << std::showpoint;
      out << Name << ' ' << p.value().R << ' ' << p.value().G << ' ' << p.value().B << '\n';
   }
}

void ShowProperty(std::ostream& out, std::string Name, Property<Appearance::ShadingTypes> const& p)
{
   if (p.is_active())
   {
      out << Name << ' ';
      switch (p.value())
      {
         case Appearance::Constant : out << "constant"; break;
         case Appearance::Flat : out << "flat"; break;
         case Appearance::Smooth : out << "smooth"; break;
         case Appearance::CSmooth : out << "csmooth"; break;
      }
      out << '\n';
   }
}

std::ostream& operator<<(std::ostream& out, Appearance const& a)
{
   ShowProperty(out, "face", a.ShowFace);
   ShowProperty(out, "edge", a.ShowEdge);
   ShowProperty(out, "vect", a.ShowVECT);
   ShowProperty(out, "transparent", a.ShowTransparent);
   ShowProperty(out, "normal", a.ShowNormal);
   ShowProperty(out, "normscale", a.NormalLength);
   ShowProperty(out, "evert", a.Evert);
   ShowProperty(out, "texturing", a.ShowTextures);
   ShowProperty(out, "backcull", a.BackCull);
   ShowProperty(out, "concave", a.HandleConcave);
   ShowProperty(out, "shadelines", a.ShadeLines);
   ShowProperty(out, "shading", a.Shading);
   ShowProperty(out, "linewidth", a.LineWidth);
   ShowProperty(out, "patchdice", a.BezierDivide);

   // material properties
   if (!a.Material.empty())
   {
      out << "material {\n";
      ShowProperty(out, "ka", a.Material.AmbientReflection);
      ShowProperty(out, "ambient", a.Material.AmbientColor);
      ShowProperty(out, "kd", a.Material.DiffuseReflection);
      ShowProperty(out, "diffuse", a.Material.DiffuseColor);
      ShowProperty(out, "ks", a.Material.SpecularReflection);
      ShowProperty(out, "specular", a.Material.SpecularColor);
      ShowProperty(out, "shininess", a.Material.Shininess);
      ShowProperty(out, "backdiffuse", a.Material.BackDiffuseColor);
      ShowProperty(out, "alpha", a.Material.Alpha);
      ShowProperty(out, "edgecolor", a.Material.EdgeColor);
      ShowProperty(out, "normalcolor", a.Material.NormalColor);
      out << "}\n";
   }

   return out;
}

} // namespace oogl
