// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// misc/oogl/appearance.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
// class Appearance and related properties for OOGL
//

#if !defined(APPEARANCE_H_SJCHRUY438578Y)
#define APPEARANCE_H_SJCHRUY438578Y

#include <utility>

namespace oogl
{

template <typename T>
class Property
{
   public:
      Property() : Value(), Active(false) {}
      Property(T const& x) : Value(x), Active(true) {}

      bool is_active() const { return Active; }
      bool is_empty() const { return !Active; }

      T& operator=(T const& x) { Active = true; Value = x; return Value; }

      T value() const { return Value; }

   private:
      T Value;
      bool Active;

   friend void disable(Property<T>& p) { p.Active = false; }
};

template <typename T>
bool operator==(Property<T> const& p1, Property<T> const& p2)
{
   return ((p1.is_active() == p2.is_active()) && (!p1.is_active() || (p1.value() == p2.value())));
}

// An Appearance is a component of all oogl-objects
class Appearance
{
   public:
      bool empty() const;

      // Shading types: Constant = fixed colour, no lighting effects
      //                Flat = Apply lighing for each face
      //                Smooth = gouraud shading with lighing
      //                CSmooth = gouraud shading without light
      enum ShadingTypes {Constant, Flat, Smooth, CSmooth};

      Property<bool>   ShowFace;        // true = draw faces
      Property<bool>   ShowEdge;        // true = draw edges
      Property<bool>   ShowVECT;        // true = show VECTs
      Property<bool>   ShowTransparent; // true = show transparency
      Property<bool>   ShowNormal;      // true = show surface normals
      Property<double> NormalLength;    // length of surface normals in object coords
      Property<bool>   Evert;           // true = every normals to always face camera
      Property<bool>   ShowTextures;    // true = show textures
      Property<bool>   BackCull;        // true = cull backfaces
      Property<bool>   HandleConcave;   // true = handle concave shapes
      Property<bool>   ShadeLines;      // true = shade lines as if they were cylinders (GL only)
      Property<ShadingTypes> Shading;
      Property<int>    LineWidth;       // width in pixels of lines
      Property<std::pair<int, int> > BezierDivide; // divide patches this finely in u,v

      struct MaterialProperties
      {
         Property<double> AmbientReflection;  // coefficient of ambient reflection
         Property<Color>  AmbientColor;      // ambient color - alpha is ignored
         Property<double> DiffuseReflection; // coefficient of diffuse reflection
         Property<Color>  DiffuseColor;      // diffuse color - alpha is ignored
         Property<double> SpecularReflection; // coefficient of specular reflection
         Property<Color>  SpecularColor;      // specular color - alpha is ignored
         Property<double> Shininess;          // specular exponent - larger values = sharper highlights
         Property<Color>  BackDiffuseColor;   // back-face color for two-sided surfaces
         Property<double> Alpha;              // transparancy.  1 = opaque, 0 = transparent

         Property<Color>  EdgeColor;          // line & edge color
         Property<Color>  NormalColor;        // color for surface normal vectors

         bool empty() const;
      };

      MaterialProperties Material;
};

bool operator==(Appearance const& a, Appearance const& a2);

std::ostream& operator<<(std::ostream& out, Appearance const& a);

} // namespace oogl

#endif
