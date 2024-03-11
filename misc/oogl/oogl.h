// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// misc/oogl/oogl.h
//
// Copyright (C) 2002-2016 Ian McCulloch <ian@qusim.net>
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

/*
  classes to write OOGL files for use with geomview.

  Created 2002-02-13 Ian McCulloch

  In OOGL, the usual CG convention is used;
  vertices are row-vectors, and transformations are applied by applying the
  transform to the right, ie v' = v * A.
  This allows the operation v *= A.

  This also means that transformations apply in the opposite order, ie
  v * (A * B) is equivalent to (v * A) * B, where transformation A is applied first.
  In maths/physics notation, this would be (B * A) * v = B * (A * v).

  Does OOGL also use a left-handed coordinate system?
*/

#if !defined(OOGL_H_EYUIREWUICDH4W75467FY4376UIHRQ4)
#define OOGL_H_EYUIREWUICDH4W75467FY4376UIHRQ4

#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <utility>
#include "vertex.h"
#include "appearance.h"

namespace oogl
{

//
// OoglPrimitive is a base class for all oogl-symbols.
// The 'primitive' is only one component of an oogl-object.
// An oogl-object also includes a possible symbol definition and
// and appearance.
//
// OoglPrimitive objects contained within OoglObject have deep-copy semantics.
// (perhaps we need a null primitive?)
//

class OoglPrimitive
{
   public:
      OoglPrimitive() {}

      virtual OoglPrimitive* clone() const = 0;

      virtual void write(std::ostream& out) const = 0;

      virtual ~OoglPrimitive() {}
};

inline
std::ostream& operator<<(std::ostream& out, OoglPrimitive const& p)
{
   p.write(out);
   return out;
}

class OoglObject
{
   public:
      OoglObject();
      OoglObject(std::string const& SymbolName_, Appearance const& Appear_, OoglPrimitive const& Data_);
      OoglObject(Appearance const& Appear_, OoglPrimitive const& Data_);
      OoglObject(std::string const& SymbolName_, OoglPrimitive const& Data_);
      OoglObject(OoglPrimitive const& Data_); // not explicit, can be used as conversion ctor

      OoglObject(OoglObject const& Obj);

      OoglObject& operator=(OoglObject const& Obj);

      ~OoglObject();

      std::string const& symbol_name() const { return SymbolName; }

      Appearance const& appearance() const { return Appear; }

      Appearance& appearance() { return Appear; }

      OoglPrimitive const& primitive() const { return *Data; }
      void set_primitive(OoglPrimitive const& Data_);

   private:
      std::string SymbolName;
      Appearance Appear;
      OoglPrimitive* Data;
};

std::ostream& operator<<(std::ostream& out, OoglObject const& Obj);

//
// derived classes of OoglPrimitive
//

class TransformList : public OoglPrimitive
{
   private:
      typedef std::list<Transform> ContainerType;

   public:
      typedef Transform value_type;
      typedef ContainerType::const_iterator const_iterator;
      typedef ContainerType::iterator       iterator;

      TransformList() {}

      virtual OoglPrimitive* clone() const;

      virtual void write(std::ostream& out) const;

      iterator begin() { return Data.begin(); }
      iterator end() { return Data.end(); }

      const_iterator begin() const { return Data.begin(); }
      const_iterator end() const { return Data.end(); }

      int size() const { return Data.size(); }

      void push_back(Transform const& t) { Data.push_back(t); }

   private:
      ContainerType Data;
};

// Polyhedra represents the OFF oogl-object type.

// PolyhedraBase implements common functionality between variants.

class PolyhedraBase : public OoglPrimitive
{
   public:
      typedef Face                            value_type;
      typedef std::list<Face>::iterator       iterator;
      typedef std::list<Face>::const_iterator const_iterator;

      PolyhedraBase();

      void append_face(int c1, int c2, int c3);
      void append_face(int c1, int c2, int c3, int c4);
      void append_face(int c1, int c2, int c3, int c4, int c5);
      void append_face(int c1, int c2, int c3, int c4, int c5, int c6);
      void append_face(int c1, int c2, int c3, int c4, int c5, int c6, int c7);
      void append_face(int c1, int c2, int c3, int c4, int c5, int c6, int c7, int c8);

      template <class Iter>
      void append_face(Iter first, Iter last);

      void append_face(Face const& f) { Faces.push_back(f); }

      int size() const { return Faces.size(); }

      int num_faces() const { return size(); }
      int num_edges() const;

      iterator begin() { return Faces.begin(); }
      iterator end() { return Faces.end(); }

      const_iterator begin() const { return Faces.begin(); }
      const_iterator end() const { return Faces.end(); }

   private:
      std::list<Face> Faces;
};

template <class Iter>
void
PolyhedraBase::append_face(Iter first, Iter last)
{
   Faces.push_back(Face(first, last));
}


// BasicPolyhedra implements a polyhedra for an arbitary type of vertex list

template <class VertexListType>
class BasicPolyhedra : public PolyhedraBase
{
   public:
      BasicPolyhedra();

      BasicPolyhedra(VertexListType const& VList_);

      virtual OoglPrimitive* clone() const;

      virtual void write(std::ostream& out) const;

      void set_vertex_list(VertexListType const& VList_);

      VertexListType const& vertex_list() const { return VList; }

      int num_vertices() const { return vertex_list().size(); }

      BasicPolyhedra& operator*=(Transform const& t);

   private:
      VertexListType VList;
};

// typedefs for the common Polyhedra types

// Polyhedra is the basic OFF type
typedef BasicPolyhedra<VertexList> Polyhedra;

// Polyhedra is the COFF type, for vertex colors
typedef BasicPolyhedra<ColorVertexList> ColorPolyhedra;

// VectorList encapsulates the VECT format
class VectorList : public OoglPrimitive
{
   public:
      VectorList();

      virtual OoglPrimitive* clone() const;

      virtual void write(std::ostream& out) const;

      void append_polyline(VertexList const& VList);
      void append_polyline_closed(VertexList const& VList);

      void append_polyline(VertexList const& VList, Color const& c);
      void append_polyline_closed(VertexList const& VList, Color const& c);

      void append_polyline(ColorVertexList const& VList);
      void append_polyline_closed(ColorVertexList const& VList);

   private:
      std::list<Vertex> Vertices;
      std::vector<int> NVertices;      // needs to be supplied; -ve indicates closed polyline
      std::vector<int> NColors;
      std::list<Color> Colors;
};

// this encapsulates the INST format, for applying a transformation
// list to an object.
class ObjectTransform : public OoglPrimitive
{
   public:
      // The location types are:
      // global: The global coordinate system
      // camera: Relative to the camera position
      // ndc: Normalized unit cube of the camera projection
      // screen: screen coordinates, X,Y are in pixels, (0,0) is lower left of the screen.
      //         Z coordinate is normalized to -1 -> near clipping plane, +1 -> far clipping plane
      // local (default): relative to the parent object
      enum location_type {global, camera, ndc, screen, local};

      static char location_string[5][7];

      ObjectTransform();

      ObjectTransform(location_type Locat, location_type Orig, Vertex const& OrigP,
                      OoglObject const& G, TransformList const& TL);

      ObjectTransform(location_type Locat, OoglObject const& G, TransformList const& TL);

      ObjectTransform(location_type Orig, Vertex const& OrigP,
                      OoglObject const& G, TransformList const& TL);

      ObjectTransform(OoglObject const& G, TransformList const& TL);

      virtual OoglPrimitive* clone() const;

      virtual void write(std::ostream& out) const;

      location_type location() const { return Location; }
      location_type origin() const { return Origin; }

      Vertex origin_point() const { return OriginPoint; }

      OoglObject const& geom() const { return Geom; }

      TransformList const& transformation_list() const { return TList; }

   private:
      location_type Location, Origin;
      Vertex OriginPoint;
      OoglObject Geom;
      TransformList TList;
};

class List : public OoglPrimitive
{
   public:
      typedef OoglObject                              value_type;
      typedef value_type&                             reference;
      typedef value_type*                             pointer;
      typedef std::list<value_type>::iterator         iterator;
      typedef std::list<value_type>::const_iterator   const_iterator;

      List() {}
      ~List() {}
      List(List const& L);
      List& operator=(List const& L);

      virtual OoglPrimitive* clone() const;

      virtual void write(std::ostream& out) const;

      iterator begin() { return Data.begin(); }
      iterator end() { return Data.end(); }

      const_iterator begin() const { return Data.begin(); }
      const_iterator end() const { return Data.end(); }

      int size() const { return Data.size(); }

      void push_back(OoglObject const& Obj) { Data.push_back(Obj); }

   private:
      std::list<OoglObject> Data;
};

class External : public OoglPrimitive
{
   public:
      External(std::string const& FileName_);

      virtual OoglPrimitive* clone() const;

      virtual void write(std::ostream& out) const;

   private:
      std::string FileName;
};

class DefinedSymbol : public OoglPrimitive
{
   public:
      DefinedSymbol(std::string const& SymbolName_);

      virtual OoglPrimitive* clone() const;

      virtual void write(std::ostream& out) const;

   private:
      std::string SymbolName;
};

} // namespace



#endif
