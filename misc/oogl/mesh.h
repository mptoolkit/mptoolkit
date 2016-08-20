// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// misc/oogl/mesh.h
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
// This encapsulates the ZMESH format
//

#if !defined(MESH_H_DSJHFEWRY84YPO9HO)
#define MESH_H_DSJHFEWRY84YPO9HO

#include "oogl.h"
#include "linearalgebra/matrix.h"

namespace oogl
{

class ZMesh : public OoglPrimitive
{
   public:
      ZMesh() : WrapU(false), WrapV(false) {}
      ZMesh(unsigned usz, unsigned vsz, bool wrapU = false, bool wrapV = false)
         : Data(usz, vsz, 0.0),  WrapU(wrapU), WrapV(wrapV) {}

      ZMesh(LinearAlgebra::Matrix<double> const& M, bool wrapU = false, bool wrapV = false)
      : Data(M), WrapU(wrapU), WrapV(wrapV) {}

      std::size_t size1() const { return Data.size1(); }
      std::size_t size2() const { return Data.size2(); }

      // true if the mesh is wrapped in the u or v direction
      bool is_wrap_u() const { return WrapU; }
      bool is_wrap_v() const { return WrapV; }

      void wrap_u(bool w = true) { WrapU = w; }
      void wrap_v(bool w = true) { WrapV = w; }

      double operator()(unsigned u, unsigned v) const { return Data(u,v); }
      double& operator()(unsigned u, unsigned v) { return Data(u,v); }

      LinearAlgebra::Matrix<double> const& data() const { return Data; }
      LinearAlgebra::Matrix<double>& data() { return Data; }

      virtual OoglPrimitive* clone() const
      {
         return new ZMesh(Data, WrapU, WrapV);
      }

      virtual void write(std::ostream& out) const;

   private:
      LinearAlgebra::Matrix<double> Data;
      bool WrapU, WrapV;
};

class ColorZMesh : public ZMesh
{
   public:
      ColorZMesh() {}
      ColorZMesh(unsigned usz, unsigned vsz, bool wrapU = false, bool wrapV = false)
         : ZMesh(usz, vsz, wrapU, wrapV) {}

      ColorZMesh(LinearAlgebra::Matrix<double> const& M,
                 LinearAlgebra::Matrix<Color> const& C,
                 bool wrapU = false, bool wrapV = false)
         : ZMesh(M, wrapU, wrapV), Colors(C) {}

      Color color(unsigned u, unsigned v) const { return Colors(u,v); }
      Color& color(unsigned u, unsigned v) { return Colors(u,v); }

      virtual OoglPrimitive* clone() const
      {
         return new ColorZMesh(this->data(), Colors, this->is_wrap_u(), this->is_wrap_v());
      }

      virtual void write(std::ostream& out) const;

   private:
      LinearAlgebra::Matrix<Color> Colors;
};


// inlines

inline
void ZMesh::write(std::ostream& out) const
{
   out << 'Z';
   if (this->is_wrap_u()) out << 'u';
   if (this->is_wrap_v()) out << 'v';
   out << "MESH\n";
   out << this->size1() << ' ' << this->size2() << '\n';
   for (unsigned v = 0; v < this->size2(); ++v)
   {
      for (unsigned u = 0; u < this->size1(); ++u)
      {
         out << Data(u,v) << ' ';
      }
      out << '\n';
   }
}

inline
void ColorZMesh::write(std::ostream& out) const
{
   out << "CZ";
   if (this->is_wrap_u()) out << 'u';
   if (this->is_wrap_v()) out << 'v';
   out << "MESH\n";
   out << this->size1() << ' ' << this->size2() << '\n';
   for (unsigned v = 0; v < this->size2(); ++v)
   {
      for (unsigned u = 0; u < this->size1(); ++u)
      {
         out << (*this)(u,v) << ' ' << this->color(u,v) << ' ';
      }
      out << '\n';
   }
}

} // namespace oogl

#endif
