// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// misc/oogl/arrow.h
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

#if !defined(ARROW_H_SDFHWERUIOYF34978Y89UJF89HS89P)
#define ARROW_H_SDFHWERUIOYF34978Y89UJF89HS89P

#include "oogl.h"

oogl::ColorPolyhedra CreateColorArrow(double BackLength, 
                                      double FrontLength, 
                                      double InnerRadius, 
                                      double OuterRadius, 
                                      int NumFaces,
                                      oogl::Color BackColor, 
                                      oogl::Color MidColor, 
                                      oogl::Color PointColor);

oogl::ColorPolyhedra Create2DArrow(double BackLength, 
                                   double FrontLength, 
                                   double MinWidth, 
                                   double MaxWidth,
				   oogl::Color BackColor, 
                                   oogl::Color MidColor,
                                   oogl::Color PointColor);

oogl::ColorPolyhedra CreateBox(double XLength, 
                               double YLength, 
                               oogl::Color Color);

#endif
