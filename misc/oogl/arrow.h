// -*- C++ -*- $Id$

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
