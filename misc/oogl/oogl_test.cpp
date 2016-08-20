// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// misc/oogl/oogl_test.cpp
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

#include "oogl.h"

using namespace oogl;
using namespace std;

int main()
{
   Vertex v1(1,0,0);
   Vertex v2(0,1,0);
   Vertex v3(0,0,1);

   Transform r1 = Transform::rot_x(3.1415927 / 4.0);
   Transform r2 = Transform::rot_y(3.1415927 / 4.0);
   Transform r3 = Transform::rot_z(3.1415927 / 4.0);

   Transform t1 = Transform::trans_x(1.5);
   Transform t2 = Transform::trans_y(2);
   Transform t3 = Transform::trans_z(2.5);

   cout.precision(12);

   std::cout << "v1 = " << v1 << '\n';
   std::cout << "v2 = " << v2 << '\n';
   std::cout << "v3 = " << v3 << '\n';

   std::cout << "r1 = \n" << r1 << '\n';
   std::cout << "r2 = \n" << r2 << '\n';
   std::cout << "r3 = \n" << r3 << '\n';

   std::cout << "t1 = \n" << t1 << '\n';
   std::cout << "t2 = \n" << t2 << '\n';
   std::cout << "t3 = \n" << t3 << '\n';

   std::cout << "v1*r1 = " << (v1*r1) << '\n';
   std::cout << "v1*r2 = " << (v1*r2) << '\n';
   std::cout << "v1*r3 = " << (v1*r3) << '\n';

   std::cout << "v2*r1 = " << (v2*r1) << '\n';
   std::cout << "v2*r2 = " << (v2*r2) << '\n';
   std::cout << "v2*r3 = " << (v2*r3) << '\n';

   std::cout << "v3*r1 = " << (v3*r1) << '\n';
   std::cout << "v3*r2 = " << (v3*r2) << '\n';
   std::cout << "v3*r3 = " << (v3*r3) << '\n';

   std::cout << "v1*t1 = " << (v1*t1) << '\n';
   std::cout << "v1*t2 = " << (v1*t2) << '\n';
   std::cout << "v1*t3 = " << (v1*t3) << '\n';

   std::cout << "v2*t1 = " << (v2*t1) << '\n';
   std::cout << "v2*t2 = " << (v2*t2) << '\n';
   std::cout << "v2*t3 = " << (v2*t3) << '\n';

   std::cout << "v3*t1 = " << (v3*t1) << '\n';
   std::cout << "v3*t2 = " << (v3*t2) << '\n';
   std::cout << "v3*t3 = " << (v3*t3) << '\n';

   std::cout << "r1*r2 = \n" << (r1*r2) << '\n';
   std::cout << "r1*r3 = \n" << (r1*r3) << '\n';

   std::cout << "t1*t2 = \n" << (t1*t2) << '\n';
   std::cout << "t1*t3 = \n" << (t1*t3) << '\n';

   std::cout << "r1*t1 = \n" << (r1*t1) << '\n';
   std::cout << "r1*t2 = \n" << (r1*t2) << '\n';

   VertexList VList;
   VList.push_back(Vertex(-1,-1,-1));
   VList.push_back(Vertex(-1,-1,1));
   VList.push_back(Vertex(-1,1,-1));
   VList.push_back(Vertex(-1,1,1));
   VList.push_back(Vertex(1,-1,-1));
   VList.push_back(Vertex(1,-1,1));
   VList.push_back(Vertex(1,1,-1));
   VList.push_back(Vertex(1,1,1));

   std::cout << "VList = \n" << VList << '\n';

   Polyhedra MyCube(VList);
   MyCube.append_face(0,1,3,2);
   MyCube.append_face(0,1,5,4);
   MyCube.append_face(0,2,6,4);
   MyCube.append_face(7,6,4,5);
   MyCube.append_face(7,6,2,3);
   MyCube.append_face(7,5,1,3);

   std::cout << "MyCube = \n" << MyCube << '\n';

   TransformList TList;
   TList.push_back(r1);
   TList.push_back(r2);
   TList.push_back(r3);
   //   TList.push_back(t1);
   //   TList.push_back(t2);
   //   TList.push_back(t3);

   std::cout << "TList = \n" << TList << '\n';

   ObjectTransform Trans(MyCube, TList);

   std::cout << "Trans = \n" << Trans << '\n';

   return 0;
}
