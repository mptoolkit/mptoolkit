// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// pheap/test/testpheap.cpp
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

#include "pheap/pheap.h"
#include "pheap/pvalueptr.h"

typedef pvalue_ptr<std::pair<int, int> > ppair;

int main()
{
   pheap::PHeapLog.SetStream(std::cerr);
   pheap::PHeapLog.SetThreshold(100);

   pheap::Initialize("testdata", 1, 1024*1024, 16*1024*1024);

   ppair N = new std::pair<int, int>(100, 200);

   pheap::ShutdownPersistent(N);
}
