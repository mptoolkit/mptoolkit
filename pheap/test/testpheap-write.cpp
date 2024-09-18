// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// pheap/test/testpheap-write.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#include "pheap/pheap.h"
#include "pheap/pvalueptr.h"

typedef std::pair<long, long> mypair;
typedef pvalue_ptr<mypair> ppair;

int main()
{
   pheap::PHeapLog.SetStream(std::cerr);
   pheap::PHeapLog.SetThreshold(100);

   pheap::Initialize("testdata", 1, 16384, 16*1024*1024);

   ppair N = new mypair(100, 200);

   pheap::ShutdownPersistent(N);
}
