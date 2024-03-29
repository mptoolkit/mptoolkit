// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// interface/inittemp.cpp
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

#include "inittemp.h"
#include "common/environment.h"
#include "common/proccontrol.h"
#include "pheap/pheap.h"

namespace mp_pheap
{

int PageSize()
{
   return getenv_or_default("MP_PAGESIZE", DEFAULT_PAGE_SIZE);
}

long CacheSize()
{
   return getenv_or_default("MP_CACHESIZE", long(DEFAULT_PAGE_CACHE_SIZE));
}

void InitializeTempPHeap(bool Verbose)
{
   std::string TempFile = getenv_or_default("MP_BINPATH", std::string());
   if (Verbose > 0)
   {
      if (TempFile.empty())
         std::cerr << "Using temporary file";
      else
         std::cerr << "Using temporary filespec " << TempFile;
      std::cerr << ", with page size " << PageSize()
                << ", cache size " << CacheSize() << '\n';
   }
   bool ShouldUnlink = getenv_or_default("MP_UNLINKTEMP", 1);
   int TempFileDesc = ProcControl::CreateUniqueFile(TempFile, ShouldUnlink);
   CHECK(TempFileDesc != -1);  // if this happens, then we failed to create a temporary file
   pheap::Initialize(TempFile, 1, PageSize(), CacheSize(), ShouldUnlink);
}

} // namespace mp_pheap
