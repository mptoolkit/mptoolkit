// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// pheap/pheapallocator.cpp
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

#include "pheapallocator.h"

namespace PHeapFileSystem
{

void Descriptor::Write(BlockFileSystem* FS_, PStream::opstream& out) const
{
   out.write<uint32>(BlockList.size());
   for (BlockListType::const_iterator I = BlockList.begin(); I != BlockList.end(); ++I)
   {
      FS_->WriteBlockRecord(out, *I);
   }
   out << IdList << uint32(DataSize);
}

Descriptor::Descriptor(BlockFileSystem* FS_, PStream::ipstream& in)
{
   size_t Size = in.read<uint32>();
   for (size_t i = 0; i < Size; ++i)
   {
      BlockList.push_back(FS_->ReadBlockRecord(in));
   }
   in >> IdList;
   DataSize = in.read<uint32>();
}

} // namespace PHeapFileSystem
