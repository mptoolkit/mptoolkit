// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// pheap/pheapallocator.h
//
// Copyright (C) 2002-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
  pheapallocator.h

  Maintains a list of BlockRecord objects, and a list of Id's.  These Id's are
  the pheap ID's of the nested objects.

  Created 2002-07-16 Ian McCulloch
*/

#if !defined(PHEAPALLOCATOR_H_DSKFJIORU3489UYR8YER89FH3489H89)
#define PHEAPALLOCATOR_H_DSKFJIORU3489UYR8YER89FH3489H89

#include "pheapfsblock.h"
#include "pstream/pstreamfwd.h"
#include <iostream>

namespace PHeapFileSystem
{

using PStream::id_type;
using namespace pheap;

class Descriptor
{
   private:
      typedef std::list<BlockRecord> BlockListType;
      typedef std::vector<id_type> IdListType;

   public:
      Descriptor() : DataSize(0) {}
      ~Descriptor() { }

      // constructs a descriptor by reading it from the given stream.  The desciptor is
      // associated with the specified filesystem.
      Descriptor(BlockFileSystem* FS_, PStream::ipstream& in);

      typedef BlockListType::const_iterator block_iterator;
      typedef IdListType::const_iterator id_iterator;
      
      size_t data_size() const { return DataSize; }

      // returns the number of data blocks
      size_t block_size() const { return BlockList.size(); }

      block_iterator block_begin() const { return BlockList.begin(); }
      block_iterator block_end() const { return BlockList.end(); }

      void append_block(BlockRecord const& Rec) { BlockList.push_back(Rec); DataSize += Rec.size(); }

      // returns the number of nested object ID's
      size_t id_size() const { return IdList.size(); }

      id_iterator id_begin() const { return IdList.begin(); }
      id_iterator id_end() const { return IdList.end(); }

      void append_id(id_type Id) { IdList.push_back(Id); }

      void Write(BlockFileSystem* FS_, PStream::opstream& out) const;

   private:
      BlockListType BlockList;
      IdListType IdList;
      size_t DataSize;

   friend std::ostream& operator<<(std::ostream& out, Descriptor const& d);
};

inline
std::ostream& operator<<(std::ostream& out, std::list<BlockRecord> const& r)
{
   std::copy(r.begin(), r.end(), std::ostream_iterator<BlockRecord>(out, ", "));
   return out;
}

inline
std::ostream& operator<<(std::ostream& out, std::vector<id_type> const& r)
{
   std::copy(r.begin(), r.end(), std::ostream_iterator<id_type>(out, ", "));
   return out;
}
inline
std::ostream& operator<<(std::ostream& out, Descriptor const& d)
{
   out << "Descriptor size=" << d.DataSize
       << "\nBlocks: " << d.BlockList 
       << "\nNested identifiers: " << d.IdList << '\n';
   return out;
}

} // namespace PHeapFileSystem

#endif
