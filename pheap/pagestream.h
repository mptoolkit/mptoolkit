// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// pheap/pagestream.h
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

/*
  A low-level pstream for FileSystem's.

  Created 2004-01-16 Ian McCulloch

  This defines an [i|o]pagestream that writes to a linked list of PageId's,
  interfacing directly with a PHeapFileSystem::FileSystem.

  The opagestream class defines a function commit() which assembles the linked list
  and returns the PageId of the head of the list.  (If the opagestream is destructed
  before calling commit(), the pages will be reclaimed safely, although this
  is almost certainly a bug and a warning is issued).

  The ipagestream constructor takes the PageId that is the head of the list.

*/

#if !defined(PAGESTREAM_H_HJYRUY43985Y8FH4H383298)
#define PAGESTREAM_H_HJYRUY43985Y8FH4H383298

#include "pstream/pstream.h"
#include "pheapfsv4.h"

namespace PHeapFileSystem
{

class opagestream : public PStream::opstream
{
   public:
      opagestream(PHeapFileSystem::FileSystem* FS_, int Format = PStream::format::XDR);

      ~opagestream();

      virtual void flush();

      // commits the buffers to the page file, and returns the first page number in the linked list.
      // This can be supplied to the ipagestream constructor to read back the same stream.
      PageId commit();

   private:
      virtual void overflow();

      PHeapFileSystem::FileSystem* MyFS;
      std::list<WriteBuffer> BufferList;    // buffer list is built up in reverse order
      bool Committed;
};

class ipagestream : public PStream::ipstream
{
 public:
      ipagestream(PHeapFileSystem::FileSystem* FS_, PageId FirstPage, int Format = PStream::format::XDR);

      ~ipagestream();

      // returns the list of all pages used so far by this stream
      std::list<PageId> pages();

      // deallocates the pages, so they can be reused by the filesystem.  This must be done explicitly,
      // simply destructing the ipagestream is not enough.
      void free();

   private:
      virtual void underflow();

      PHeapFileSystem::FileSystem* MyFS;
      PageId NextPage;                          // The next page to read on an underflow
      std::list<PageId> Pages;                  // the list of pages that we have read
      ReadBuffer Buf;
};

} // namespace PHeapFileSystem

#endif
