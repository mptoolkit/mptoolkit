// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// pheap/rawpagestream.h
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

/*
  A low-level pstream for PageFile's.

  Created 2004-01-16 Ian McCulloch
*/

#if !defined(RAWPAGESTREAM_H_HJYRUY43985Y8FH4H383298)
#define RAWPAGESTREAM_H_HJYRUY43985Y8FH4H383298

#include "pstream/pstream.h"
#include "pagefile.h"

namespace PHeapFileSystem
{

class orawpagestream : public PStream::opstream
{
   public:
      orawpagestream(PHeapFileSystem::PageFile* MyFile_, int Format = PStream::format::XDR);

      ~orawpagestream();

      virtual void flush();

      // commits the buffers to the page file, and returns the first page number in the linked list.
      // This can be supplied to the irawpagestream constructor to read back the same stream.
      size_t commit();

   private:
      virtual void overflow();

      PHeapFileSystem::PageFile* MyPageFile;
      std::list<unsigned char*> BufferList;    // buffer list is built up in reverse order
      bool Committed;
};

class irawpagestream : public PStream::ipstream
{
 public:
      irawpagestream(PHeapFileSystem::PageFile* MyFile_, size_t FirstPage, int Format = PStream::format::XDR);

      ~irawpagestream();

      // deallocates the pages, so they can be reused by the PageFile.
      // This must be done explicitly,
      // simply destructing the irawpagestream is not enough.
      void free();

      // deallocates memory, but instead of deallocating file pages,
      // return the list of used pages.  This is to support transactions,
      // so that the metadata pages are not overwritten immediately.
      std::list<size_t> defer_free();

   private:
      virtual void underflow();

      PHeapFileSystem::PageFile* MyPageFile;
      size_t NextPage;                          // The next page to read on an underflow
      std::list<size_t> Pages;                  // the list of pages that we have read
};

} // namespace PHeapFileSystem

#endif
