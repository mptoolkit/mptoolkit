// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// pheap/rawpagestream.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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

#include "rawpagestream.h"

namespace PHeapFileSystem
{

//
// orawpagestream
//

orawpagestream::orawpagestream(PHeapFileSystem::PageFile* MyFile_, int Format)
  : opstream(Format, NULL, NULL, NULL), MyPageFile(MyFile_), Committed(false)
{
}

orawpagestream::~orawpagestream()
{
   CHECK(Committed || BufferList.empty());
   if (this->buf_begin()) MyPageFile->deallocate_buffer(this->buf_begin());
}

void orawpagestream::flush()
{
   // Nothing to flush, since we require an explicit commit()
}

void orawpagestream::overflow()
{
   CHECK(this->buf_ptr() == this->buf_end());
   unsigned char* NewBuffer = MyPageFile->allocate_buffer();
   BufferList.push_back(NewBuffer);
   this->set_buffer(NewBuffer, NewBuffer+MyPageFile->get_page_size(), NewBuffer+4); // +4 reserved for linked list
}

size_t orawpagestream::commit()
{
   CHECK(!Committed);
   // Walk the buffer list in reverse order, assembling the linked list of page numbers as we go.
   // The tail of the linked list has a 'next page' number of zero.
   uint32 Page = 0;
   while (!BufferList.empty())
   {
      this->set_buffer(BufferList.back(), BufferList.back()+4, BufferList.back());
      this->write(Page);
      Page = MyPageFile->write(BufferList.back());
      BufferList.pop_back();
   }
   this->set_buffer(NULL, NULL, NULL);
   return Page;
}

//
// irawpagestream
//

irawpagestream::irawpagestream(PHeapFileSystem::PageFile* MyFile_, size_t FirstPage, int Format)
  : ipstream(Format, NULL, NULL, NULL), MyPageFile(MyFile_), NextPage(FirstPage)
{
}

irawpagestream::~irawpagestream()
{
   if (this->buf_begin())
   {
      MyPageFile->deallocate_buffer(this->buf_begin());
   }
}

void irawpagestream::underflow()
{
   if (this->buf_begin())
   {
      MyPageFile->deallocate_buffer(this->buf_begin());
      this->set_buffer(NULL, NULL, NULL);
   }

   // if we've underflowed, there should be pages left to read
   CHECK(NextPage != 0);

   unsigned char const* Buf = MyPageFile->read(NextPage);
   this->set_buffer(Buf, Buf+MyPageFile->get_page_size(), Buf);
   Pages.push_back(NextPage);
   NextPage = this->read<uint32>();
}

std::list<size_t>
irawpagestream::defer_free()
{
   if (this->buf_begin()) MyPageFile->deallocate_buffer(this->buf_begin());

   // make sure that we have a complete list of pages
   while (NextPage != 0)
   {
     unsigned char const* Buf = MyPageFile->read(NextPage);
     this->set_buffer(Buf, Buf+MyPageFile->get_page_size(), Buf);
     Pages.push_back(NextPage);
     NextPage = this->read<uint32>();
     MyPageFile->deallocate_buffer(Buf);
   }
   this->set_buffer(NULL, NULL, NULL);

   return Pages;
}

void irawpagestream::free()
{
   this->defer_free();

   // Walk the Pages list and deallocate the pages.
   while (!Pages.empty())
   {
      MyPageFile->deallocate(Pages.front());
      Pages.pop_front();
   }
}

} // namespace PHeapFileSystem
