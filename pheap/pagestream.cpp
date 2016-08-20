// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// pheap/pagestream.cpp
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

#include "pagestream.h"

namespace PHeapFileSystem
{

//
// opagestream
//

opagestream::opagestream(PHeapFileSystem::FileSystem* FS_, int Format)
  : opstream(Format, NULL, NULL, NULL), MyFS(FS_), Committed(false)
{
}

opagestream::~opagestream()
{
   if (!Committed && !BufferList.empty())
   {
      WARNING("Destructing an opsagestream with uncommitted data.");
   }
}

void opagestream::flush()
{
   // Nothing to flush, since we require an explicit commit()
}

void opagestream::overflow()
{
   CHECK(this->buf_ptr() == this->buf_end() && !Committed);

   WriteBuffer NewBuffer = MyFS->allocate();
   BufferList.push_back(NewBuffer);
   this->set_buffer(NewBuffer.buffer(),
                    NewBuffer.buffer()+MyFS->get_page_size(),
                    NewBuffer.buffer()+4); // +4 reserved for linked list page_id
}

PageId opagestream::commit()
{
   CHECK(!Committed);
   Committed = true;
   // Walk the buffer list in reverse order, assembling the linked list of page numbers as we go.
   // The tail of the linked list has a 'next page' number of zero.
   PageId Page;
   while (!BufferList.empty())
   {
      this->set_buffer(BufferList.back().buffer(),
                       BufferList.back().buffer()+4,
                       BufferList.back().buffer());
      MyFS->write_page_id(*this, Page);
      Page = BufferList.back().get_page();
      BufferList.pop_back();
   }
   this->set_buffer(NULL, NULL, NULL);
   return Page;
}

//
// ipagestream
//

ipagestream::ipagestream(PHeapFileSystem::FileSystem* FS_, PageId FirstPage, int Format)
  : ipstream(Format, NULL, NULL, NULL), MyFS(FS_), NextPage(FirstPage)
{
}

ipagestream::~ipagestream()
{
}

void ipagestream::underflow()
{
   // if we've underflowed, there should be pages left to read
   //   CHECK(NextPage != PageId());

   //   std::cout << "ipagestream::underflow(), NextPage = " << NextPage << std::endl;
   Buf = ReadBuffer(NextPage);
   this->set_buffer(Buf.buffer(), Buf.buffer()+MyFS->get_page_size(), Buf.buffer());
   //std::cout << "Buffer is " << (void*) this->buf_begin() << ", " << (void*) this->buf_ptr()
   //        << ", " << (void*) this->buf_end() << std::endl;
   Pages.push_back(NextPage);
   NextPage = MyFS->read_page_id(*this);
}

std::list<PageId>
ipagestream::pages()
{
   return Pages;
}

void ipagestream::free()
{
   // make sure that we have a complete list of pages
   while (NextPage != PageId())
   {
      Buf = ReadBuffer(NextPage);
      this->set_buffer(Buf.buffer(), Buf.buffer()+MyFS->get_page_size(), Buf.buffer());
      Pages.push_back(NextPage);
      NextPage = MyFS->read_page_id(*this);
   }
   this->set_buffer(NULL, NULL, NULL);

   Pages.clear();
   Buf = ReadBuffer();
}

} // namespace PHeapFileSystem
