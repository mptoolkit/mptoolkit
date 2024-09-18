// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// pheap/pagefile.cpp
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

#include "pagefile.h"

//#if defined(MULTITHREAD)
//#include "pagefile_mt.cpp"
//#else
#include "pagefile_st.cpp"
//#endif

namespace PHeapFileSystem
{

PageFile::PageFile()
  : Impl(new PageFileImpl())
{
}

void PageFile::create(size_t PageSize, std::string const& FileName, bool Unlink, bool AllowOverwrite)
{
   Impl->create(PageSize, FileName, Unlink, AllowOverwrite);
}

uint64 PageFile::open(std::string const& FileName, bool ReadOnly)
{
   return Impl->open(FileName, ReadOnly);
}

void PageFile::shutdown(bool Remove)
{
   Impl->shutdown(Remove);
}

void PageFile::persistent_shutdown(uint64 UserData)
{
   Impl->persistent_shutdown(UserData);
}

size_t PageFile::write(unsigned char const* Buffer)
{
   return Impl->write(Buffer);
}

unsigned char const* PageFile::read(size_t Page)
{
   return Impl->read(Page);
}

void PageFile::deallocate(size_t Page)
{
   Impl->deallocate(Page);
}

size_t PageFile::try_defragment(size_t Page)
{
   return Impl->try_defragment(Page);
}

BufferAllocator* PageFile::get_allocator() const
{
   return Impl->get_allocator();
}

std::string const& PageFile::name() const
{
   return Impl->name();
}

int PageFile::version() const
{
   return Impl->version();
}

size_t PageFile::get_page_size() const
{
   return Impl->get_page_size();
}

size_t PageFile::num_allocated_pages() const
{
   return Impl->num_allocated_pages();
}

size_t PageFile::num_free_pages() const
{
   return Impl->num_free_pages();
}

int PageFile::num_pages_written() const
{
   return Impl->num_pages_written();
}

int PageFile::num_pages_read() const
{
   return Impl->num_pages_read();
}

void PageFile::set_checkpoint_limit_kb(unsigned long Limit)
{
   Impl->set_checkpoint_limit_kb(Limit);
}

unsigned long PageFile::get_checkpoint_limit_kb() const
{
   return Impl->get_checkpoint_limit_kb();
}

void PageFile::Debug()
{
   Impl->Debug();
}

} // namespace PHeapFileSystem
