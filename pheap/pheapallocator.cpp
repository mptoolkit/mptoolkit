// -*- C++ -*- $Id$

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
