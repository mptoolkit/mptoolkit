/* -*- C++ -*- $Id$

  Stream classes to support writing to a BlockFileSystem.

  Originally created in antiquity, updated to use the new BlockFileSystem classes
  on 2002-07-18.

  Updated to the new pstreams on 2004-01-22

  2006-01-19: Bug fixed, formatting now works.  The ipheapstream reads the format byte
  automatically on construction.  For opheapstream, the format byte must be explicitly
  written.  (The bug was that automatically writing the format byte on construction
  leads to a second format entry in the case where an object is duplicated by copying
  the entire buffer.)
*/

#if !defined(PHEAPSTREAM_H_SDG236325SDFGHSVTY4RWT6SDAVG46T7SFDG)
#define PHEAPSTREAM_H_SDG236325SDFGHSVTY4RWT6SDAVG46T7SFDG

#include "pstream/pstream.h"
#include "pheapallocator.h"

namespace PHeapFileSystem
{

using namespace PStream;

class opheapstream : public PStream::opstream
{
   public:
      opheapstream(BlockFileSystem* FS_, Descriptor* Descriptor);
      ~opheapstream();

      void put_format(int Format) { this->PStream::opstream::put_format(Format); }

      virtual void put_id(id_type ID);

   private:
      virtual void overflow();

      Descriptor* MyDescriptor;
      OutputBuffer* OutBuffer;
      BlockFileSystem* MyFS;
};

class ipheapstream : public PStream::ipstream
{
   public:
      ipheapstream(BlockFileSystem* FS_, Descriptor const* Descriptor_);

      ~ipheapstream();

      virtual id_type get_id();

   private:
      virtual void underflow();

      Descriptor const* MyDescriptor;
      Descriptor::block_iterator CurrentBlock;
      Descriptor::id_iterator CurrentId;
      InputBuffer InBuffer;
      BlockFileSystem* MyFS;
};

} // namespace PHeapFileSystem

#endif




