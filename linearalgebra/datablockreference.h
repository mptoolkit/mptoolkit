/* -*- C++ -*- $Id$

  datablockreference.h

  A double-indirection scheme to allow simultaneous reference and copy-on-write semantics for data blocks.

  Created 2004-04-20 Ian McCulloch

  Idea:
  
  struct IndirectBlock
  {
    int* RefCount;
    DataBlock<T>* Data;
  };

  The initial value & subsequent references share the RefCount and Data members
  of the IndirectBlock.  When the RefCount goes to zero, the Data member is
  deleted.

  A value-semantic copy is done by constructing a new IndirectBlock
  (with initial reference count of 1), and a new Data member, but this
  initially has the same value as the old (ie, Data new DataBlock<T>(*Old.Data).
  Thus, they initially share representations.

  A copy-on-write is done by modifying the DataBlock itself, ie
  (*Data) = DataBlock<T>(size);

  In the case that Data->is_shared() is false, the copy-on-write can be elided.

  The DataBlock reference count is a count of the number of value-semantic
  *groups* that are sharing the representation.

  The BlockValue::RefCount reference count is a count of the number of
  reference-semantic members of a group.
*/

#if !defined(DATABLOCKREFERENCE_H_CSHJ4389RU89UFG8743Y98H4P3)
#define DATABLOCKREFERENCE_H_CSHJ4389RU89UFG8743Y98H4P3

#include "datablock.h"
#include "common/poolallocator.h"

namespace Private
{

template <typename T, typename BlockHeader = NoHeader, typename LocalHeader = NoHeader>
class BlockIndirection : private LocalHeader
{
   public:
      // compiler-generated copy ctor, copy assignment and dtor are OK
 
      BlockIndirection() : RefCount(1) {}

      template <typename U>
      BlockIndirection(int RefCount_, DataBlock<U, BlockHeader> const& Data_) : RefCount(RefCount_), Data(Data_) {}

      template <typename U>
      BlockIndirection(int RefCount_, DataBlock<U, BlockHeader> const& Data_, LocalHeader const& HLocal) 
	: LocalHeader(HLocal), RefCount(RefCount_), Data(Data_) {}

      BlockIndirection(int RefCount_, size_t Size) : RefCount(RefCount_), Data(Size) {}

      BlockIndirection(int RefCount_, size_t Size, T const& Fill) : RefCount(RefCount_), Data(Size, Fill) {}

      BlockIndirection(int RefCount_, size_t Size, BlockHeader const& H, 
		       LocalHeader const& HLocal = LocalHeader()) 
	: LocalHeader(HLocal), RefCount(RefCount_), Data(Size, H) {}

      BlockIndirection(int RefCount_, BlockHeader const& H, 
		       LocalHeader const& HLocal = LocalHeader()) 
	: LocalHeader(HLocal), RefCount(RefCount_), Data(H) {}

      BlockIndirection(int RefCount_, size_t Size, T const& Fill, BlockHeader const& H, 
		       LocalHeader const& HLocal = LocalHeader()) 
	: LocalHeader(HLocal), RefCount(RefCount_), Data(Size, Fill, H) {}

      DataBlock<T, BlockHeader> const& data() const { return Data; }
      DataBlock<T, BlockHeader>& data() { return Data; }

      void add_reference() { TRACE_DATABLOCK("add_reference")(this)(RefCount.value())(Data.get()); ++RefCount; }
      void sub_reference() { TRACE_DATABLOCK("sub_reference")(this)(RefCount.value())(Data.get()); if (--RefCount == 0) delete this; }

      void cow() { Data.cow(); }

      T* get() const { return Data.get(); }

      size_t size() const { return Data.size(); }

      BlockHeader const& data_header() const { return Data.header(); }
      BlockHeader& data_header() { return Data.header(); }

      LocalHeader const& local_header() const { return static_cast<LocalHeader const&>(*this); }
      LocalHeader& local_header() { return static_cast<LocalHeader&>(*this); }

      bool is_shared() const { return RefCount.value() > 1; }

      // since BlockIndirection's are allocated using new, we had  better make it fast.
      IMPLEMENT_POOL_ALLOCATOR_NEW(BlockIndirection)

   private:
      //      BlockIndirection(); // not implemented
      AtomicRefCount RefCount;
      DataBlock<T, BlockHeader>   Data;
};

} // namespace Private


template <typename T, typename Header = NoHeader, typename LocalHeader = NoHeader>
class BlockReference
{
   public:
      typedef Private::BlockIndirection<T, Header, LocalHeader> Indirector;

      BlockReference() : Block(new Indirector) {}

      explicit BlockReference(size_t Size) 
	: Block(new Indirector(1, Size)) {}

      explicit BlockReference(size_t Size, T const& Fill) 
	: Block(new Indirector(1, Size, Fill)) {}

      explicit BlockReference(size_t Size, Header const& H, LocalHeader const& HL = LocalHeader()) 
	: Block(new Indirector(1, Size, H, HL)) {}

      explicit BlockReference(Header const& H, LocalHeader const& HL = LocalHeader()) 
	: Block(new Indirector(1, H, HL)) {}

      explicit BlockReference(size_t Size, T const& Fill, Header const& H, LocalHeader const& HL = LocalHeader()) 
	: Block(new Indirector(1, Size, Fill, H, HL)) {}

      BlockReference(BlockReference const& v) : Block(v.Block)
      {
	 Block->add_reference();
      }

      ~BlockReference()
      {
         // 	 TRACE_DATABLOCK("~BlockReference()")(Block);
	 Block->sub_reference();
      }

      BlockReference& operator=(BlockReference const& v)
      {
	 v.Block->add_reference();
	 Block->sub_reference();
	 Block = v.Block;
	 return *this;
      }

      // returns value-semantic copies
      BlockReference copy() const 
        { TRACE_DATABLOCK("BlockReference::copy()")(Block); return BlockReference(1, Block->data(), this->local_header()); }
  //      BlockReference<T const, Header, LocalHeader> const_copy() const 
  //        { return BlockReference<T const, Header, LocalHeader>(1, Block->data(), this->local_header()); }

      T* get() const { return Block->get(); }

      size_t size() const { return Block->size(); }

      Header const& data_header() const { return Block->data_header(); }
      Header& data_header() { return Block->data_header(); }

      LocalHeader const& local_header() const { return Block->local_header(); }
      LocalHeader& local_header() { return Block->local_header(); }

      bool is_shared() const { return Block->is_shared(); }

      DataBlock<T, Header>& data() { return Block->data(); }
      DataBlock<T, Header> const& data() const { return Block->data(); }

      void cow() { Block->cow(); }

   private:
      BlockReference(int RefCount, DataBlock<T, Header> const& Data, LocalHeader const& LocalH) 
	: Block(new Indirector(RefCount, Data, LocalH)) {}

      Indirector* Block;

   friend class BlockReference<T const, Header, LocalHeader>;
};

template <typename T, typename Header, typename LocalHeader>
class BlockReference<T const, Header, LocalHeader>
{
   public:
      typedef Private::BlockIndirection<T const, Header, LocalHeader> Indirector;

      BlockReference() : Block(new Indirector) {}

      BlockReference(BlockReference<T const, Header, LocalHeader> const& v) : Block(v.Block)
      {
	 Block->add_reference();
      }

      BlockReference(BlockReference<T, Header, LocalHeader> const& v) : Block(v.Block)
      {
	 Block->add_reference();
      }

      ~BlockReference()
      {
	 Block->sub_reference();
      }

      BlockReference& operator=(BlockReference<T const, Header> const& v)
      {
	 v.Block->add_reference();
	 Block->sub_reference();
	 Block = v.Block;
	 return *this;
      }

      BlockReference& operator=(BlockReference<T, Header> const& v)
      {
	 v.Block->add_reference();
	 Block->sub_reference();
	 Block = v.Block;
	 return *this;
      }

      T const* get() const { return Block->get(); }

      BlockReference<T const, Header> copy() const 
      { return BlockReference<T const, Header, LocalHeader>(1, Block->data(), this->local_header()); }
      BlockReference<T const, Header> const_copy() const 
      { return BlockReference<T const, Header>(1, Block->data(), this->local_header()); }

      size_t size() const { return Block->size(); }

      Header const& data_header() const { return Block->data_header(); }
      Header& data_header() { return Block->data_header(); }

      LocalHeader const& local_header() const { return Block->local_header(); }
      LocalHeader& local_header() { return Block->local_header(); }

      DataBlock<T, Header> const& data() const { return Block->data(); }

   private:
      BlockReference(int RefCount, DataBlock<T, Header> const& Data) 
	: Block(new Private::BlockIndirection<T, Header>(RefCount, Data)) {}

      Private::BlockIndirection<T, Header, LocalHeader>* Block;

   friend class BlockReference<T, Header, LocalHeader>;
};

#endif
