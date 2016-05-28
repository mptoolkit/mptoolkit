// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// pstream/pstream.cpp
//
// Copyright (C) 1999-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

// pstream.cpp
//
// created 1999-07-13
// Ian McCulloch
//
// non-template implementations for the pstream classes.

#include "common/trace.h"
#include "pstream.h"
#include <map>
#include <stdexcept>

namespace
{

using namespace PStream;
using namespace std;

// mapping to retrieve the factory class from a type_info or UID.

struct TypeInfoLess
{
   bool operator()(const std::type_info* t1, const std::type_info* t2) const
   {
      return t1->before(*t2);
   }
};

struct UID_Factory
{
   UID_Factory() : Factory(NULL) {}
   UID_Factory(const string& UID_, StreamableFactory* Factory_) : UID(UID_), Factory(Factory_) {}
   UID_Factory(const UID_Factory& f) : UID(f.UID), Factory(f.Factory) {}
 
   string UID;
   StreamableFactory* Factory;
};

typedef map<const std::type_info*, UID_Factory, TypeInfoLess> TypeInfoFactoryListType;
typedef map<string, StreamableFactory*>                  UIDFactoryListType;

TypeInfoFactoryListType TypeInfoFactoryList;
UIDFactoryListType      UIDFactoryList;

StreamableFactory* FindFactory(const string& UID)
{
   UIDFactoryListType::const_iterator I = UIDFactoryList.find(UID);
   if (I == UIDFactoryList.end()) return NULL;
   // else
   return I->second;
}

} // namespace

namespace PStream 
{

//
// opstream
//

opstream::~opstream()
{
   delete Buffer;
}

void opstream::flush()
{
   if (Buffer->buf_ptr() != Buffer->buf_begin() || Buffer->buf_ptr() == Buffer->buf_end())
   {
      Buffer->overflow();
   }
}

void opstream::set_format(int OtherFormat)
{
   if (format() == OtherFormat) return;

   generic_opstreambuf* New = generic_opstreambuf::CopyOpstreambuf(OtherFormat, *Buffer);
   delete Buffer;
   Buffer = New;
   BaseFormat = OtherFormat;
}

void opstream::put_id(id_type)
{
   PANIC("put_id() not overridden for this stream type - persistent objects not supported!");
}

void
opstream::push_version(VersionTag const& Tag, int Version)
{
   VersionNumbers[&Tag].push(Version);
}

void
opstream::pop_version(VersionTag const& Tag)
{
   DEBUG_CHECK(!VersionNumbers[&Tag].empty());
   VersionNumbers[&Tag].pop();
}

int
opstream::version_of(VersionTag const& Tag) const
{
   std::map<VersionTag const*, VersionStackType>::const_iterator I = VersionNumbers.find(&Tag);
   if (I == VersionNumbers.end() || I->second.empty())
      return Tag.default_version();
   return I->second.top();
}

//
// generic_opstreambuf
//

generic_opstreambuf::~generic_opstreambuf()
{
}

generic_opstreambuf::generic_opstreambuf(generic_opstreambuf const& Buffer)
  : BackPointer(Buffer.BackPointer), BufBegin(Buffer.BufBegin), BufEnd(Buffer.BufEnd),
     BufPtr(Buffer.BufPtr)
{
}

generic_opstreambuf::generic_opstreambuf(opstream* BackPointer_,
					 byte_type* Beg, byte_type* End, byte_type* Ptr)
  : BackPointer(BackPointer_), BufBegin(Beg), BufEnd(End), BufPtr(Ptr)
{
}

generic_opstreambuf*
generic_opstreambuf::CopyOpstreambuf(int Format, generic_opstreambuf const& Buffer)
{
   RANGE_CHECK(Format, 0, NumConcreteStreamClasses-1);
   switch (Format)
   {
      case 0 : return new opstreambuf<0>(Buffer);
      case 1 : return new opstreambuf<1>(Buffer);
      case 2 : return new opstreambuf<2>(Buffer);
      case 3 : return new opstreambuf<3>(Buffer);
	//      case 4 : return new opstreambuf<4>(Buffer);
   }
   return NULL;  // this can never be executed
}

generic_opstreambuf*
generic_opstreambuf::MakeOpstreambuf(int Format, opstream* BackPointer,
				     byte_type* Bufbegin, byte_type* BufEnd, byte_type* BufPtr)
{
   RANGE_CHECK(Format, 0, NumConcreteStreamClasses-1);
   switch (Format)
   {
      case 0 : return new opstreambuf<0>(BackPointer, Bufbegin, BufEnd, BufPtr);
      case 1 : return new opstreambuf<1>(BackPointer, Bufbegin, BufEnd, BufPtr);
      case 2 : return new opstreambuf<2>(BackPointer, Bufbegin, BufEnd, BufPtr);
      case 3 : return new opstreambuf<3>(BackPointer, Bufbegin, BufEnd, BufPtr);
	//      case 4 : return new opstreambuf<4>(BackPointer, Bufbegin, BufEnd, BufPtr);
   }
   return NULL;  // this can never be executed
}

//
// ipstream
//

ipstream::~ipstream()
{
   delete Buffer;
}

id_type ipstream::get_id()
{
   PANIC("get_id() not overridden for this stream type - persistent objects not supported!");
   return 0;
}

void ipstream::set_format(int OtherFormat)
{
   if (format() == OtherFormat) return;

   generic_ipstreambuf* New = generic_ipstreambuf::CopyIpstreambuf(OtherFormat, *Buffer);
   delete Buffer;
   Buffer = New;
   BaseFormat = OtherFormat;
}

void
ipstream::push_version(VersionTag const& Tag, int Version)
{
   VersionNumbers[&Tag].push(Version);
}

void
ipstream::pop_version(VersionTag const& Tag)
{
   DEBUG_CHECK(!VersionNumbers[&Tag].empty());
   VersionNumbers[&Tag].pop();
}

int
ipstream::version_of(VersionTag const& Tag) const
{
   std::map<VersionTag const*, VersionStackType>::const_iterator I = VersionNumbers.find(&Tag);
   if (I == VersionNumbers.end() || I->second.empty())
      return Tag.default_version();
   return I->second.top();
}

//
// generic_opstreambuf
//

generic_ipstreambuf::~generic_ipstreambuf()
{
}

generic_ipstreambuf::generic_ipstreambuf(generic_ipstreambuf const& Buffer)
  : BackPointer(Buffer.BackPointer), BufBegin(Buffer.BufBegin), BufEnd(Buffer.BufEnd),
     BufPtr(Buffer.BufPtr)
{
}

generic_ipstreambuf::generic_ipstreambuf(ipstream* BackPointer_,
					 byte_type const* Beg,
					 byte_type const* End,
					 byte_type const* Ptr)
  : BackPointer(BackPointer_), BufBegin(Beg), BufEnd(End), BufPtr(Ptr)
{
}

generic_ipstreambuf*
generic_ipstreambuf::CopyIpstreambuf(int Format, generic_ipstreambuf const& Buffer)
{
   RANGE_CHECK(Format, 0, NumConcreteStreamClasses-1);
   switch (Format)
   {
      case 0 : return new ipstreambuf<0>(Buffer);
      case 1 : return new ipstreambuf<1>(Buffer);
      case 2 : return new ipstreambuf<2>(Buffer);
      case 3 : return new ipstreambuf<3>(Buffer);
	//      case 4 : return new opstreambuf<4>(Buffer);
   }
   return NULL;  // this can never be executed
}

generic_ipstreambuf*
generic_ipstreambuf::MakeIpstreambuf(int Format, ipstream* BackPointer,
				     byte_type const* Bufbegin, 
				     byte_type const* BufEnd,
				     byte_type const* BufPtr)
{
   RANGE_CHECK(Format, 0, NumConcreteStreamClasses-1);
   switch (Format)
   {
      case 0 : return new ipstreambuf<0>(BackPointer, Bufbegin, BufEnd, BufPtr);
      case 1 : return new ipstreambuf<1>(BackPointer, Bufbegin, BufEnd, BufPtr);
      case 2 : return new ipstreambuf<2>(BackPointer, Bufbegin, BufEnd, BufPtr);
      case 3 : return new ipstreambuf<3>(BackPointer, Bufbegin, BufEnd, BufPtr);
	//      case 4 : return new opstreambuf<4>(BackPointer, Bufbegin, BufEnd, BufPtr);
   }
   return NULL;  // this can never be executed
}

// StreamableBase

void StreamableBase::ReadStream(ipstream& in)
{
   PANIC("ReadStream not supported for this object.");
}

void StreamableBase::WriteStream(opstream&) const
{
}

StreamableBase* StreamableBase::ConstructFromStream(ipstream& in)
{
   return NULL;
}

// global functions

void RegisterClass(const std::type_info& t, const string& UID, StreamableFactory* Factory)
{
   // check for duplicate ID's, firstly for a second registration of the same type_info
   string dUID = GetUID(t);
   if (!dUID.empty())
   {
      if (dUID != UID) 
      {
         delete Factory;
	 PANIC("") << "Duplicate UID for type " << t.name() << ": " << dUID << " and " << UID;
      }
      // else dUID == UID.  replace the existing factory with the new one
      delete UIDFactoryList[UID];
      UIDFactoryList.erase(UID);
   }

   // also make sure that there are no duplicated UID strings
   if (UIDFactoryList.find(UID) != UIDFactoryList.end()) 
   {
      PANIC("") << "Duplicate UID for type " << t.name() << ": " << dUID << " is already used." << UID;
   }

   UIDFactoryList[UID] = Factory;
   TypeInfoFactoryList[&t] = UID_Factory(UID, Factory);
}

string GetUID(const std::type_info& t)
{
   TypeInfoFactoryListType::const_iterator I = TypeInfoFactoryList.find(&t);
   if (I == TypeInfoFactoryList.end()) return "";
   // else
   return I->second.UID;
}

StreamableBase* StreamExtract(ipstream& in)
{
   string UID;
   in >> UID;

   StreamableFactory* Factory = FindFactory(UID);
   if (!Factory) // is the class registered?
   {
      // the UID of "$" is a special case, it denotes a NULL object.
      if (UID != "$")
      {
         PANIC("Expecting a class factory UID.");
      }
      return NULL;
   }

   return Factory->construct(in);
}

//opstream& operator<<(opstream& stream, const std::string& str)
//{
//   return stream << str.size() << make_proxy_container(str);
//}



} // namespace PStream
