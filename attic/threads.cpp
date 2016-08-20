// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// attic/threads.cpp
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

#if !defined(MULTITHREAD)
#define MULTITHREAD
#endif

#include "common/trace.h"
#include "common/threads.h"
//#include "common/threadspecific_mt.h"
#include "common/mutex_mt.h"

namespace pthread
{

ThreadBase::ThreadBase()
{
}

ThreadBase::~ThreadBase()
{
}

#if defined(DEBUGMUTEX)
rec_mutex DebugMutex;
fast_mutex DebugMsgMutex;
#endif

extern "C"
void* StartPThread(void* arg)
{
   return static_cast<ThreadBase*>(arg)->dispatch_run();
}

template <>
void join_helper<void>(pthread_t const& Thread)
{
   void* Ret;
   if (int ErrCode = ::pthread_join(Thread, &Ret) != 0)
   {
      throw std::runtime_error("cannot join thread.");
   }
}

void throw_pthread_error(std::string Pre, int Err)
{
   throw std::runtime_error(Pre + ErrorMsg(Err));
}

char const* ErrorMsg(int Err)
{
   switch (Err)
   {
      case EINVAL  : return "not initialized";
      case EDEADLK : return "deadlock";
      case EBUSY   : return "busy";
      case EPERM   : return "not owner";
   }
   return "unknown error";
}

// thread_specific

#if 0
tsd_base::~tsd_base()
{
}

extern "C" void thread_specific_dtor(void* tsd)
{
   if (tsd != NULL) delete static_cast<tsd_base*>(tsd);
}
#endif

// mutex

#if !defined(BROKEN_DEBUGMUTEX)

std::ostream& operator<<(std::ostream& out, debug_mutex const& D)
{
   out << "debug_mutex " << (void*) &D << ", &pthread_mutex_t = " << D.get_mutex();
#if 0
#if defined(__linux)
   out << ", __m_count = " << D.get_mutex()->__m_count
       << ", __m_kind = " << D.get_mutex()->__m_kind
       << ", __m_lock.__status = " << D.get_mutex()->__m_lock.__status
       << ", __m_lock.__spinlock = " << D.get_mutex()->__m_lock.__spinlock
       << ", __m_owner = " << (void*) D.get_mutex()->__m_owner;
#endif
#endif
   return out;
}

debug_mutex::debug_mutex()
{
   //CHECK(Magic != MagicInit)(this); // this check is broken for multiprocessors!
   pthread_mutexattr_t myAttr;
   ::pthread_mutexattr_init(&myAttr);
   ::pthread_mutexattr_settype(&myAttr, PTHREAD_MUTEX_ERRORCHECK);
   ::pthread_mutex_init(&M, &myAttr);
   ::pthread_mutexattr_destroy(&myAttr);
   Owner = 0;
   Magic = MagicInit;
   //TRACE("mutex initialized")(this);
}

debug_mutex::~debug_mutex()
{
   CHECK(Magic == MagicInit)(this);
   Magic = 0;
   //TRACE("mutex destroyed")(this);
   ::pthread_mutex_destroy(&M);
}

void debug_mutex::lock()
{
   //   CHECK(Magic == MagicInit)(this);
   if (int Err = ::pthread_mutex_lock(&M) != 0) PANIC("debug_mutex::lock(): ")(ErrorMsg(Err));
   if (Owner != 0)
   {
      PANIC(" attempted lock by thread_id")(pthread_self())(Owner);
   }
   //   DEBUG_TRACE(*this << " locked by thread_id " << pthread_self());
   Owner = pthread_self();
}

void debug_mutex::unlock()
{
   //   CHECK(Magic == MagicInit)(this);
   //   DEBUGMSG(*this << " unlocked by thread_id " << pthread_self()
   //       << ", owner is thread_id " << Owner);
   pthread_t OldOwner = Owner;
   Owner = 0;
   if (int Err = ::pthread_mutex_unlock(&M) != 0)
   {
      //      PANIC("debug_mutex::unlock(): " << ErrorMsg(Err)
      //            << " this = " << (void*) this
      //            << " pthread_self = " << ::pthread_self()
      //            << " Owner = " << OldOwner);
   }
}

void wait(condition& c, debug_mutex& m)
{
   //   DEBUGMSG("wait on condition " << (void*) &c << " with " << m << " by thread_id " << pthread_self());
   //   CHECK(m.Magic == debug_mutex::MagicInit);
   CHECK(m.Owner == pthread_self());
   m.Owner = 0;
   if (int Err = pthread_cond_wait(c.get_cond(), m.get_mutex()) != 0)
      PANIC("condition")(ErrorMsg(Err));
   //   DEBUGMSG("condition triggered " << (void*) &c << " with " << m << " by thread_id " << pthread_self());
   m.Owner = pthread_self();
}

#endif

rec_mutex::rec_mutex()
{
   pthread_mutexattr_t myAttr;
   ::pthread_mutexattr_init(&myAttr);
   ::pthread_mutexattr_settype(&myAttr, PTHREAD_MUTEX_RECURSIVE);
   ::pthread_mutex_init(&M, &myAttr);
   ::pthread_mutexattr_destroy(&myAttr);
}

} // namespace pthread
