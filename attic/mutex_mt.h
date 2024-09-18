// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/mutex_mt.h
//
// Copyright (C) 2004-2021 Ian McCulloch <ian@qusim.net>
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
  mutex_mt.h

  pthreads implementation of various locking primitives.

  Split from threads.h on 2002-07-22.
*/

#if !defined(MUTEX_MT_H_JH34UIY78934YR78Y89FU43P89H0)
#define MUTEX_MT_H_JH34UIY78934YR78Y89FU43P89H0

#include <semaphore.h>
#include "common/trace.h"

#include <iostream>

namespace pthread
{

// DMUTEX can be used as a quick-n-dirty serialization lock for debugging.
#if defined(NDEBUG)
#define DMUTEX
#else
#define DMUTEX pthread::rec_mutex::sentry Debug_Mutex_Lock_43YTDYT87348GUYDG(pthread::DebugMutex);
#endif

//
// some debugging check macros.
//
// CHECK_OWNER can be applied to a mutex, and
// verifies that the mutex is locked by the calling thread.
// This facility is only available for debug_mutex's.
//
// DEBUG_CHECK_OWNER works the same way, but
// expands to a no-operation if NDEBUG is defined.
//

#define CHECK_OWNER(x) CHECK(x.is_owner())
#define DEBUG_CHECK_OWNER(x) DEBUG_CHECK(x.is_owner())

class condition;

// some older glibc's have broken debug mutexes.

#if !defined(BROKEN_DEBUGMUTEX)
class debug_mutex
{
   public:
      debug_mutex();
      ~debug_mutex();

      void lock();
      void unlock();

      class sentry
      {
         public:
            sentry(debug_mutex& m) : M(m) { M.lock(); }
            ~sentry() { M.unlock(); }
         private:
            debug_mutex& M;
      };

      ::pthread_mutex_t const* get_mutex() const { return &M; }
      ::pthread_mutex_t* get_mutex() { return &M; }

      // is_owner() returns true if the mutex is locked by the calling thread,
      // or if the ownership of the mutex cannot be established.
      // If it is known that the mutex is unlocked or is locked by another thread,
      // is_owner() returns false.
      bool is_owner() { return Owner == pthread_self(); }

   private:
      debug_mutex(debug_mutex const&);            // not implemented
      debug_mutex& operator=(debug_mutex const&); // not implemented

      int Magic;   // a magic number set to MagicInit in the ctor, set to 0 in the dtor.  This detects
                   // attempts to use the mutex before construction or after destruction.
      ::pthread_mutex_t M;  // the actual mutex object
      pthread_t Owner;      // The owning thread id, or 0 if there is no owner.

      static int const MagicInit = 1837465;

   friend std::ostream& operator<<(std::ostream& out, debug_mutex const& D);
   friend void wait(condition& c, debug_mutex& m);
};
#endif

class fast_mutex
{
   public:
      fast_mutex() { ::pthread_mutex_init(&M, NULL);  }
      ~fast_mutex() { ::pthread_mutex_destroy(&M); }

      void lock();
      void unlock();

      class sentry
      {
         public:
            sentry(fast_mutex& m) : M(m) { M.lock(); }
            ~sentry() { M.unlock(); }
         private:
            fast_mutex& M;
      };

      ::pthread_mutex_t* get_mutex() { return &M; }

      // is_owner() returns true if the mutex is locked by the calling thread,
      // or if the ownership of the mutex cannot be established.
      // If it is known that the mutex is unlocked or is locked by another thread,
      // is_owner() returns false.
      bool is_owner() { return true; }

   private:
      fast_mutex(fast_mutex const&);            // not implemented
      fast_mutex& operator=(fast_mutex const&); // not implemented

      ::pthread_mutex_t M;
};

class rec_mutex
{
   public:
      rec_mutex();

      ~rec_mutex() { ::pthread_mutex_destroy(&M); }

      void lock();
      void unlock();

      class sentry
      {
         public:
            sentry(rec_mutex& m) : M(m) { M.lock(); }
            ~sentry() { M.unlock(); }
         private:
            rec_mutex& M;
      };

      ::pthread_mutex_t* get_mutex() { return &M; }

      // is_owner() returns true if the mutex is locked by the calling thread,
      // or if the ownership of the mutex cannot be established.
      // If it is known that the mutex is unlocked or is locked by another thread,
      // is_owner() returns false.
      bool is_owner() { return true; }

   private:
      rec_mutex(fast_mutex const&);            // not implemented
      rec_mutex& operator=(fast_mutex const&); // not implemented

      ::pthread_mutex_t M;
};

extern rec_mutex DebugMutex;

#if !defined(DEBUGMUTEX) || defined(BROKEN_DEBUGMUTEX)
typedef fast_mutex mutex;
#else
typedef debug_mutex mutex;
#endif

#if defined(BROKEN_DEBUGMUTEX)
typedef fast_mutex debug_mutex;
#endif

class condition
{
   public:
      condition();
      ~condition();

      pthread_cond_t* get_cond() { return &C; }

   private:
      pthread_cond_t C;

      condition(condition const&);            // not implemented
      condition& operator=(condition const&); // not implemented
};

// free functions for signalling and waiting on a condition.
// Note that we do not allow waiting on a recursive mutex.  This is
// allowed in pthreads but is extremely error-prone - the caller
// would have to be able to guarantee that the mutex is only locked
// ONCE before waiting.

void signal(condition& c);
void broadcast(condition& c);
#if !defined(BROKEN_DEBUGMUTEX)
void wait(condition& c, debug_mutex& m);
#endif
void wait(condition& c, fast_mutex& m);

class semaphore
{
   public:
      semaphore();
      explicit semaphore(unsigned int Initial);
      ~semaphore();

      void up();
      void down();  // this should be renamed to down_wait()

   private:
      sem_t Sem;

      semaphore(semaphore const&);            // not implemented
      semaphore& operator=(semaphore const&); // not implemented
};

class readers_writers
{
   public:
      readers_writers() : readers_pending(0), writer_count(0) {}

      void write_lock();
      void write_unlock();

      void read_lock();
      void read_unlock();

      class write_sentry
      {
         public:
            write_sentry(readers_writers& rw_) : rw(rw_) { rw.write_lock(); }
            ~write_sentry() { rw.write_unlock(); }
         private:
            write_sentry(write_sentry const&);            // not implemented
            write_sentry& operator=(write_sentry const&); // not implemented
            readers_writers& rw;
      };

      class read_sentry
      {
         public:
            read_sentry(readers_writers& rw_) : rw(rw_) { rw.read_lock(); }
            ~read_sentry() { rw.read_unlock(); }
         private:
            read_sentry(read_sentry const&);            // not implemented
            read_sentry& operator=(read_sentry const&); // not implemented
            readers_writers& rw;
      };

   private:
      int volatile readers_pending;  // number of pending or active readers
      int volatile writer_count;    // number of active writers (0 or 1)
      mutex M;
      condition LockFree;

   friend class read_sentry;
   friend class write_sentry;
};

//
// inlines
//

// fast_mutex

inline
void
fast_mutex::lock()
{
   if (int Err = ::pthread_mutex_lock(&M) != 0) throw_pthread_error("fast_mutex::lock(): ", Err);
}

inline
void
fast_mutex::unlock()
{
   if (int Err = ::pthread_mutex_unlock(&M) != 0) throw_pthread_error("fast_mutex::unlock(): ", Err);
}

// rec_mutex

inline
void
rec_mutex::lock()
{
   if (int Err = ::pthread_mutex_lock(&M) != 0) throw_pthread_error("rec_mutex::lock(): ", Err);
}

inline
void
rec_mutex::unlock()
{
   if (int Err = ::pthread_mutex_unlock(&M) != 0) throw_pthread_error("rec_mutex::unlock(): ", Err);
}

// condition

inline
condition::condition()
{
   if (int Err = pthread_cond_init(&C, NULL) != 0) throw_pthread_error("condition ", Err);
}

inline
condition::~condition()
{
   if (int Err = pthread_cond_destroy(&C) != 0) throw_pthread_error("condition ", Err);
}

inline
void signal(condition& c)
{
   //   DEBUGMSG("signal of condition " << (void*) &c << " by thread_id " << pthread_self());
   if (int Err = pthread_cond_signal(c.get_cond()) != 0) throw_pthread_error("condition ", Err);
}

inline
void broadcast(condition& c)
{
   //   DEBUGMSG("broadcast of condition " << (void*) &c << " by thread_id " << pthread_self());
   if (int Err = pthread_cond_broadcast(c.get_cond()) != 0) throw_pthread_error("condition ", Err);
}

inline
void wait(condition& c, fast_mutex& m)
{
   if (int Err = pthread_cond_wait(c.get_cond(), m.get_mutex()) != 0) throw_pthread_error("condition ", Err);
}

// semaphore

inline
semaphore::semaphore()
{
   sem_init(&Sem, 0, 0);
}

inline
semaphore::semaphore(unsigned int Initial)
{
   sem_init(&Sem, 0, Initial);
}

inline
semaphore::~semaphore()
{
   sem_destroy(&Sem);
}

inline
void semaphore::up()
{
   sem_post(&Sem);
}

inline
void semaphore::down()
{
   sem_wait(&Sem);
}

// readers_writers

inline
void readers_writers::write_lock()
{
   M.lock();
   while (readers_pending > 0 || writer_count > 0)
   {
      wait(LockFree, M);
   }
   DEBUG_CHECK(readers_pending == 0 && writer_count == 0);
   writer_count = 1;
   //   M.unlock();
}

inline
void readers_writers::write_unlock()
{
  //   M.lock();
   DEBUG_CHECK(writer_count == 1);
   writer_count = 0;
   broadcast(LockFree);
   M.unlock();
}

inline
void readers_writers::read_lock()
{
   M.lock();
   ++readers_pending;
   while (writer_count > 0)
   {
      wait(LockFree, M);
   }
   M.unlock();
}

inline
void readers_writers::read_unlock()
{
   M.lock();
   DEBUG_PRECONDITION(readers_pending > 0);
   if (--readers_pending == 0) signal(LockFree);
   M.unlock();
}

} // namespace pthread

#endif
