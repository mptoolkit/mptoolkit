// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// attic/threads.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER
/*
  threads.h

  First attempt at a wrapper around pthreads.

  Created 2000-10-06 Ian McCulloch

  TODO: there is no error checking so far.  Some exceptions would be nice.
        Need to document which pthreads functions are safe to call direct.  
        eg pthread_exit() is safe.
	These functions should then be included in the interface.
*/

#if !defined(THREADS_H_FDHJSF89RF7Y34TER89U3U8989U)
#define THREADS_H_FDHJSF89RF7Y34TER89U3U8989U

#include </usr/include/pthread.h>
#if defined(__DECCXX)
#define SEM_NAME_MAX 255
#endif
#include <sched.h>
#include <errno.h>
#include <stdexcept>
#include <string>

#if !defined(MULTITHREAD)
#error "threads.h cannot be used without -DMULTITHREAD!";
#endif

namespace pthread
{

class pthread_attributes
{
   public:
      pthread_attributes() { pthread_attr_init(&Attr); }

      ~pthread_attributes() { pthread_attr_destroy(&Attr); }

   private:
      pthread_attr_t Attr;
};

// returns an error message given the corresponding error code; only pthreads-related
// errors are used.  There is probably a libc function to do this?
char const* ErrorMsg(int Err);

void throw_pthread_error(std::string Pre, int Err);

extern "C"
void* StartPThread(void* arg);

//
// thread
//
// base class for creating a thread function.  This is basically an
// asynchronous version of unary_function<Arg, Result>,
// where the argument is passed to run(), and the result
// is returned by join().
//

extern "C" void* StartPThread(void* arg);  // for internal use

class ThreadBase
{
   public:
      ThreadBase();
      virtual ~ThreadBase() = 0;
   private:
      ThreadBase(ThreadBase const&); // not implemented
      ThreadBase& operator=(ThreadBase const&); // not implemented
      virtual void* dispatch_run() = 0;
   friend void* StartPThread(void* arg);
};

template <class Arg, class Result>
class thread : private ThreadBase
{
   public:
      typedef Arg    argument_type;
      typedef Result result_type;

      thread() {}

      // starts the thread, returns immediately.
      void run(argument_type const& arg);

      // a hack to create a detached thread, the correct solution is to make a detached_thread class.
      void run_detached(argument_type const& arg);

      result_type join();

      pthread_t const& get_thread() const { return Thread; }

   private:
      pthread_t Thread;

      argument_type* Argument;
      void* dispatch_run();

      // the actual thread function.
      // Define this in the derived class to run the threaded code.
      virtual result_type do_run(argument_type const& arg) = 0;
};

template <class Arg>
class thread<Arg, void> : private ThreadBase
{
   public:
      typedef Arg    argument_type;
      typedef void result_type;

      thread() {}

      // starts the thread, returns immediately.
      void run(argument_type const& arg);

      // a hack to create a detached thread, the correct solution is to make a detached_thread class.
      void run_detached(argument_type const& arg);

      result_type join();

      pthread_t const& get_thread() const { return Thread; }

   private:
      pthread_t Thread;

      argument_type* Argument;
      void* dispatch_run();

      // the actual thread function.
      // Define this in the derived class to run the threaded code.
      virtual result_type do_run(argument_type const& arg) = 0;
};

// partial specialization for argument_type of void.  This
// changes the interface.
template <class Result>
class thread<void, Result> : private ThreadBase
{
   public:
      typedef void   argument_type;
      typedef Result result_type;

      thread() {}

      // starts the thread, returns immediately
      void run();

      // a hack to create a detached thread, the correct solution is to make a detached_thread class.
      void run_detached();

      result_type join();

      pthread_t const& get_thread() const { return Thread; }

   private:
      pthread_t Thread;

      void* dispatch_run();

      // the actual thread function.
      // Define this in the derived class to run the threaded code.
      virtual result_type do_run() = 0;
};

template <>
class thread<void, void> : private ThreadBase
{
   public:
      typedef void   argument_type;
      typedef void result_type;

      thread() {}

      // starts the thread, returns immediately
      void run();

      // a hack to create a detached thread, the correct solution is to make a detached_thread class.
      void run_detached();

      result_type join();

      pthread_t const& get_thread() const { return Thread; }

   private:
      pthread_t Thread;

      void* dispatch_run();

      // the actual thread function.
      // Define this in the derived class to run the threaded code.
      virtual result_type do_run() = 0;
};

// a free-standing version of join.
template <class Arg, class Result>
Result join(thread<Arg, Result>& Thread);

// yields the processor for the remainder of the timeslice.
inline void yield()
{
   ::sched_yield();
}

inline
pthread_t self()
{
   return ::pthread_self();
}

} // namespace pthread

#include "threads.cc"

#endif
