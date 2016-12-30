// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// attic/threads.cc
//
// Copyright (C) 2015-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

namespace pthread
{

// helper for joining threads

template <class Result>
Result join_helper(pthread_t const& Thread)
{
   void* Ret;
   if (int ErrCode = ::pthread_join(Thread, &Ret) != 0)
   {
      throw std::runtime_error("cannot join thread.");
   }
   Result RetVal = *static_cast<Result*>(Ret);
   delete static_cast<Result*>(Ret);
   return RetVal;
}

template <>
void join_helper<void>(pthread_t const& Thread);

//
// thread
//

template <class Arg, class Result>
void* thread<Arg, Result>::dispatch_run()
{
   Result* Res = new Result(this->do_run(*this->Argument));
   delete Argument;
   Argument = NULL;
   return static_cast<void*>(Res);
}

template <class Arg>
void* thread<Arg, void>::dispatch_run()
{
   this->do_run(*this->Argument);
   delete Argument;
   Argument = NULL;
   return NULL;
}

template <class Result>
void* thread<void, Result>::dispatch_run()
{
   Result* Res = new Result(this->do_run());
   return static_cast<void*>(Res);
}

inline
void* thread<void, void>::dispatch_run()
{
   this->do_run();
   return NULL;
}

template <class Arg, class Result>
void thread<Arg, Result>::run(Arg const& arg)
{
   Argument = new Arg(arg);
   if (int ECode = ::pthread_create(&Thread, NULL,
                                    StartPThread,
                                    static_cast<void*>(static_cast<ThreadBase*>(this))) != 0)
   {
      throw std::runtime_error("not enough resources to create thread.");
   }
}

template <class Arg, class Result>
void thread<Arg, Result>::run_detached(Arg const& arg)
{
   Argument = new Arg(arg);
   if (int ECode = ::pthread_create(&Thread, NULL,
                                    StartPThread,
                                    static_cast<void*>(static_cast<ThreadBase*>(this))) != 0)
   {
      throw std::runtime_error("not enough resources to create thread.");
   }
   ::pthread_detach(Thread);
}

template <class Arg, class Result>
inline
Result
thread<Arg, Result>::join()
{
   return join_helper<Result>(this->get_thread());
}

// the partial specialization for no-arg threads

template <class Result>
void thread<void, Result>::run()
{
   if (int ECode = ::pthread_create(&Thread, NULL, StartPThread,
                                    static_cast<void*>(static_cast<ThreadBase*>(this))) != 0)
   {
      throw std::runtime_error("not enough resources to create thread.");
   }
}

inline
void thread<void, void>::run()
{
   if (int ECode = ::pthread_create(&Thread, NULL, StartPThread,
                                    static_cast<void*>(static_cast<ThreadBase*>(this))) != 0)
   {
      throw std::runtime_error("not enough resources to create thread.");
   }
}

template <class Result>
void thread<void, Result>::run_detached()
{
   if (int ECode = ::pthread_create(&Thread, NULL, StartPThread,
                                    static_cast<void*>(static_cast<ThreadBase*>(this))) != 0)
   {
      throw std::runtime_error("not enough resources to create thread.");
   }
   ::pthread_detach(Thread);
}

inline
void thread<void, void>::run_detached()
{
   if (int ECode = ::pthread_create(&Thread, NULL, StartPThread,
                                    static_cast<void*>(static_cast<ThreadBase*>(this))) != 0)
   {
      throw std::runtime_error("not enough resources to create thread.");
   }
   ::pthread_detach(Thread);
}


template <class Result>
inline
Result
thread<void, Result>::join()
{
   return join_helper<Result>(this->get_thread());
}

// free functions

template <class Arg, class Result>
inline
Result join(thread<Arg, Result>& Thread)
{
   return join_helper<Result>(Thread.get_thread());
}

inline
void thread<void, void>::join()
{
   return join_helper<void>(this->get_thread());
}

} // namespace pthread
