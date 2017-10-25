// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// cuda/cuda.h
//
// Copyright (C) 2017 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#if !defined(MPTOOLKIT_ASYNC_ASYNC_H)
#define MPTOOLKIT_ASYNC_ASYNC_H

namespace async
{

class event;

// A queue is a list of tasks.  A task is a function that returns a bool.
// If the result is true, then the task was executed and can be removed
// from the queue.  If the result is false, then the task couldn't be executed
// yet (eg, bdecause it would block.
using task_type = std::function<bool()>;

class queue
{
   public:
      // records an event at the end of the queue
      event record() const;

      // force the queue to wait for event e to complete
      void wait(event const& e);

      void add_task(task_type const& t);

   private:
      std::queue<task_type> TaskList;
      mutable event Sync;
};

struct TaskWaitEvent
{
   TaskWaitEvent(event const& e_) : e(e_) {}

   bool operator() const
   {
      return e.is_complete();
   }

   event e;
};

void
queue::wait(event const& e)
{
   TaskList.push_back(TaskWaitEvent(e));
}

void
queue::add_task(queue::task_type const& t)
{
   Senc.clear();
   TaskList.push_back(t);
}

struct TaskTriggerEvent
{
   TaskTriggerEvent(event const& e_) : e(e_) {}

   bool operator() const
   {
      e.trigger();
      return true;
   }

   event e;
};


event
queue::record() const
{
   if (Sync.is_null())
   {
      TaskList.push_back(TaskTriggerEvent(Sync));
   }
}

class event
{
   public:
      event();
      event(event const& other);
      event(event&& other);
      event& operator=(event const& other);
      event& operator=(event&& other);
      ~event();

      // precondition: is_null().  sets the event to false
      void initialize();

      // triggers the event
      void trigger()
      {
         *trigger_ = true;
         std::atomic_thread_fence(std::memory_order_acquire);

      // clears the event - equivalent to *this = event()
      void clear();

      // returns true if this is a null event
      bool is_null() const { return event_ == nullptr; }

      // returns true if work has been sucessfully completed
      bool is_complete() const;

   private:
      std::atomic<bool>* trigger_;
      shared_counter count_;
};

template <typename T>
class async_matrix : public blas::BlasMatrix<typename T::value_type, async_matrix<T>, async_tag>
{
   public:
      using base_type          = T;
      using value_type         = typename base_type::value_type;
      using storage_type       = typename base_type::storage_type;
      using const_storage_type = typename base_type::const_storage_type;

      async_matrix(async_matrix&& Other) = default;

      template <typename U>
      async_matrix(blas::MatrixRef<value_type, U, async_tag> const& E)
         : async_matrix(E.rows(), E.cols())
      {
         assign(*this, E.as_derived());
      }

      // check here that the queue is empty, and there are no other objects waiting for us.
      // we can possibly work around other objects waiting for us with reference counting.
      ~async_matrix();

      void wait_for(async::queue const& x)
      {
         Queue.wait(x.record());
      }

   private:
      base_type Base;
      async::queue Queue;
};

// BLAS-like functions

template <typename T, typename U, typename V, typename W>
inline
void
gemv(T alpha, blas::BlasMatrix<T, U, async_tag> const& A,
     blas::BlasVector<T, V, async_tag> const& x, T beta,
     async_matrix<W>& y)
{
   DEBUG_CHECK_EQUAL(A.cols(), x.size());
   DEBUG_CHECK_EQUAL(y.size(), A.rows());
   y.wait_for(A.queue());
   y.wait_for(x.queue());
   y.add_task([&]{gemv(alpha, A.base(), x.base(), beta, y.base()); return true;});
}

} // namespace async

#endif
