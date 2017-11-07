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

// thread pool design:
//
// we try to keep one thread spinning.  (maybe if there are no queues
// we could sleep.)  The thread goes through the list of queues
// and tries to find a queue that has a runnable task.  To do that,
// we first grab a spinlock from the queue (eg using an atomic_flag,
// see http://en.cppreference.com/w/cpp/atomic/atomic_flag).  We then
// see if the queue is runnable.  If not, then drop the spinlock and try the
// next queue.
// If we can run a task then signal the thread pool, which will wake up another
// thread.  Suggest that we only ever have one thread reading the queue list
// at a time.  In fact it can be protected by a mutex.  If the thread finds
// a queue that ia defunct then it can remove it immediately.  So threads in the
// thread pool will be blocked on that mutex, unless they are actually running
// or can get the mutex and search the queue list.
//
// Adding tasks to a queue will typically be done by a different thread.  So
// this needs to be threadsafe.  It can be a simple std::list, queue or boost
// lock-free queue.

namespace async
{

class event;

struct task
{
   std::function<bool()> can_run;
   std::function<void()> task;
};

class task
{
   public:
      task() {}
      virtual ~task() = default;

      virtual bool can_run() const = 0;
      virtual void run() = 0;
};

class task_simple : public task
{
   public:
      task_simple(std::function<void()> const& f_) : f(f_) {}
      virtual bool can_run() const { return true; }
      virtual void run() { f(); }
   private:
      std::function<void()> f;
};

// we store the events via pointer, which is effectively a double
// indirection.  There is a danger however that we end up with a
// dangling pointer.  Maybe we need to incorporate the double indirection
// into the cuda::event itself?
// But we should be able to pass by value here?!?! No, because we
// might schedule some task before the event has been recorded.  Although
// the user must ensure that there is enough synchronization to ensure that
// the event is recorded before the task executes, the task object itself
// may be created arbitarily far before the event is initialized.
//
// Note that we cannot schedule an event before it is recorded because that
// would violate normal ordering.  That is, we need to do a thread synchronization
// surrounding blocking calls.  Typical scenario:
//
// 1. A = B*C
// 2. D = eigenvalues(A)
// 3. E = B*A
// If we run asynchronously then step 3 cannot be scheduled until after
// step 2 has completed, and we cannot use cuda events to synchronize this.
// When we schedule step 3, we need to keep in mind that the event is not yet
// recorded.  The async gpu matrix needs to store the event by reference.
// Is it sufficient to do reference-counted events?  Probably not needed -
// we don't need to reference count memory buffers etc because the destructor
// is itself an async event so normal scheduling will ensure that the destructor
// doesn't run early (this means that the task must take ownership via move
// semantics when destroying).  Is all of this solved by simply storing a
// pointer to a gpu matrix?  No, we must need some amount of additional syncronization,
// consider the case where some queue (thread) is running and waiting for some
// other queue.  eg
//
// Queue 1:
// A = B;
// A *= 2;
//
// Queue 2:
// C = B*A;
//
// Queue 2 must be ordered with respect to the operations in queue 1,
// in order to get the right semantics and get the correct cuda event pointer.
// That is, C waits for A, so A must have already recorded the event.  That is,
// we cannot schedule the cuda kernel for C until after we have scheduled the
// appropriate operations in queue 1.  That is, the cuda calls need to
// be serialized appropriately.
// How efficient does the queue synchronization need to be?  We can use a similar
// mechanism, but it should be fairly efficient.  Is polling good enough?
//
// If we want to have the thread hold a mutex on the task list while running,
// we could have a separate list for newly added tasks.

class task_wait_gpu : public task
{
   public:
      explicit task_wait_gpu(std::function<void()> const& f_,
                             std::initializer_list<cuda::event const&> events)
         : f(f_)
      {
         event_list.reserve(events.size());
         for (auto const& e : events)
         {
            event_list.push_back(e);
         }
      }

      virtual void can_run() const
      {
         auto I = event_list.begin();
         while (I != event_list.end())
         {
            if (I->is_complete())
            {
               auto J = I;
               ++I;
               event_list.erase(J);
            }
            else
            {
               return false;
            }
         }
         return true;
      }

      void run()
      {
         f();
      }
   private:
      std::function<void()> f;
      mutable std::vector<cuda::event_ref> event_list;
};

// a task to asynchronously destroy a gpu_buffer.  We move the buffer into
// the task, which runs once the stream has finshed.
class task_destroy_gpu_buffer : public task
{
   public:
      template <typename T>
      explicit task_destroy_gpu_buffer(gpu_buffer<T>&& Buf)
      {
         std::tie(Ptr, ByteSize, Stream, Arena) = std::move(Buf).move_buffer();
      }

      virtual void can_run() const
      {
         return !Stream.is_running();
      }

      void run()
      {
         Stream.synchronize();
         Arena.free(Ptr, ByteSize);
      }

   private:
      void* Ptr;
      std::size_t ByteSize;
      mutable cuda::stream Stream;
      blas::arena Arena;
};

// A queue is a list of tasks. Typically a queue is constructed on the heap,
// and destroys itself with the self_destruct() method.

// A task is
// If the result is true, then the task was executed and can be removed
// from the queue.  If the result is false, then the task couldn't be executed
// yet (eg, bdecause it would block.
using task_type = std::function<bool()>;

class queue
{
   public:
      queue();
      queue(Queue const&) = delete;
      queue(queue&& Other);

      // Queue destructor blocks until the task list is empty
      ~queue();

      // The normal mode is to allocate a queue on the heap and
      // let the queue control its own lifetime.  This function
      // requires that the queue was allocated with operator new,
      // and schedules the deallocation of the queue once the
      // task list is empty.
      void self_destruct();

      // records an event at the end of the queue
      event record() const;

      // force the queue to wait for event e to complete
      void wait(event const& e);

      void add_task(task_type const& t);

      void add_task(std::function<void()> const& t)
      {
         this->add_task([t]->bool { t(); return true; });
      }

   private:
      std::queue<std::unique_ptr<task>> TaskList;
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
         std::atomic_thread_fence(std::memory_order_release);
         std::atomic_store_explicit(trigger_, true, std::memory_order_release);
      }

      // clears the event - equivalent to *this = event()
      void clear();

      // returns true if this is a null event
      bool is_null() const { return event_ == nullptr; }

      // returns true if work has been sucessfully completed
      bool is_complete() const
      {
         bool x = std::atomic_load_explicit(trigger_, std::memory_order_acquire);
         if (x)
         {
            std::atomic_thread_fence(std::memory_order_acquire);
         }
         return x;
      }


   private:
      std::atomic<bool>* trigger_;
      shared_counter count_;
};

} // namespace async

#endif
