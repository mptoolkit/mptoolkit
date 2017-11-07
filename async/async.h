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
         Srream.synchronize();
         Arena.free(Ptr, ByteSize);
      }

   private:
      void* Ptr;
      std::size_t ByteSize;
      mutable cuda::stream Stream;
      blas::arena Arena;
};

// A queue is a list of tasks.  A task is a function that returns a bool.
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
      // order here is important, we need to destroy the queue before the base
      // matrix, so construct in order Base, Queue, destroy in opposite order.
      base_type Base;
      async::queue Queue;
};

// BLAS-like functions

// **NOTE**:
// This sketch is erroneous, because the nested gemv() call here may involve temporary
// objects which might be destroyed prior to the task getting executed.
// The tasks must involve only raw memory pointers, no proxies.  Therefore we need
// separate async_gpu_matrix and async_matrix classes.

template <typename T, typename U, typename V, typename W>
inline
void
gemv(T alpha, blas::BlasMatrix<T, U, async_tag> const& A,
     blas::BlasVector<T, V, async_tag> const& x, T beta,
     async_vector<W>& y)
{
   DEBUG_CHECK_EQUAL(A.cols(), x.size());
   DEBUG_CHECK_EQUAL(y.size(), A.rows());
   y.wait_for(A.queue());
   y.wait_for(x.queue());
   y.queue().add_task([&]{gemv(alpha, A.base(), x.base(), beta, y.base()); return true;});
   A.wait_for(y.queue());
   x.wait_for(y.queue());
}

// version for an async_gpu_tag
template <typename T, typename U, typename V>
inline
void
gemv(T alpha, blas::BlasMatrix<T, U, async_gpu_tag> const& A,
     blas::BlasVector<T, V, async_gpu_tag> const& x, T beta,
     async_gpu_vector<T>& y)
{
   DEBUG_CHECK_EQUAL(A.cols(), x.size());
   DEBUG_CHECK_EQUAL(y.size(), A.rows());
   y.wait_for(A.queue());
   y.wait_for(x.queue());
   y.queue().add_task(std::bind(cublas::gemv, A.trans(), A.rows(), A.cols(),
                                alpha, A.storage(),
                                A.leading_dimension(), x.storage(), x.stride(),
                                beta, y.storage(), y.stride()));
   A.wait_for(y.queue());
   x.wait_for(y.queue());
}

// another attempt at a generic version; assuming cublas calls have a
// compatible signature to standard blas.
// This is a bit tedious because we need to pass lots of parameters to the lambda function
// by value with automatic variables.

// In async mode we are guaranteed that the objects live longer than this function, because we
// put that intelligence into the destructors of the buffers.  So we don't need reference counting,
// but do need careful handling of temporaries.
template <typename T, typename U, typename V, typename W>
inline
void
gemv(T alpha, blas::BlasMatrix<T, U, async_tag> const& A,
     blas::BlasVector<T, V, async_tag> const& x, T beta,
     async_vector<W>& y)
{
   DEBUG_CHECK_EQUAL(A.cols(), x.size());
   DEBUG_CHECK_EQUAL(y.size(), A.rows());
   y.wait_for(A.queue());
   y.wait_for(x.queue());
   auto Atrans = A.trans();
   auto Arows = A.rows();
   auto Acols = A.cols();
   auto Astorage = A.storage();
   auto Aleading_dimension = A.leading_dimension();
   auto xstorage = x.storage();
   auto xstride = x.stride();
   auto ystorage = y.storage();
   auto ystride = y.stride();
   y.queue().add_task([=](){ gemv(Atrans, Arows, Acols, alpha, Astorage, Aleading_dimension,
                                  xstorage, xstride, beta, ystorage, ystride);});
   A.wait_for(y.queue());
   x.wait_for(y.queue());
}


} // namespace async

#endif
