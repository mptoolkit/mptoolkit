//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// cuda/cuda.ipp
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

namespace cuda
{

//
// stream
//

inline
stream::stream() : stream_(stream::Allocate())
{
}

inline
stream::stream(stream&& other) : stream_(other.stream_)
{
   other.stream_ = nullptr;
}

inline
stream&
stream::operator=(stream&& other)
{
   if (stream_)
   {
      std::lock_guard<std::mutex> lock(stream::FreeListMutex);
      stream::FreeList.push_back(stream_);
   }
   stream_ = other.stream_;
   other.stream_ = nullptr;
   return *this;
}

inline
stream::~stream()
{
   if (stream_)
   {
      std::lock_guard<std::mutex> lock(stream::FreeListMutex);
      stream::FreeList.push_back(stream_);
      stream_ = nullptr;
   }
}

inline
cudaStream_t
stream::Allocate()
{
   std::lock_guard<std::mutex> lock(stream::FreeListMutex);
   if (FreeList.empty())
   {
      for (int i = 0; i < 10; ++i)
      {
	 cudaStream_t s;
	 check_error(cudaStreamCreate(&s));
	 FreeList.push_back(s);
      }
   }
   cudaStream_t s = FreeList.top();
   FreeList.pop();
   return s;
}

//
// event
//

inline
event::event() : event_(event::Allocate())
{
}

inline
event::event(event&& other) : event_(other.event_)
{
   other.event_ = nullptr;
}

inline
event&
event::operator=(event&& other)
{
   if (event_)
   {
      std::lock_guard<std::mutex> lock(event::FreeListMutex);
      event::FreeList.push_back(event_);
   }
   event_ = other.event_;
   other.event_ = nullptr;
   return *this;
}

inline
event::~event()
{
   if (event_)
   {
      std::lock_guard<std::mutex> lock(event::FreeListMutex);
      event::FreeList.push_back(event_);
      event_ = nullptr;
   }
}

inline
cudaEvent_t
event::Allocate()
{
   std::lock_guard<std::mutex> lock(event::FreeListMutex);
   if (FreeList.empty())
   {
      for (int i = 0; i < 10; ++i)
      {
	 cudaEvent_t s;
	 // for debugging add cudaEventBlockingSync to the flags here
	 check_error(cudaEventCreateWithFlags(&s, cudaEventDisableTiming));
	 FreeList.push_back(s);
      }
   }
   cudaEvent_t s = FreeList.top();
   FreeList.pop();
   return s;
}


bool
event::is_complete() const
{
   cudaError_t e = cudaEventQuery(event_);
   if (e == cudaSuccess)
      return true;
   if (e == cudaErrorNotReady)
      return false;
   throw error(e);
}

//
// timer
//

inline
timer::timer() : start_(Allocate()), stop_(Allocate())
{
}

inline
timer::timer() : timer_(timer::Allocate())
{
}

inline
timer::timer(timer&& other) : start_(other.start_), stop_(other.stop_)
{
   other.start_ = nullptr;
   other.stop_ = nullptr;
}

inline
timer&
timer::operator=(timer&& other)
{
   if (start_ || stop)
   {
      std::lock_guard<std::mutex> lock(timer::FreeListMutex);
      if (start)
	 timer::FreeList.push_back(start_);
      if (stop_)
	 timer::FreeList.push_back(stop_);
   }
   start_ = other.start_;
   other.start_ = nullptr;
   stop_ = other.stop_;
   other.stop_ = nullptr;
   return *this;
}

inline
timer::~timer()
{
   if (start_ || stop_)
   {
      std::lock_guard<std::mutex> lock(timer::FreeListMutex);
      if (start_)
      {
	 timer::FreeList.push_back(start_);
	 start_ = nullptr;
      }
      if (stop_)
      {
	 timer::FreeList.push_back(stop_);
	 stop_ = nullptr;
      }
   }
}

inline
cudaEvent_t
timer::Allocate()
{
   std::lock_guard<std::mutex> lock(timer::FreeListMutex);
   if (FreeList.empty())
   {
      for (int i = 0; i < 10; ++i)
      {
	 cudaEvent_t s;
	 check_error(cudaEventCreateWithFlags(&s, cudaEventDefault));
	 FreeList.push_back(s);
      }
   }
   cudaEvent_t s = FreeList.top();
   FreeList.pop();
   return s;
}

inline
float
timer::elapsed_time_ms()
{
   float Result;
   cudaError_t e = cudaEventElapsedTime(&Result, start_, stop_);
   if (e == cudaSuccess)
      return Result;
   throw error(e);
}

bool
timer::is_started() const
{
   cudaError_t e = cudaEventQuery(start_);
   if (e == cudaSuccess)
      return true;
   if (e == cudaErrorNotReady)
      return false;
   throw error(e);
}

bool
timer::is_complete() const
{
   if (!this->is_started())
      return false;
   cudaError_t e = cudaEventQuery(stop_);
   if (e == cudaSuccess)
      return true;
   if (e == cudaErrorNotReady)
      return false;
   throw error(e);
}

} // namespace cuda
