// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// attic/threadspecific_mt.cc
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

extern "C" void thread_specific_dtor(void* tsd);

class tsd_base
{
   public:
      virtual ~tsd_base() = 0;
};

template <class T>
struct tsd_item : public tsd_base
{
   public:
      tsd_item() { item = new T(); }
      ~tsd_item() { delete item; }

      T* get() { return item; }

   private:
      T* item;
};

template <class T>
thread_specific<T>::thread_specific()
{
   ::pthread_key_create(&key, &thread_specific_dtor);
}

template <class T>
thread_specific<T>::thread_specific(T const& InitialValue)
  : Init(InitialValue)
{
   ::pthread_key_create(&key, &thread_specific_dtor);
}

template <class T>
thread_specific<T>::~thread_specific()
{
   ::pthread_key_delete(key);
}

template <class T>
T&
thread_specific<T>::operator T&()
{
   void* Val = ::pthread_getspecific(key);
   if (Val == NULL)
   {
      Val = static_cast<tsd_base*>(new tsd_item<T>(Init));
      ::pthread_setspecific(key, Val);
   }
   return *static_cast<tsd_item<T>*>(static_cast<tsd_base*>(Val))->get();
}

template <class T>
T const&
thread_specific<T>::operator T const&() const
{
   void* Val = ::pthread_getspecific(key);
   if (Val == NULL)
   {
      Val = static_cast<tsd_base*>(new tsd_item<T>(Init));
      ::pthread_setspecific(key, Val);
   }
   return *static_cast<tsd_item<T>*>(static_cast<tsd_base*>(Val))->get();
}

} // namespace pthread
