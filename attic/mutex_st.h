// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/mutex_st.h
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

namespace pthread
{

class mutex
{
   public:
      void lock() {}
      void unlock() {}

      class sentry
      {
         public:
            sentry(mutex& m) {}
      };

      void* get_mutex() { return 0; }

      bool is_owner() const { return true; }
};

typedef mutex debug_mutex;
typedef mutex fast_mutex;

#define CHECK_OWNER(x) /* */

class rec_mutex
{
   public:
      void lock() {}
      void unlock() {}

      class sentry
      {
         public:
            sentry(rec_mutex& m) {}
      };

      void* get_rec_mutex() { return 0; }
};

class condition
{
   public:
      void* get_cond() { return 0; }
};

class semaphore
{
   public:
      explicit semaphore(unsigned int Initial) {}

      void up() {}
      void down() {}
};

class readers_writers
{
   public:
      void write_lock() {}
      void write_unlock() {}

      void read_lock() {}
      void read_unlock() {}

      class write_sentry
      {
         public:
            write_sentry(readers_writers& rw_) {}
         private:
            write_sentry(write_sentry const&);            // not implemented
            write_sentry& operator=(write_sentry const&); // not implemented
      };

      class read_sentry
      {
         public:
            read_sentry(readers_writers& rw_) {}
         private:
            read_sentry(read_sentry const&);            // not implemented
            read_sentry& operator=(read_sentry const&); // not implemented
      };
};

inline
void signal(condition& c)
{
}

inline
void broadcast(condition& c)
{
}

inline
void wait(condition& c, mutex& m)
{
}

} // namespace pthread
