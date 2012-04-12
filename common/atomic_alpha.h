// -*- C++ -*- $Id$

/*
  Tru64/alpha version of atomic.h

  Created 2002-08-12 Ian McCulloch
*/

#include <c_asm.h>
#include <sys/machine/builtins.h>

// specializations of atomic for int and long.
// 

template <>
class atomic<int>
{
   public:
      atomic() {}
      atomic(int t) : value(t) {}

      int read() const { return value; }
      void operator=(int t) { value = t; }

      int operator++() { return __ATOMIC_INCREMENT_LONG(&value)+1; }
      int operator++(int) { return __ATOMIC_INCREMENT_LONG(&value); }

      int operator--() { return __ATOMIC_DECREMENT_LONG(&value)-1; }
      int operator--(int) { return __ATOMIC_DECREMENT_LONG(&value); }

      int operator+=(int t) { return __ATOMIC_ADD_LONG(&value, t); }
      int operator-=(int t) { return __ATOMIC_ADD_LONG(&value, -t); }

      int exchange(int t) { return __ATOMIC_EXCH_LONG(&value, t); }

   private:
      atomic(atomic const&);  // not implemented
      atomic& operator=(atomic const&); // not implemented

      int volatile value;
};

template <>
class atomic<unsigned int>
{
   public:
      atomic() {}
      atomic(unsigned int t) : value(t) {}

      unsigned int read() const { return value; }
      void operator=(unsigned int t) { value = t; }

      unsigned int operator++() { return static_cast<unsigned int>(__ATOMIC_INCREMENT_LONG(&value)+1); }
      unsigned int operator++(int) { return static_cast<unsigned int>(__ATOMIC_INCREMENT_LONG(&value)); }

      unsigned int operator--() { return static_cast<unsigned int>(__ATOMIC_DECREMENT_LONG(&value)-1); }
      unsigned int operator--(int) { return static_cast<unsigned int>(__ATOMIC_DECREMENT_LONG(&value)); }

      unsigned int operator+=(unsigned int t) { return static_cast<unsigned int>(__ATOMIC_ADD_LONG(&value, t)); }
      unsigned int operator-=(unsigned int t) { return static_cast<unsigned int>(__ATOMIC_ADD_LONG(&value, -t)); }

      unsigned int exchange(unsigned int t) { return static_cast<unsigned int>(__ATOMIC_EXCH_LONG(&value, t)); }

   private:
      atomic(atomic const&);  // not implemented
      atomic& operator=(atomic const&); // not implemented

      unsigned int volatile value;
};

template <>
class atomic<long>
{
   public:
      atomic() {}
      atomic(long t) : value(t) {}

      long read() const { return value; }
      void operator=(long t) { value = t; }

      long operator++() { return __ATOMIC_INCREMENT_QUAD(&value)+1; }
      long operator++(int) { return __ATOMIC_INCREMENT_QUAD(&value); }

      long operator--() { return __ATOMIC_DECREMENT_QUAD(&value)-1; }
      long operator--(int) { return __ATOMIC_DECREMENT_QUAD(&value); }

      long operator+=(long t) { return __ATOMIC_ADD_QUAD(&value, t); }
      long operator-=(long t) { return __ATOMIC_ADD_QUAD(&value, -t); }

      long exchange(long t) { return __ATOMIC_EXCH_QUAD(&value, t); }

   private:
      atomic(atomic const&);  // not implemented
      atomic& operator=(atomic const&); // not implemented

      long volatile value;
};

template <>
class atomic<unsigned long>
{
   public:
      atomic() {}
      atomic(unsigned long t) : value(t) {}

      unsigned long read() const { return value; }
      void operator=(unsigned long t) { value = t; }

      unsigned long operator++() { return static_cast<unsigned long>(__ATOMIC_INCREMENT_QUAD(&value)+1); }
      unsigned long operator++(int) { return static_cast<unsigned long>(__ATOMIC_INCREMENT_QUAD(&value)); }

      unsigned long operator--() { return static_cast<unsigned long>(__ATOMIC_DECREMENT_QUAD(&value)-1); }
      unsigned long operator--(int) { return static_cast<unsigned long>(__ATOMIC_DECREMENT_QUAD(&value)); }

      unsigned long operator+=(unsigned long t) { return static_cast<unsigned long>(__ATOMIC_ADD_QUAD(&value, t)); }
      unsigned long operator-=(unsigned long t) { return static_cast<unsigned long>(__ATOMIC_ADD_QUAD(&value, -t)); }

      unsigned long exchange(unsigned long t) { return static_cast<unsigned long>(__ATOMIC_EXCH_QUAD(&value, t)); }

   private:
      atomic(atomic const&);  // not implemented
      atomic& operator=(atomic const&); // not implemented

      unsigned long volatile value;
};

class AtomicRefCount
{
   public:
      AtomicRefCount();
      AtomicRefCount(int Initial);

      // atomically increments the counter
      void operator++();

      // atomically decrements the counter.  Acts as a memory barrier.
      int operator--();

      // returns true if the counter is zero.  Should this act as a memory barrier?
      bool is_zero() const;

      int value() const { return Count; }  // reads are atomic on alpha

   private:
      AtomicRefCount(AtomicRefCount const&); // not implemented
      AtomicRefCount& operator=(AtomicRefCount const&); // not implemented

      int volatile Count;
};

inline
AtomicRefCount::AtomicRefCount() : Count(0)
{
}

inline
AtomicRefCount::AtomicRefCount(int Initial) : Count(Initial)
{
}

inline
void AtomicRefCount::operator++()
{
   __ATOMIC_INCREMENT_LONG(&Count);
}

inline
int AtomicRefCount::operator--()
{
   asm("wmb");
   if (__ATOMIC_DECREMENT_LONG(&Count) == 1)
   {
      asm("mb");
      return 0;
   }
   return 1;
}

inline
bool AtomicRefCount::is_zero() const
{
   // this is rather conservative, we can probably get away with a simple "return Count == 0".
   // This version forces a concurrent ATOMIC_INC/DECREMENT to retry.
   return !__ATOMIC_ADD_LONG(const_cast<int volatile*>(&Count), 0);
}

inline void memory_barrier()
{
   asm("mb");
}

inline void write_memory_barrier()
{
   asm("wmb");
}

inline void read_memory_barrier()
{
   asm("mb");
}
