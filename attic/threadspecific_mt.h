
#if !defined(MULTITHREAD)
#error "Multithread header threadspecific_mt.h cannot be included without MULTITHREAD."
#endif

#include "threads.h"

namespace pthread
{

template <class T>
class thread_specific
{
   public:
      thread_specific();
      explicit thread_specific(T const& InitialValue);
      ~thread_specific();

      operator T&();
      operator T const&() const;

   private:
      thread_specific(thread_specific<T> const&); // not implemented
      thread_specific<T>& operator=(thread_specific<T> const&); // not implemented

      ::pthread_key_t key;
      T Init;
};

} // namespace pthread

#include "threadspecific_mt.cc"

