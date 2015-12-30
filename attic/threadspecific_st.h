/*
  single-threaded (dummy) version of threadspecific.h

  Created 2002-10-26 Ian McCulloch
*/

#if defined(MULTITHREAD)
#error "Cannot use single thread header threadspecific_st.h with MULTITHREAD."
#endif

namespace pthread
{

template <class T>
class thread_specific
{
   public:
      thread_specific() {}
      thread_specific(T const& Init) : data(Init) {}

      operator T&() { return data; }
      operator T const&() const { return data; }

   private:
      thread_specific(thread_specific<T> const&); // not implemented
      thread_specific<T>& operator=(thread_specific<T> const&); // not implemented

      T data;
};

} // namespace pthread

