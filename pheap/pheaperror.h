// -*- C++ -*-

// Error classes for pheap

#if !defined(MPTOOLKIT_PHEAP_PHEAPERROR_H)

#include <exception>
#include <string>

class PHeapError : public std::runtime_error
{
   public:
      PHeapError(std::string const& s) : std::runtime_error(s) {}
};

class PHeapVersionMismatch : public PHeapError
{
   public:
      PHeapVersionMismatch(std::string const& s) : PHeapError(s) {}
};

#endif
