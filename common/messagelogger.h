// -*- C++ -*- $Id$

/*
  Message logging mechanism. 

  Created 12/2/99
  Ian McCulloch

  Modified 24/3/99: Added the StreamType template parameter, which defaults to 
  basic_ostream<char>.

  Modified 2000-05-05: changed msg_log and notify_log to use explicit template function
  parameters.  Added the CompileTimeInfo<Level>::Accept nested static bool const
  to give a compile time guess as to whether Accept(Level) would pass.  The compile
  time version only gives false if the run time version would also.
  +DisabledLogger now stores the stream name.
  +panic_log added.

  The Stream of the log objects should always point to something valid, fallback to std::cerr.

  Rewritten 2000-08-09.  Removed the Logger classes.  The interface is now via integers.  If 
  the integer is a compile time constant then the compiler will (hopefully) optimize out
  the I/O, but this infrastructure is now mainly targetted at runtime logging.

  Note that the multiple log functions will get very confused if multiple logs write to the 
  same stream,
  the output will be overlapped.  This could be fixed by separate buffering, but do we care 
  that much?
  A better idea would be to detect duplicate streams and write the thing once only.

  Wish list:
  It would be really cool to have an option to flush the buffer after every output 
  (like the cerr behaviour) for debugging purposes.

*/

#if !defined(ERRORLOGGER_H_DSFHJFDSH3778H897HI37Y37YHUIERUIH)
#define ERRORLOGGER_H_DSFHJFDSH3778H897HI37Y37YHUIERUIH

#include <ostream>
#include <string>
#include <set>

namespace MessageLogger
{

// ShouldFlush returns true if we want to flush the output streams after every write.
// By default, this is true in debug mode, false in non-debug mode.
bool ShouldFlush();

// and we can set it manually.
void ShouldFlush(bool Flush);

// the logger class
 
class LogImpl;

class Logger
{
   public:
      typedef std::ostream ostream_t;

      // the default threshold for run-time logging.
      enum { DefaultRunTimeThreshold = -1 };

      // implicit conversion from string
      Logger(std::string const& s);
      Logger(char const* s);

      Logger(const std::string& LogName_, ostream_t& stream_, 
	     int Threshold_ = DefaultRunTimeThreshold);

      std::string const& Name() const;

      void SetThreshold(int Thresh_) const;
      int GetThreshold() const;

      ostream_t& GetStream() const;
      void SetStream(ostream_t& stream_) const;

      bool Accept(int Level) const;

   private:
      LogImpl* pImpl;
      Logger(); // not implemented
};

// the logging mechanism

class OStreamWrapper
{
   public:
      typedef std::ostream ostream_t;
      typedef std::set<ostream_t*> StreamListType;
      typedef StreamListType::const_iterator const_iterator;

      OStreamWrapper(int N_) : Threshold(N_) {}

      OStreamWrapper(int N_, Logger const& L) : Threshold(N_) { this->push_back(L); }
      OStreamWrapper(int N_, Logger const& L1, Logger const& L2);
      OStreamWrapper(int N_, Logger const& L1, Logger const& L2, Logger const& L3);
      OStreamWrapper(int N_, Logger const& L1, Logger const& L2, Logger const& L3, 
		     Logger const& L4);
      OStreamWrapper(int N_, Logger const& L1, Logger const& L2, Logger const& L3, 
		     Logger const& L4, Logger const& L5);

      void push_back(Logger const& L);

      int GetThreshold() const { return Threshold; }

      const_iterator begin() const { return StreamList.begin(); }
      const_iterator end() const { return StreamList.end(); }

      // I think we only need the ostream_t version, but KAI disagrees.  If we have _only_ the
      // ios_base version, DEC disagrees.

      const OStreamWrapper& operator<<(ostream_t& (*pf)(ostream_t&)) const;
      const OStreamWrapper& operator<<(std::ios_base& (*pf)(std::ios_base&)) const;

   private:
      int Threshold;
      StreamListType StreamList;
};

template <class T>
const OStreamWrapper& 
operator<<(const OStreamWrapper& stream, const T& Thing);

//
// msg_log function is the basic log function.  The template version is for 
// backward compatability,
// it doesn't offer any speed advantages over the non-template.
//

template <int N>
inline
OStreamWrapper msg_log(const Logger& L)
{
   return OStreamWrapper(N, L); 
}

inline
OStreamWrapper msg_log(int N, const Logger& L)
{
   return OStreamWrapper(N, L); 
}

inline
OStreamWrapper msg_log(int N, Logger const& L1, Logger const& L2)
{
   return OStreamWrapper(N, L1, L2); 
}

inline
OStreamWrapper msg_log(int N, Logger const& L1, Logger const& L2, Logger const& L3)
{
   return OStreamWrapper(N, L1, L2, L3); 
}

inline
OStreamWrapper msg_log(int N, Logger const& L1, Logger const& L2, Logger const& L3, 
		       Logger const& L4)
{
   return OStreamWrapper(N, L1, L2, L3, L4); 
}

inline
OStreamWrapper msg_log(int N, Logger const& L1, Logger const& L2, Logger const& L3, 
		       Logger const& L4, Logger const& L5)
{
   return OStreamWrapper(N, L1, L2, L3, L4, L5); 
}

//
// notify_log also writes the log name
//


template <int N>
inline
OStreamWrapper notify_log(Logger const& L)
{
   return OStreamWrapper(N, L) << L.Name() << ": ";
}

inline
OStreamWrapper notify_log(int N, Logger const& L)
{
   return OStreamWrapper(N, L) << L.Name() << ": ";
}

inline
OStreamWrapper notify_log(int N, Logger const& L1, Logger const& L2)
{
   if (L1.Accept(N)) L1.GetStream() << L1.Name() << ": ";
   if (L2.Accept(N)) L2.GetStream() << L2.Name() << ": ";
   return OStreamWrapper(N, L1, L2);
}

inline
OStreamWrapper notify_log(int N, Logger const& L1, Logger const& L2, Logger const& L3)
{
   if (L1.Accept(N)) L1.GetStream() << L1.Name() << ": ";
   if (L2.Accept(N)) L2.GetStream() << L2.Name() << ": ";
   if (L3.Accept(N)) L3.GetStream() << L3.Name() << ": ";
   return OStreamWrapper(N, L1, L2, L3);
}

inline
OStreamWrapper notify_log(int N, Logger const& L1, Logger const& L2, Logger const& L3, 
			  Logger const& L4)
{
   if (L1.Accept(N)) L1.GetStream() << L1.Name() << ": ";
   if (L2.Accept(N)) L2.GetStream() << L2.Name() << ": ";
   if (L3.Accept(N)) L3.GetStream() << L3.Name() << ": ";
   if (L4.Accept(N)) L4.GetStream() << L4.Name() << ": ";
   return OStreamWrapper(N, L1, L2, L3, L4);
}

inline
OStreamWrapper notify_log(int N, Logger const& L1, Logger const& L2, Logger const& L3, 
			  Logger const& L4, Logger const& L5)
{
   if (L1.Accept(N)) L1.GetStream() << L1.Name() << ": ";
   if (L2.Accept(N)) L2.GetStream() << L2.Name() << ": ";
   if (L3.Accept(N)) L3.GetStream() << L3.Name() << ": ";
   if (L4.Accept(N)) L4.GetStream() << L4.Name() << ": ";
   if (L5.Accept(N)) L5.GetStream() << L5.Name() << ": ";
   return OStreamWrapper(N, L1, L2, L3, L4, L5);
}

// 
// panic log doesn't check the threshold
//

inline
std::ostream& panic_log(Logger const& L)
{
   return L.GetStream() << L.Name() << " PANIC! "; 
}

// other inlines

inline
void OStreamWrapper::push_back(Logger const& L)
{
   if (L.Accept(Threshold)) StreamList.insert(&L.GetStream());
}

} // namespace MessageLogger

// export the log functions to the global namespace, since DEC doesn't seem to have 
// Koenig lookup ???

using MessageLogger::notify_log;
using MessageLogger::msg_log;
using MessageLogger::panic_log;

#include "messagelogger.cc"

#endif

