// -*- C++ -*- $Id$

#include <map>
#include <iostream>

namespace MessageLogger
{

namespace Private
{

typedef std::map<std::string, LogImpl*> LogMapType;

struct Cleanup
{
   ~Cleanup();
};

template <typename T = void>
struct Variables
{
   static std::map<std::string, LogImpl*> LogMap;
   static bool ShouldFlush;
   static Cleanup Cleaner;
};

template <typename T>
std::map<std::string, LogImpl*> Variables<T>::LogMap;

template <typename T>
bool Variables<T>::ShouldFlush;

template <typename T>
Cleanup Variables<T>::Cleaner;

inline
LogImpl* Lookup(std::string const& Name)
{
   LogMapType::const_iterator I = Variables<>::LogMap.find(Name);
   return (I != Variables<>::LogMap.end()) ? I->second : NULL;
}

inline
void Register(LogImpl* Log);

} // namespace Private

class LogImpl
{
   public:
      typedef Logger::ostream_t ostream_t;

      LogImpl(std::string const& LogName_) : LogName(LogName_), stream(&std::cerr), Threshold(-1)
      {
         Private::Register(this);
      }

      LogImpl(const std::string& LogName_, ostream_t& stream_, 
              int Threshold_ = Logger::DefaultRunTimeThreshold)
         : LogName(LogName_), stream(&stream_), Threshold(Threshold_) 
      {
         Private::Register(this);
      }

      std::string const& Name() const { return LogName; }

      void SetThreshold(int Thresh_) { Threshold = Thresh_; }
      int GetThreshold() const { return Threshold; }

      ostream_t& GetStream() const { return *stream; }
      void SetStream(ostream_t& stream_) { stream = &stream_; }

      bool Accept(int Level) const { return this->GetThreshold() >= Level; }

   private:
      std::string LogName;
      ostream_t* stream;
      int Threshold;
};


namespace Private
{

inline
void Register(LogImpl* Log)
{
   Variables<>::LogMap[Log->Name()] = Log;
}

inline
Cleanup::~Cleanup()
{
   for (LogMapType::const_iterator I = Variables<>::LogMap.begin();
        I != Variables<>::LogMap.end(); ++I)
   {
      I->second->GetStream().flush();
   }
}

}

inline
Logger::Logger(std::string const& s)
{
   pImpl = Private::Lookup(s);
   if (!pImpl)
      pImpl = new LogImpl(s);
}

inline
Logger::Logger(char const* s)
{
   pImpl = Private::Lookup(s);
   if (!pImpl)
      pImpl = new LogImpl(s);
}

inline
Logger:: Logger(const std::string& LogName_, ostream_t& stream_, int Threshold_)
   : pImpl(new LogImpl(LogName_, stream_, Threshold_))
{
}

inline
std::string const& Logger::Name() const
{
   return pImpl->Name();
}

inline
void Logger::SetThreshold(int Thresh_) const
{ 
   pImpl->SetThreshold(Thresh_);
}

inline
int Logger::GetThreshold() const
{
   return pImpl->GetThreshold();
}

inline
Logger::ostream_t& Logger::GetStream() const
{
   return pImpl->GetStream();
}

inline
void Logger::SetStream(ostream_t& stream_) const
{
   pImpl->SetStream(stream_);
}

inline
bool Logger::Accept(int Level) const
{
   return pImpl->Accept(Level);
}

inline
bool ShouldFlush() 
{ 
   return Private::Variables<>::ShouldFlush; 
}

inline
void ShouldFlush(bool Flush)
{
   Private::Variables<>::ShouldFlush = Flush;
}

inline
const OStreamWrapper&
OStreamWrapper::operator<<(ostream_t& (*pf)(ostream_t&)) const
{
   for (const_iterator I = begin(); I != end(); ++I)
   {
      (**I) << pf;
      if (ShouldFlush()) (**I).flush();
   }
   return *this;
}

inline
const OStreamWrapper&
OStreamWrapper::operator<<(std::ios_base& (*pf)(std::ios_base&)) const
{
   for (const_iterator I = begin(); I != end(); ++I)
   {
      (**I) << pf;
      if (ShouldFlush()) (**I).flush();
   }
   return *this;
}

inline
OStreamWrapper::OStreamWrapper(int N_, Logger const& L1, Logger const& L2) 
  : Threshold(N_) 
{ 
   this->push_back(L1); 
   this->push_back(L2); 
}

inline
OStreamWrapper::OStreamWrapper(int N_, Logger const& L1, Logger const& L2, Logger const& L3) 
  : Threshold(N_) 
{ 
   this->push_back(L1); 
   this->push_back(L2); 
   this->push_back(L3);
}

inline
OStreamWrapper::OStreamWrapper(int N_, Logger const& L1, Logger const& L2, Logger const& L3, 
                               Logger const& L4) 
  : Threshold(N_) 
{ 
   this->push_back(L1); 
   this->push_back(L2); 
   this->push_back(L3);
   this->push_back(L4);
}

inline
OStreamWrapper::OStreamWrapper(int N_, Logger const& L1, Logger const& L2, Logger const& L3, 
                               Logger const& L4, Logger const& L5) 
  : Threshold(N_) 
{ 
   this->push_back(L1); 
   this->push_back(L2); 
   this->push_back(L3);
   this->push_back(L4);
   this->push_back(L5);
}

template <class T>
const OStreamWrapper& 
operator<<(const OStreamWrapper& stream, const T& Thing)
{
   for (OStreamWrapper::const_iterator I = stream.begin(); I != stream.end(); ++I)
   {
      (**I) << Thing;
      if (ShouldFlush()) (**I).flush();
   }
   return stream;
}

} // namespace MessageLogger
