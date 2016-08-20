// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/trace.h
//
// Copyright (C) 1997-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

/*
   trace.h

   trace/Precondition checking classes and macros.

   Created: Ian McCulloch   1997

   2004-04-19: Updated to use the preprocessor tricks described in Andrei Alexandrescu and
               John Torjo, C/C++ Users Journal Experts Forum August 2003,
               http://www.cuj.com/documents/s=8464/cujcexp0308alexandr/

               Also updated the logging mechanism so it is customizable via
               set_panic_handler() and set_trace_handler().

               added PANIC() macro.  This replaces THROW()

               Eliminated dependency on trace.cpp by making the global variables static
               members of inline functions.  This is standard C++, the compiler is
               required to just figure it out.
               All functions are either inlined or templates.

               Defining NDEBUG disables DEBUG_*() macros.

               Defining TRACER_DISABLE_ALL_CHECKS disables everything except PANIC()

   2004-05-12: Updated documentation.

   2004-06-28: Added some hacks to (maybe) compile on pre-standard IOStreams

   2004-07-25: Added PRECONDITION_EQUAL, CHECK_EQUAL

   2005-02-11: Added typeid_name<T>() function, to demangle typeid(T).name() on gcc.
               It really has nothing to do with the rest of the header, but there is no
               better place.

   2005-03-02: Modified the string conversion logic to call operator<< in the scope of the
               caller.  No longer requires ADL in namespace ::tracer to find the relevant
               overload.

   2007-02-19: Added set_precision(int p) function to set the precision for floating point output.
               Default set to 16.

   2007-07-18: Added MULTITHREAD support

   2016      : C++-11 support

  BUGS: String messages appear before all variables, rather than in the order they appear.
        (This may be a feature.)

        DUMMY_TRACE might not be optimized away completely on all compilers, need to check this.

  Sample Usage:

  All macros allow additional parameters in brackets, for example given integers i=5, j=1, k=0,
  TRACE(i)(j)(k);

  emits

  TRACE in file tracetest.cpp at line 7:
        i = 5
        j = 1
        k = 0

  Additional parameters can be string literals.  The first such string message appears on
  the first line of the output.  For example,
  TRACE("Hello, World!")(i);

  emits

  TRACE in file tracetest.cpp at line 8: Hello, World!
        i = 5

  The msg() function can be used to print the value of a string variable.  Again, the first
  string message appears on the first line.  For example,
  string S = "Hello";
  string T = "Goodbye";
  TRACE(i).msg(S).msg(T);

  emits

  TRACE in file tracetest.cpp at line 11: Hello
        Goodbye
        i = 5

  The stream insertion operator<< can be used after the final parameter, the combined result
  is treated like a string message.  For example,

  TRACE(i) << "1/3 to 16 figures is " << std::setprecision(16) << (1.0 / 3.0);

  emits

  TRACE in file tracetest.cpp at line 14: 1/3 to 16 figures is 0.3333333333333333
        i = 5

  Unfortunately it is not possible to have a TRACE that consists only of
  a stream insertion operation: TRACE() [with nothing inside the brackets] is not allowed.
  But TRACE("") is OK.

  The available macros are:

  PANIC(...)
     invokes panic() with the specified message

  PRECONDITION(cond) ...
     precondition check: if cond is false, then print a message and invoke panic()

  PRECONDITION_COMPARE(lhs, comp, rhs)
     precondition check: for comp a relational operator,
     if (lhs comp rhs) is false, then print a message and invoke panic()
     shortcut for PRECONDITION(lhs comp rhs)(lhs)(rhs)

  PRECONDITION_EQUAL(lhs, rhs)
     precondition check: if (lhs == rhs) is false, then print a message and invoke panic()
     shortcut for PRECONDITION(lhs == rhs)(lhs)(rhs)
     or PRECONDITION_COMPARE(lhs, ==, rhs)

  CHECK(cond) ...
     assertion check: if cond is false, then print a message and invoke panic()

  CHECK_COMPARE(lhs, comp, rhs)
     precondition check: for comp a relational operator,
     if (lhs comp rhs) is false, then print a message and invoke panic()

  CHECK_EQUAL(lhs, rhs)
     assertion check: if (lhs == rhs) is false, then print a message and invoke panic()
     shortcut for CHECK(lhs == rhs)(lhs)(rhs)

  RANGE_CHECK(variable, lower, upper) ...
     checks that lower <= variable <= upper

  RANGE_CHECK_OPEN(variable, lower, upper) ...
     checks that lower <= variable < upper
     (useful to avoid (upper-1) which has problems if upper is unsigned)

  WARNING(...)
     emits a WARNING: message to the trace stream

  TRACE(...)
     emits a TRACE: message to the trace stream

  TRACE_IF(cond) ...
     If cond evaluates to true, emits a TRACE_IF: message to the trace stream

  DUMMY_TRACE(...)
     expands to nothing - useful for defining conditional checks

  In general it is a good idea to leave precondition checks in the program,
  until it is demonstrated that it affects the performance.
  For cases where the checks are not wanted in the production version,
  rhere are versions that expand to nothing if NDEBUG is defined:

  DEBUG_PRECONDITION
  DEBUG_PRECONDITION_EQUAL
  DEBUG_CHECK
  DEBUG_CHECK_EQUAL
  DEBUG_RANGE_CHECK
  DEBUG_RANGE_CHECK_OPEN
  DEBUG_WARNING
  DEBUG_TRACE
  DEBUG_TRACE_IF

  For living dangerously, if TRACER_DISABLE_ALL_CHECKS is defined then
  ALL of the precondition and trace/warning macros are disabled
  (except PANIC() can never be disabled).
  Use this only as a performance check to see if any checks are slowing
  down the program, or, in rare cases, the precondition checks themselves
  may have unintended side-effects (for example, CHECK(i++ < 1) is guaranteed to
  behave differently if TRACER_DISABLE_ALL_CHECKS is defined :-).

  The assertion checks invoke the tracer::panic(char const* Msg) function, which then invokes
  the 'panic handler'.  The default 'panic handler' writes the message to standard error
  and calls abort().

  Similarly, the trace macros invoke the tracer::trace(char const* Msg) function, which then
  invokes the 'trace handler'.  The default 'trace handler' writes the message to standard error
  and returns.

  It is possible to customize this behaviour via:

  typedef void (*assert_handler)(char const*);

  assert_handler get_panic_handler();
  assert_handler get_trace_handler();

  assert_handler set_panic_handler(assert_handler H);
  assert_handler set_trace_handler(assert_handler H);

  The get_xxx functions return the current handler.
  The set_xxx functions set the handler function to H and return the old handler.

  These functions are in namespace tracer.

  For example,

  void my_panic_handler(char const* Msg)
  {
     system((string("mail -s \"Concerning your software\" billg@microsoft.com <<END\n")
            + Msg + "\nEND\n").c_str());
     invoke_blue_screen();
  }

  int main()
  {
    tracer::set_panic_handler(my_panic_handler);
    // ....
  }

  On systems using GNU libc, a more interesting panic handler is

#include <execinfo.h>
#include <cstdlib>
#include <algorithm>
#include <iterator>

void show_backtrace_handler(char const* Msg)
{
   std::cerr << Msg;

   void* array[20];

   int size = backtrace(array, 20);
   char** strings = backtrace_symbols(array, size);

   std::cerr << "\nStack trace: obtained " << size << " stack frames.\n";
   std::copy(strings, strings+size, std::ostream_iterator<char*>(std::cerr, "\n"));
   std::free(strings);
   std::abort();
}

  To define a conditional check, use the following scheme:

#if defined(MY_LIBRARY_DETAILED_TRACE)
#define TRACE_MY_LIBRARY(x) TRACE(x)
#else
#define TRACE_MY_LIBRARY(x) DUMMY_TRACE(x)
#endif

  And use TRACE_MY_LIBRARY(...) as you would use TRACE.
  If MY_LIBRARY_DETAILED_TRACE is not defined, then it expands to nothing.

*/

#if !defined(MPTOOLKIT_COMMON_TRACE_H)
#define MPTOOLKIT_COMMON_TRACE_H

#if defined(HAVE_CONFIG_H)
#include "config.h"
#if !defined(FUNC_NORETURN)
#error "config.h does not define FUNC_NORETURN"
#endif
#if !defined(__PRETTYFUNC__)
#error "config.h does not define __PRETTYFUNC__"
#endif
#else
#if !defined(FUNC_NORETURN)
#define FUNC_NORETURN
#endif
#if !defined(__PRETTYFUNC__)
#define __PRETTYFUNC__ __func__
#endif
#endif

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <boost/type_traits.hpp>

#ifdef __GNUC__
#include <cxxabi.h>
#endif

#include <boost/lexical_cast.hpp>

#if defined(MULTITHREAD)
#include <pthread.h>
#endif

namespace tracer
{

typedef void (*assert_handler)(char const*);

assert_handler get_panic_handler();
assert_handler get_trace_handler();

assert_handler set_panic_handler(assert_handler H);
assert_handler set_trace_handler(assert_handler H);

// invokes the panic handler
void panic(char const* Msg) FUNC_NORETURN;

// invokes the trace handler
void trace(char const* Msg);

inline
std::string indent_typeid(std::string In)
{
#if defined(INDENT_TYPEID_NAME)
   std::ostringstream Out;
   std::string::const_iterator I = In.begin();
   int Indent = 0;
   while (I != In.end())
   {
      if (*I == '<')
      {
         ++Indent;
         Out << "<\n";
         for (int i = 0; i < Indent*3; ++i) Out << ' ';
      }
      else if (*I == '>')
      {
         --Indent;
         Out << '\n';
         for (int i = 0; i < Indent*3; ++i) Out << ' ';
         Out << ">";
      }
      else if (*I == ',')
      {
         Out << '\n';
         for (int i = 0; i < (Indent-1)*3; ++i) Out << ' ';
         Out << " ,";
      }
      else Out << (*I);
      ++I;
   }
   return Out.str();
#else
   return In;
#endif
}

template<class T>
inline
std::string typeid_name()
{
#ifdef __GNUC__
   int         stat;
   char* Result = abi::__cxa_demangle(typeid(T).name(),0,0,&stat);
   std::string Str = Result ? Result : typeid(T).name();
   std::free(Result);
   return indent_typeid(Str);
#else
   return typeid(T).name();
#endif
}

template<class T>
inline
std::string typeid_name(T const& x)
{
#ifdef __GNUC__
   int         stat;
   char* Result = abi::__cxa_demangle(typeid(x).name(),0,0,&stat);
   std::string Str = Result ? Result : typeid(T).name();
   std::free(Result);
   return indent_typeid(Str);
#else
   return typeid(x).name();
#endif
}

template <class T = void>
struct TracerDummyClass
{
   static int Precision;
};

template <class T>
int TracerDummyClass<T>::Precision = 16;

inline
void set_precision(int p)
{
   TracerDummyClass<>::Precision = p;
}

template <typename T>
std::string ToString(T const& x)
{
  std::ostringstream Buf;
  Buf.precision(TracerDummyClass<>::Precision);
  Buf << x;
  return Buf.str();
}

// weird problem with gcc 3.3.4 and std::boolalpha, we only need it with
// bool anyway
inline
std::string ToBool(bool x)
{
   return x ? "true" : "false";
}

inline
std::string ToString(char const* x)
{
   return x ? x : "null";
}

inline
std::string ToString(char* x)
{
   return x ? x : "null";
}

template <typename T>
inline
T const& fudge_value(T const& x)
{
   return x;
}

inline
char const* fudge_value(char const* x)
{
   return x ? x : "null";
}

inline
char const* fudge_value(char* x)
{
   return x ? x : "null";
}

inline
char const* fudge_value(bool x)
{
   return x ? "true" : "false";
}

inline
std::string fudge_value(float x)
{
   std::ostringstream Buf;
  Buf.precision(12);
  Buf << x;
  return Buf.str();
}

inline
std::string fudge_value(double x)
{
   std::ostringstream Buf;
  Buf.precision(16);
  Buf << x;

  return Buf.str();
}

template <typename T>
T& make_lvalue(T const& x)
{
   return const_cast<T&>(x);
}

inline
std::ostream&& PrepareStream(std::ostream&& out)
{
   out.precision(TracerDummyClass<>::Precision);
   return static_cast<std::ostream&&>(out);
}

#if 1
#define TRACER_CONVERT_TO_STRING(x)                                             \
(static_cast<std::ostringstream&>(::tracer::PrepareStream(std::ostringstream()) << (x)).str())
#else

#define TRACER_CONVERT_TO_STRING(x)                                             \
(static_cast<std::ostringstream&&>(::tracer::PrepareStream(std::ostringstream().flush()) << (x)).str())
#endif

template <int Dummy = 0>
class Assert
{
   public:
      // helpers to enable the recursive macros
      Assert& SMART_ASSERT_A;
      Assert& SMART_ASSERT_B;

      Assert(Assert const& Other);

      ~Assert();

      Assert& msg(std::string const& Msg);

      // Adds a (Name, Value) pair to the list of variables to display.
      // This attempts to disambiguate string literals from other expressions.
      // if Name starts with '"' then it is assumed to be a string literal, in which case
      // msg(Value) is invoked instead.
      // Additionally, if T is an arithmetic type and the string conversion coincides with
      // Name, then nothing is displayed.  This avoids silly messages like "0 = 0".
   //      template <typename T>
      Assert& print_value(char const* Name, //T const& Value,
                          std::string const& ValueStr);

      template <typename T>
      Assert& operator<<(T const& x);

      static Assert MakeAssert(char const* Preamble_, char const* File_,
                               int Line_, char const* Func_,
                               char const* Message_, assert_handler Handler);

  private:
      Assert(char const* Preamble_, char const* File_, int Line_, char const* Func_,
             char const* Message_, assert_handler Handler_)
        : SMART_ASSERT_A(*this), SMART_ASSERT_B(*this),
          ShouldHandle(true), Preamble(Preamble_), File(File_), Line(Line_), Func(Func_),
          Message(Message_),
        Handler(Handler_) {}

      Assert& operator=(Assert const&); // not implemented

      typedef std::pair<std::string, std::string> NameValuePair;
      typedef std::vector<NameValuePair>          VariableListType;

      mutable bool ShouldHandle;
      std::string Preamble, File;
      int Line;
      std::string Func;
      std::string Message;
      assert_handler Handler;
      VariableListType VariableList;
      std::ostringstream ExtraBuf;
};

template <int Dummy>
template <typename T>
Assert<Dummy>& Assert<Dummy>::operator<<(T const& x)
{
   ExtraBuf << ToString(x);
   return *this;
}

template <int Dummy>
inline
Assert<Dummy> Assert<Dummy>::MakeAssert(char const* Preamble_, char const* File_,
                                        int Line_, char const* Func_, char const* Message_,
                                        assert_handler Handler)
{
#if defined(MULTITHREAD)
   std::ostringstream obuf;
   obuf << pthread_self();
   return Assert(("THREAD " + obuf.str() + ": " + Preamble_).c_str(), File_, Line_, Message_, Handler);
#else
   return Assert(Preamble_, File_, Line_, Func_, Message_, Handler);
#endif
}

template <int Dummy>
inline
Assert<Dummy>::Assert(Assert const& Other)
  : SMART_ASSERT_A(*this), SMART_ASSERT_B(*this),
    ShouldHandle(Other.ShouldHandle), Preamble(Other.Preamble),
    File(Other.File), Line(Other.Line),
    Message(Other.Message), Handler(Other.Handler),
     VariableList(Other.VariableList)
{
   ExtraBuf << Other.ExtraBuf.rdbuf();
   Other.ShouldHandle = false;
}

template <int Dummy>
inline
Assert<Dummy>& Assert<Dummy>::msg(std::string const& Msg)
{
   if (!Message.empty() && !Msg.empty())
      Message += '\n' + std::string(Preamble.size()+1, ' ');
   Message += Msg;
   return *this;
}

template <int Dummy>
//template <typename T>
Assert<Dummy>& Assert<Dummy>::print_value(char const* Name,// T const& Value,
                                          std::string const& ValueStr)
{
   if (Name[0] == '"')
     this->msg(ValueStr);
   else
   {
      // do nothing for (literal,literal) pairs
      //      if (boost::is_fundamental<T>::value && ValueStr == Name) return *this;
      if (ValueStr == Name) return *this;

      VariableList.push_back(NameValuePair(Name, ValueStr));
   }
   return *this;
}

template <int Dummy>
Assert<Dummy>::~Assert()
{
   if (!this->ShouldHandle) return;

   this->msg(ExtraBuf.str());
   std::string FullMessage = Preamble + " in file " + File + " in function " + Func
     + " at line " + ToString(Line) + ": " + Message + '\n';
   if (!VariableList.empty())
   {
      int PadSize = int(Preamble.size() + 1);
      if (PadSize < 2) PadSize = 2;
      std::string Padding(PadSize, ' ');
      FullMessage += Padding + VariableList[0].first + " = " + VariableList[0].second + '\n';
      for (int i = 1; i < int(VariableList.size()); ++i)
      {
         FullMessage += Padding + VariableList[i].first + " = " + VariableList[i].second + '\n';
      }
   }
   this->Handler(FullMessage.c_str());
}

inline
void DefaultPanicHandler(char const* msg)
{
   std::cerr << msg;
   std::abort();
}

inline
assert_handler& GetPanicHandler()
{
   static assert_handler Handler = DefaultPanicHandler;
   return Handler;
}

inline
void DefaulttraceHandler(char const* msg)
{
   std::cerr << msg;
}

inline
assert_handler& GettraceHandler()
{
   static assert_handler Handler = DefaulttraceHandler;
   return Handler;
}

inline assert_handler get_panic_handler()
{
   return GetPanicHandler();
}

inline assert_handler get_trace_handler()
{
   return GettraceHandler();
}

inline
assert_handler set_panic_handler(assert_handler H)
{
   assert_handler Old = GetPanicHandler();
   GetPanicHandler() = H;
   return Old;
}

inline
void panic(char const* msg)
{
   GetPanicHandler()(msg);
   abort(); // this is redundant, the panic handler should not return
}

// trace

inline
assert_handler set_trace_handler(assert_handler H)
{
   assert_handler Old = GettraceHandler();
   GettraceHandler() = H;
   return Old;
}

inline
void trace(char const* msg)
{
   GettraceHandler()(msg);
}

struct DummyAssert
{
   DummyAssert& SMART_ASSERT_A;
   DummyAssert& SMART_ASSERT_B;

   DummyAssert() : SMART_ASSERT_A(*this), SMART_ASSERT_B(*this) {}
   DummyAssert(DummyAssert const&) : SMART_ASSERT_A(*this), SMART_ASSERT_B(*this) {}

   static DummyAssert MakeAssert() { return DummyAssert(); }

   template <typename T>
   DummyAssert& msg(T const&) { return *this; }

   //   template <typename T>
   DummyAssert& print_value(char const*, std::string const&) { return *this; }

   template <typename T>
   DummyAssert& operator<<(T const&) { return *this; }
};

template <typename T = int>
struct DummyAssertHandler
{
   static DummyAssert value;
};

template <typename T>
DummyAssert DummyAssertHandler<T>::value;

#define INVOKE_ASSERT(Preamble, Message)                                        \
   ::tracer::Assert<>::MakeAssert((Preamble), __FILE__, __LINE__, __PRETTYFUNC__, (Message), \
                                ::tracer::GetPanicHandler()).SMART_ASSERT_A     \
   /**/

#define INVOKE_PANIC(Preamble, Name, Value)                                     \
   ::tracer::Assert<>::MakeAssert((Preamble), __FILE__, __LINE__, __PRETTYFUNC__, "",           \
                                ::tracer::GetPanicHandler()).                   \
       print_value(Name, TRACER_CONVERT_TO_STRING((Value))).SMART_ASSERT_A      \
   /**/

#define INVOKE_TRACE(Preamble, Name, Value)                                     \
   ::tracer::Assert<>::MakeAssert((Preamble), __FILE__, __LINE__,  __PRETTYFUNC__, "",           \
                                ::tracer::GettraceHandler()).                   \
      print_value(Name, TRACER_CONVERT_TO_STRING((Value))).SMART_ASSERT_A       \
   /**/

#define INVOKE_TRACE_IF(Preamble, Message)                                      \
   ::tracer::Assert<>::MakeAssert((Preamble), __FILE__, __LINE__,  __PRETTYFUNC__, (Message),    \
                                ::tracer::GettraceHandler()).SMART_ASSERT_A     \
   /**/

#if 0

#define INVOKE_DUMMY(junk)                              \
   ::tracer::DummyAssertHandler<>::value.SMART_ASSERT_A \
   /**/

#else

struct DummyAssertType {};

DummyAssertType const DUMMY_ASSERT_A = {};
DummyAssertType const DUMMY_ASSERT_B = {};

template <typename T>
DummyAssertType operator<<(DummyAssertType, T const&)
{
   return DummyAssertType();
}

//int const DUMMY_ASSERT_A = 0;
//int const DUMMY_ASSERT_B = 0;

#define INVOKE_DUMMY(junk)                      \
   ::tracer:: DUMMY_ASSERT_A            \
   /* */

#define DUMMY_ASSERT_A(x) DUMMY_ASSERT_OP(B)
#define DUMMY_ASSERT_B(x) DUMMY_ASSERT_OP(A)

#define DUMMY_ASSERT_OP(next) DUMMY_ASSERT_ ## next     \
   /* */

#endif

#define SMART_ASSERT_A(x) SMART_ASSERT_OP(x, B)
#define SMART_ASSERT_B(x) SMART_ASSERT_OP(x, A)

#define SMART_ASSERT_OP(x, next)                                                                \
   SMART_ASSERT_A.print_value(#x, TRACER_CONVERT_TO_STRING((x))).SMART_ASSERT_ ## next  \
   /**/


// some simple range-check functions.

template <class T, class U, class V>
inline
bool RangeCheck(T const& Var, U const& Min, V const& Max)
{
   return Var < Min || Max < Var;
}

// range check against a half-open range; useful shortcut to avoid problems
// with subtracting 1 from unsigned types; indeed this does not
// require that type V even supports subtraction.
template <class T, class U, class V>
inline
bool RangeCheckOpen(T const& Var, U const& Min, V const& Max)
{
   return Var < Min || !(Var < Max);
}

// dummy version of TRACE so that users can define their own conditional trace functions
#define DUMMY_TRACE(Message)                    \
if (true) ;                                     \
else INVOKE_DUMMY("")

// Debugging versions of the precondition and checking macros.  They expand to nothing
// if NDEBUG is defined.

#if defined(NDEBUG) || defined(TRACER_DISABLE_ALL_CHECKS)

#define DEBUG_CHECK(Condition)                  \
if (true) ;                                     \
else INVOKE_DUMMY("")

#define DEBUG_CHECK_EQUAL(x,y)                  \
if (true) ;                                     \
else INVOKE_DUMMY("")

#define DEBUG_CHECK_COMPARE(x,comp,y)           \
if (true) ;                                     \
else INVOKE_DUMMY("")

#define DEBUG_PRECONDITION(Condition)           \
if (true) ;                                     \
else INVOKE_DUMMY("")

#define DEBUG_PRECONDITION_EQUAL(x,y)           \
if (true) ;                                     \
else INVOKE_DUMMY("")

#define DEBUG_PRECONDITION_COMPARE(x,comp,y)    \
if (true) ;                                     \
else INVOKE_DUMMY("")

#define DEBUG_RANGE_CHECK(Var, Min, Max)        \
if (true) ;                                     \
else INVOKE_DUMMY("")

#define DEBUG_RANGE_CHECK_OPEN(Var, Min, Max)   \
if (true) ;                                     \
else INVOKE_DUMMY("")

#define DEBUG_TRACE(Message)                    \
if (true) ;                                     \
else INVOKE_DUMMY("")

#define DEBUG_TRACE_IF(Message)                 \
if (true) ;                                     \
else INVOKE_DUMMY("")

#define DEBUG_WARNING(Message)                  \
if (true) ;                                     \
else INVOKE_DUMMY("")

#else // !( defined(NDEBUG) || defined(TRACER_DISABLE_ALL_CHECKS) )

#define DEBUG_CHECK(Condition)                          \
if ((Condition)) ;                                      \
else INVOKE_ASSERT("CHECK", #Condition " is false")

#define DEBUG_CHECK_EQUAL(x, y)                                 \
if ((x) == (y)) ;                                               \
else INVOKE_ASSERT("CHECK", #x " == " #y " is false")(x)(y)

#define DEBUG_CHECK_COMPARE(x, comp, y)                                 \
if ((x) comp (y)) ;                                                     \
else INVOKE_ASSERT("CHECK", #x " " #comp " " #y " is false")(x)(y)

#define DEBUG_PRECONDITION(Condition)                           \
if ((Condition)) ;                                              \
else INVOKE_ASSERT("PRECONDITION", #Condition " is false")

#define DEBUG_PRECONDITION_EQUAL(x, y)                                  \
if ((x) == (y)) ;                                                       \
else INVOKE_ASSERT("PRECONDITION", #x " == " #y " is false")(x)(y)

#define DEBUG_PRECONDITION_COMPARE(x, comp, y)                                  \
if ((x) comp (y)) ;                                                             \
else INVOKE_ASSERT("PRECONDITION", #x " " #comp " " #y " is false")(x)(y)

#define DEBUG_RANGE_CHECK(Var, Min, Max)                                        \
if (!::tracer::RangeCheck((Var), (Min), (Max))) ;                               \
else INVOKE_ASSERT("RANGE CHECK",                                               \
                   std::string(#Var " = "                                       \
                               + TRACER_CONVERT_TO_STRING(Var)                  \
                               + " is not in range ["                           \
                               + TRACER_CONVERT_TO_STRING(Min) + ", "           \
                               + TRACER_CONVERT_TO_STRING(Max) + "]").c_str())

#define DEBUG_RANGE_CHECK_OPEN(Var, Min, Max)                                   \
if (!::tracer::RangeCheckOpen((Var), (Min), (Max))) ;                           \
else INVOKE_ASSERT("RANGE CHECK",                                               \
                   std::string(#Var " = "                                       \
                               + TRACER_CONVERT_TO_STRING(Var)                  \
                               + " is not in half-open range ["                 \
                               + TRACER_CONVERT_TO_STRING(Min) + ", "           \
                               + TRACER_CONVERT_TO_STRING(Max) + ")").c_str())

#define DEBUG_TRACE(Msg)                        \
if (false) ;                                    \
else INVOKE_TRACE("TRACE", #Msg, (Msg))

#define DEBUG_TRACE_IF(Cond)                    \
if (!bool(Cond)) ;                              \
else INVOKE_TRACE_IF("TRACE_IF", #Cond " is true")

#define DEBUG_WARNING(Msg)                      \
if (false) ;                                    \
else INVOKE_TRACE("WARNING", #Msg, (Msg))

#endif // ( defined(NDEBUG) || defined(TRACER_DISABLE_ALL_CHECKS) ) else

#if defined(TRACER_DISABLE_ALL_CHECKS)

#define CHECK(Condition)                        \
if (true) ;                                     \
else INVOKE_DUMMY("")

#define CHECK_EQUAL(x,y)                        \
if (true) ;                                     \
else INVOKE_DUMMY("")

#define CHECK_COMPARE(x,comp,y)         \
if (true) ;                                     \
else INVOKE_DUMMY("")

#define PRECONDITION(Condition)                 \
if (true) ;                                     \
else INVOKE_DUMMY("")

#define PRECONDITION_EQUAL(x,y)                 \
if (true) ;                                     \
else INVOKE_DUMMY("")

#define PRECONDITION_COMPARE(x,comp,y)  \
if (true) ;                                     \
else INVOKE_DUMMY("")

#define RANGE_CHECK(Var, Min, Max)              \
if (true) ;                                     \
else INVOKE_DUMMY("")

#define RANGE_CHECK_OPEN(Var, Min, Max)         \
if (true) ;                                     \
else INVOKE_DUMMY("")

#define TRACE(Message)                          \
if (true) ;                                     \
else INVOKE_DUMMY("")

#define TRACE_IF(Message)                       \
if (true) ;                                     \
else INVOKE_DUMMY("")

#define WARNING(Message)                        \
if (true) ;                                     \
else INVOKE_DUMMY("")

#else // !( defined(TRACER_DISABLE_ALL_CHECKS) )

#define CHECK(Condition)                                \
if ((Condition)) ;                                      \
else INVOKE_ASSERT("CHECK", #Condition " is false")

#define CHECK_EQUAL(x, y)                                       \
if ((x) == (y)) ;                                               \
else INVOKE_ASSERT("CHECK", #x " == " #y " is false")(x)(y)

#define CHECK_COMPARE(x, comp, y)                                       \
if ((x) comp (y)) ;                                                     \
else INVOKE_ASSERT("CHECK", #x " " #comp " " #y " is false")(x)(y)

#define PRECONDITION(Condition)                                 \
if ((Condition)) ;                                              \
else INVOKE_ASSERT("PRECONDITION", #Condition " is false")

#define PRECONDITION_EQUAL(x, y)                                        \
if ((x) == (y)) ;                                                       \
else INVOKE_ASSERT("PRECONDITION", #x " == " #y " is false")(x)(y)

#define PRECONDITION_COMPARE(x, comp, y)                                        \
if ((x) comp (y)) ;                                                             \
else INVOKE_ASSERT("PRECONDITION", #x " " #comp " " #y " is false")(x)(y)

#define RANGE_CHECK(Var, Min, Max)                                              \
if (!::tracer::RangeCheck((Var), (Min), (Max))) ;                               \
else INVOKE_ASSERT("RANGE CHECK",                                               \
                   std::string(#Var " = "                                       \
                               + TRACER_CONVERT_TO_STRING(Var)                  \
                               + " is not in range ["                           \
                               + TRACER_CONVERT_TO_STRING(Min) + ", "           \
                               + TRACER_CONVERT_TO_STRING(Max) + "]").c_str())

#define RANGE_CHECK_OPEN(Var, Min, Max)                                         \
if (!::tracer::RangeCheckOpen((Var), (Min), (Max))) ;                           \
else INVOKE_ASSERT("RANGE CHECK",                                               \
                   std::string(#Var " = "                                       \
                               + TRACER_CONVERT_TO_STRING(Var)                  \
                               + " is not in half-open range ["                 \
                               + TRACER_CONVERT_TO_STRING(Min) + ", "           \
                               + TRACER_CONVERT_TO_STRING(Max) + ")").c_str())

#define TRACE(Msg)                              \
if (false) ;                                    \
else INVOKE_TRACE("TRACE", #Msg, (Msg))

#define TRACE_IF(Cond)                          \
if (!bool(Cond)) ;                              \
else INVOKE_TRACE_IF("TRACE_IF", #Cond " is true")

#define WARNING(Msg)                            \
if (false) ;                                    \
else INVOKE_TRACE("WARNING", #Msg, (Msg))

#endif // ( TRACER_DISABLE_ALL_CHECKS ) else

// even if TRACER_DISABLE_ALL_CHECKS is defined, PANIC still functions

#define PANIC(Msg)                              \
if (false) ;                                    \
else INVOKE_PANIC("PANIC", #Msg, (Msg))

} // namespace tracer

#endif
