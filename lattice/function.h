// -*- C++ -*- $Id$

// Prototypes and argument handling for functions
// Currently the only function arguments we allow are c-numbers
//
// Default values are allowed for arguments.  The default should be a string that is parsed,
// so that it can itself be a function of the other parameters.  Eg,
// x=1/y, y=1/x as two parameters, only one of which needs to be specified.
//
// Parameters may be named.  Unnamed parameters are taken in left-to-right
// order.  Default parameters are parsed in right-to-left order.
//

#if !defined(MPTOOLKIT_LATTICE_FUNCTION_H)
#define MPTOOLKIT_LATTICE_FUNCTION_H

#include "common/trace.h"
#include <pstream/pstream.h>
#include <string>
#include <vector>
#include <ostream>
#include <boost/optional.hpp>
#include "common/formatting.h"

namespace Function
{

// FormalArgument - a Name,Default pair that represents
// a named formal argument of a function

struct FormalArgument
{
   FormalArgument() {}
   FormalArgument(std::string const& arg) : Name(arg) {}
   FormalArgument(char const* arg) : Name(arg) {}
   FormalArgument(std::string const& arg, std::string const& def)
      : Name(arg), Default(def) {}

   std::string Name;
   std::string Default;  // empty for no default
};

std::ostream& operator<<(std::ostream& out, FormalArgument const& arg);

PStream::opstream& operator<<(PStream::opstream& out, FormalArgument const& Arg);
PStream::ipstream& operator>>(PStream::ipstream& in, FormalArgument& Arg);

struct ArgHelper
{
   ArgHelper(std::string Name_) : Name(Name_) {}

   operator FormalArgument() const { return FormalArgument(Name); }

   FormalArgument operator=(std::string const& Default) const
   {
      return FormalArgument(Name, Default);
   }

   FormalArgument operator=(std::complex<double> const& c) const
   {
      return FormalArgument(Name, format_complex(c));
   }

   FormalArgument operator=(double c) const
   {
      return FormalArgument(Name, format_complex(c));
   }

   FormalArgument operator=(int c) const
   {
      return FormalArgument(Name, format_complex(double(c)));
   }

   std::string Name;
};

inline
ArgHelper arg(std::string const& Name)
{
   return ArgHelper(Name);
}

// FormalArgumentList

struct FormalArgumentList : public std::vector<FormalArgument>
{
   FormalArgumentList() {}
};

PStream::opstream& operator<<(PStream::opstream& out, FormalArgumentList const& Args);
PStream::ipstream& operator>>(PStream::ipstream& in, FormalArgumentList& Args);

std::ostream& operator<<(std::ostream& out, FormalArgumentList const& args);

// Parses an argument list of the form
// identifier [= default], [identifier [= default]], ...
FormalArgumentList ParseFormalArguments(std::string const& Args);

// Parameter

struct Parameter
{
   Parameter(std::string const& Name_, std::complex<double> const& Value_)
      : Name(Name_), Value(Value_) {}

   // conversion constructor
   Parameter(std::complex<double> const& v)
      : Value(v) {}

   Parameter(double const& v)
      : Value(v) {}

   std::string Name;
   std::complex<double> Value;
};

struct ParameterHelper
{
   ParameterHelper(std::string const& Name_) : Name(Name_) {}

   Parameter operator=(std::complex<double> const& x) const
   {
      return Parameter(Name, x);
   }

   Parameter operator=(double x) const
   {
      return Parameter(Name, x);
   }

   std::string Name;
};

inline
ParameterHelper par(std::string const& Name)
{
   return ParameterHelper(Name);
}

std::ostream& operator<<(std::ostream& out, Parameter const& param);

typedef std::vector<Parameter> ParameterList;

std::ostream& operator<<(std::ostream& out, ParameterList const& params);

typedef std::pair<std::string, std::complex<double> > Argument;

// a type for holding actual arguments
typedef std::map<std::string, std::complex<double> > ArgumentList;

// function to convert a ParameterList into actual arguments.
// A function must be supplied that parses default arguments,
// which must have the signature 
// complex ParseDefault(ArgumentList, std::string Def)
// which will be called for each parameter that uses the default value.
template <typename Func>
ArgumentList
GetArguments(FormalArgumentList const& Args, ParameterList const& Params,
	     Func ParseDefaultArgument);

class OperatorFunction
{
   public:
      OperatorFunction() {}

      OperatorFunction(FormalArgumentList const& Arguments, std::string const& Definition)
	 : Args(Arguments), Def(Definition) {}

      FormalArgumentList const& Arguments() const { return Args; }

      std::string const& Definition() const { return Def; }

      void set_description(std::string const& s)
      {
	 Desc = s;
      }

      std::string const& description() const { return Desc; }

      // operator() is a helper for defining a function - it allows notation
      // OperatorFunction f;
      // f("arg1", arg("arg2")=x) = y;
      OperatorFunction& operator()(FormalArgument const& arg1);

      OperatorFunction& operator()(FormalArgument const& arg1,
				   FormalArgument const& arg2);

      OperatorFunction& operator()(FormalArgument const& arg1,
				   FormalArgument const& arg2,
				   FormalArgument const& arg3);

      OperatorFunction& operator()(FormalArgument const& arg1,
				   FormalArgument const& arg2,
				   FormalArgument const& arg3,
				   FormalArgument const& arg4);

      OperatorFunction& operator()(FormalArgument const& arg1,
				   FormalArgument const& arg2,
				   FormalArgument const& arg3,
				   FormalArgument const& arg4,
				   FormalArgument const& arg5);

      OperatorFunction& operator()(FormalArgument const& arg1,
				   FormalArgument const& arg2,
				   FormalArgument const& arg3,
				   FormalArgument const& arg4,
				   FormalArgument const& arg5,
				   FormalArgument const& arg6);

      // not the assignment operator, but sets the function definition.
      void operator=(std::string const& Def_) { Def = Def_; }

      OperatorFunction& operator=(OperatorFunction const& f)
      {
	 Args = f.Args;
	 Def = f.Def;
	 return *this;
      }

      FormalArgumentList Args;
      std::string Def;
      std::string Desc;
};

std::ostream& operator<<(std::ostream& out, OperatorFunction const& f);

inline
OperatorFunction&
OperatorFunction::operator()(FormalArgument const& arg1)
{
   Args.clear();
   Args.push_back(arg1);
   return *this;
}

inline
OperatorFunction&
OperatorFunction::operator()(FormalArgument const& arg1,
			     FormalArgument const& arg2)
{
   Args.clear();
   Args.push_back(arg1);
   Args.push_back(arg2);
   return *this;
}

inline
OperatorFunction&
OperatorFunction::operator()(FormalArgument const& arg1,
			     FormalArgument const& arg2,
			     FormalArgument const& arg3)
{
   Args.clear();
   Args.push_back(arg1);
   Args.push_back(arg2);
   Args.push_back(arg3);
   return *this;
}

inline
OperatorFunction&
OperatorFunction::operator()(FormalArgument const& arg1,
			     FormalArgument const& arg2,
			     FormalArgument const& arg3,
			     FormalArgument const& arg4)
{
   Args.clear();
   Args.push_back(arg1);
   Args.push_back(arg2);
   Args.push_back(arg3);
   Args.push_back(arg4);
   return *this;
}

inline
OperatorFunction&
OperatorFunction::operator()(FormalArgument const& arg1,
			     FormalArgument const& arg2,
			     FormalArgument const& arg3,
			     FormalArgument const& arg4,
			     FormalArgument const& arg5)
{
   Args.clear();
   Args.push_back(arg1);
   Args.push_back(arg2);
   Args.push_back(arg3);
   Args.push_back(arg4);
   Args.push_back(arg5);
   return *this;
}

inline
OperatorFunction&
OperatorFunction::operator()(FormalArgument const& arg1,
			     FormalArgument const& arg2,
			     FormalArgument const& arg3,
			     FormalArgument const& arg4,
			     FormalArgument const& arg5,
			     FormalArgument const& arg6)
{
   Args.clear();
   Args.push_back(arg1);
   Args.push_back(arg2);
   Args.push_back(arg3);
   Args.push_back(arg4);
   Args.push_back(arg5);
   Args.push_back(arg6);
   return *this;
}

PStream::opstream& operator<<(PStream::opstream& out, OperatorFunction const& f);
PStream::ipstream& operator>>(PStream::ipstream& in, OperatorFunction& f);

template <typename Func>
ArgumentList
GetArguments(FormalArgumentList const& FormalArgs, ParameterList const& Params,
	     Func ParseDefaultArgument)
{
   Function::ArgumentList Args;

   // keep an index into the current argument number for unnamed arguments
   unsigned CurrentAnonArg = 0;

   for (unsigned i = 0; i < Params.size(); ++i)
   {
      if (Params[i].Name.empty())
      {
	 // we have an anonymous parameter, determine which 
	 // argument it corresponds to
	 while (CurrentAnonArg < FormalArgs.size() 
		&& Args.find(FormalArgs[CurrentAnonArg].Name) != Args.end())
	 {
	    ++CurrentAnonArg;
	 }
	 CHECK(CurrentAnonArg < FormalArgs.size())
	    ("Too many parameters supplied to function")
	    (Params);
	 Args[FormalArgs[CurrentAnonArg++].Name] = Params[i].Value;
      }
      else
      {
	 // named parameter.
	 // Check that it corresponds to an actual argument
	 bool Found = false;
	 for (auto const& x : FormalArgs)
	 {
	    if (x.Name == Params[i].Name)
	    {
	       Found = true;
	       break;
	    }
	 }
	 CHECK(Found)("Named parameter to function does not exist - check spelling!")(Params[i].Name);
	 // and also check that it isn't duplicated
	 if (Args.find(Params[i].Name) != Args.end())
	 {
	    PANIC("Value of named parameter has already been set!")(Params[i].Name);
	 }
	 Args[Params[i].Name] = Params[i].Value;
      }
   }

   // For the remaining arguments, we take the default values.
   // This requires invoking the parser with the existing arguments.
   // We work from right to left, so that default arguments set on the right hand side
   // can be used in expressions further to the left.

   for (int i = int(FormalArgs.size())-1; i >= 0; --i)
   {
      // is the argument not set yet?
      if (Args.find(FormalArgs[i].Name) == Args.end())
      {
	 // make sure that it has a default
	 CHECK(!FormalArgs[i].Default.empty())
	    ("Function argument has no default value")
	    (FormalArgs[i].Name);
	 // Make sure that we don't add the argument until 
	 // we've finished parsing it, so make a temporary before calling Args[...]
	 std::complex<double> c = ParseDefaultArgument(Args, 
						       FormalArgs[i].Default);
	 Args[FormalArgs[i].Name] = c;
      }
      DEBUG_TRACE("value for argument")(FormalArgs[i].Name)(Args[FormalArgs[i].Name]);
   }

   return Args;
}

} // namespace Function

using Function::arg;
using Function::par;

#endif
