// -*- C++ -*- $Id$

#include "function.h"
#include "parser/parser.h"
#include <boost/spirit/include/classic_lists.hpp>

namespace Function
{

using namespace Parser;

std::ostream& operator<<(std::ostream& out, FormalArgument const& arg)
{
   out << arg.Name;
   if (!arg.Default.empty())
      out << '=' << (arg.Default);
   return out;
}

PStream::opstream& operator<<(PStream::opstream& out, FormalArgument const& Arg)
{
   out << Arg.Name;
   out << Arg.Default;
   return out;
}

PStream::ipstream& operator>>(PStream::ipstream& in, FormalArgument& Arg)
{
   in >> Arg.Name;
   in >> Arg.Default;
   return in;
}

std::ostream& operator<<(std::ostream& out, FormalArgumentList const& args)
{
   out << '(';
   if (!args.empty())
   {
      out << args[0];
      for (unsigned i = 1; i < args.size(); ++i)
      {
	 out << ", " << args[i];
      }
   }
   out << ')';
   return out;
}

PStream::opstream& operator<<(PStream::opstream& out, FormalArgumentList const& Arg)
{
   out << static_cast<std::vector<FormalArgument> const&>(Arg);
   return out;
}

PStream::ipstream& operator>>(PStream::ipstream& in, FormalArgumentList& Arg)
{
   in >> static_cast<std::vector<FormalArgument>&>(Arg);
   return in;
}

std::ostream& operator<<(std::ostream& out, Parameter const& param)
{
   if (!param.Name.empty())
      out << (param.Name) << '=';
   out << param.Value;
   return out;
}

std::ostream& operator<<(std::ostream& out, ParameterList const& params)
{
   out << '(';
   if (!params.empty())
   {
      out << params[0];
      for (unsigned i = 1; i < params.size(); ++i)
      {
	 out << ", " << params[i];
      }
   }
   out << ')';
   return out;
}

PStream::opstream& operator<<(PStream::opstream& out, OperatorFunction const& f)
{
   out << f.Args;
   out << f.Def;
   return out;
}

PStream::ipstream& operator>>(PStream::ipstream& in, OperatorFunction& f)
{
   in >> f.Args;
   in >> f.Def;
   return in;
}

#if 0
std::map<std::string, std::complex<double> >
GetArguments(ParameterList const& Params, FormalArgumentList const& FormalArgs)
{
   CHECK(Params.size() <= FormalArgs.size())("Too many arguments to function!")
      (FormalArgs)(Params);

   std::vector<ActualArgument> Result;

   // Add each parameter to the arguments
   for (unsigned i = 0; i < Params.size(); ++i)
   {
      if (Params[i].Name)
      {
	 // make sure the Name actually corresponds to an argument
	 int n = 0;
	 while (n < FormalArgs.size() && *Params[i].Name != FormalArgs[n].Name)
	    ++n;
	 CHECK(n != FormalArgs.size())("Unknown parameter name!")
	    (*Params[p].Name)(FormalArgs);

	 Result[*Params[i].Name] = Params[i].Value;
      }
      else
      {
	 // search for the first unset argument
	 FormalArgumentList::const_iterator I = FormalArgs.begin();
	 while (I != FormalArgs.end() && Result.find(I->Name) != Result.end())
	    ++I;
	 CHECK(I != FormalArgs.end())("Too many parameters "
				      "supplied to function!")
	    (FormalArgs)(Params);
	 Result[I->Name] = Params[i].Value;
      }
   }

   // set any remaining parameters to default, or error if there is
   // no default
   for (unsigned i = 0; i < FormalArgs.size(); ++i)
   {
      if (Result.find(FormalArgs[i].Name == Result.end()))
      {
	 CHECK(FormalArgs[i].Value)("Parameters with no default must be set")
	    (FormalArgs[i].Name);

	 Result[FormalArgs[i].Name] = *FormalArgs[i].Value;
      }
   }
   return Result;
}
#endif

struct add_arg
{
   add_arg(std::stack<std::string>& IdentStack_, FormalArgumentList& ArgList_) :
      IdentStack(IdentStack_), ArgList(ArgList_) {}

   void operator()(char const*, char const*) const
   {
      CHECK(!IdentStack.empty());
      ArgList.push_back(FormalArgument(IdentStack.top()));
      IdentStack.pop();
      CHECK(IdentStack.empty());
   }

   std::stack<std::string>& IdentStack;
   FormalArgumentList& ArgList;
};

struct add_arg_default
{
   add_arg_default(std::stack<std::string>& IdentStack_, FormalArgumentList& ArgList_) :
      IdentStack(IdentStack_), ArgList(ArgList_) {}

   void operator()(char const* Start, char const* End) const
   {
      CHECK(!IdentStack.empty());
      ArgList.push_back(FormalArgument(IdentStack.top(), std::string(Start, End)));
      IdentStack.pop();
      CHECK(IdentStack.empty());
   }

   std::stack<std::string>& IdentStack;
   FormalArgumentList& ArgList;
};

struct ArgumentListParser : public grammar<ArgumentListParser>
{
   typedef std::stack<std::string>    IdentStackType;

   ArgumentListParser(IdentStackType& identifier_stack_,
		      FormalArgumentList& arg_list_)
      : identifier_stack(identifier_stack_), arg_list(arg_list_)
   {}
   
   template <typename ScannerT>
   struct definition
   {
      definition(ArgumentListParser const& self)
      {
	 identifier = lexeme_d[alpha_p >> *(alnum_p | '_')]
	    [Parser::push_identifier(self.identifier_stack)];

	 expression_string = lexeme_d[+((anychar_p - chset<>("(){}[]"))
					| (ch_p('(') >> expression_string >> ch_p(')'))
					| (ch_p('[') >> expression_string >> ch_p(']'))
					| (ch_p('{') >> expression_string >> ch_p('}'))
					)];

	 argument = identifier >> 
	    (('=' >> expression_string[add_arg_default(self.identifier_stack, self.arg_list)])
	     | eps_p[add_arg(self.identifier_stack, self.arg_list)]);

	 argument_list = !list_p(argument, ',');
      }
      
      rule<ScannerT> identifier, expression_string, argument, argument_list;
      rule<ScannerT> const&
      start() const { return argument_list; }
   };
   
   IdentStackType& identifier_stack;
   FormalArgumentList& arg_list;
};

FormalArgumentList ParseFormalArguments(std::string const& Args)
{
   ArgumentListParser::IdentStackType IdentStack;
   FormalArgumentList Result;

   ArgumentListParser Parser(IdentStack, Result);

   parse_info<> info = parse(Args.c_str(), Parser, space_p);
   if (!info.full)
   {
      PANIC("Argument list parser failed, stopped at")(info.stop);
   }

   CHECK(IdentStack.empty());

   return Result;
}

std::ostream& operator<<(std::ostream& out, OperatorFunction const& f)
{
   out << f.Args << " = " << f.Def;
   return out;
}

} // namespace Function

