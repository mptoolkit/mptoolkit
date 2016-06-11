// ENDHEADER

#include "parser.h"
#include "color-output.h"

namespace Parser
{

//
// Color support
//

bool Color = is_cout_terminal();

void SetShowColors(ShowColors c)
{
   if (c == ColorNever)
      Color = false;
   else if (c == ColorAuto)
      Color = is_cout_terminal();
   else if (c == ColorAlways)
      Color = true;
   else
      Color = false;
}

std::string ConditionalColorText(std::string s, TerminalColor c)
{
   return Color ? ColorText(s, c) : s;
}

std::string ConditionalColorText(std::string s, TerminalColor c1, TerminalColor c2)
{
   return Color ? ColorText(s, c1, c2) : s;
}

std::string ColorHighlight(std::string s)
{
   return ConditionalColorText(s, TerminalColor::Bold);
}

std::string ColorHighlight(char c)
{
   return ConditionalColorText(std::string(1,c), TerminalColor::Bold);
}

std::string ColorError(std::string s)
{
   return ConditionalColorText(s, TerminalColor::Red, TerminalColor::Bold);
}

std::string ColorWarning(std::string s)
{
   return ConditionalColorText(s, TerminalColor::Magenta, TerminalColor::Bold);
}

std::string ColorNote(std::string s)
{
   return ConditionalColorText(s, TerminalColor::Cyan, TerminalColor::Bold);
}

std::string ColorHint(std::string s)
{
   return ConditionalColorText(s, TerminalColor::Blue, TerminalColor::Bold);
}

std::string ColorPrompt(std::string s)
{
   return ConditionalColorText(s, TerminalColor::Green, TerminalColor::Bold);
}

std::string ColorQuote(std::string s)
{
   return "'" + ColorHighlight(s) + "'";
}

std::string Spaces(int Size)
{
   return std::string(Size, ' ');
}

//
// ParserError
//

ParserError::ParserError(std::string const& Why)
   : CallStack(1, Why), Pos(NULL), End(NULL)
{
   this->AssembleMessage();
}

ParserError::ParserError(ParserError const& Prev, std::string const& Why)
   : CallStack(Prev.CallStack), Pos(Prev.Pos), End(Prev.End), Hint(Prev.Hint)
{
   CallStack.push_back(Why);
   this->AssembleMessage();
}

ParserError::ParserError(std::exception const& Prev, std::string const& Why)
   : CallStack(1, Prev.what()), Pos(NULL), End(NULL)
{
   CallStack.push_back(Why);
   this->AssembleMessage();
}

ParserError::ParserError(std::list<std::string> const& CallStack_, char const* Position,
			 std::string const& Hint_)
   : CallStack(CallStack_), Pos(Position), End(NULL), Hint(Hint_)
{
   this->AssembleMessage();
}

ParserError::ParserError(std::list<std::string> const& CallStack_, char const* Position,
			 char const* End_, std::string const& Hint_)
   : CallStack(CallStack_), Pos(Position), End(End_), Hint(Hint_)
{
   this->AssembleMessage();
}

ParserError::ParserError(std::list<std::string> const& CallStack_, 
			 std::string const& Why, char const* Position, char const* End_,
			 char const* beg, char const* end,
			 std::string const& Hint_)
   : CallStack(CallStack_), Pos(NULL), End(NULL), Hint(Hint_)
{
   std::string Next = Why + "\n" + std::string(beg, end);
   if (Position)
   {
      int Count = End_ ? std::distance(Position, End_) : 1;
      Next = Next + "\n" + std::string(std::distance(beg, Position), ' ')
	 + ColorPrompt(std::string(Count, '^'));
   }
   CallStack.push_back(Next);
   this->AssembleMessage();
}

ParserError 
ParserError::AtPosition(std::string const& Why, char const* Position)
{
   return ParserError(std::list<std::string>(1, Why), Position);
}

ParserError 
ParserError::AtPosition(ParserError const& Prev, char const* Position)
{
   return ParserError(Prev.CallStack, Position, Prev.hint());
}

ParserError 
ParserError::AtPosition(std::exception const& Prev, char const* Position)
{
   return ParserError(std::list<std::string>(1, Prev.what()), Position);
}

ParserError
ParserError::AtRange(std::string const& Why, char const* Start, char const* End)
{
   return ParserError(std::list<std::string>(1, Why), Start, End);
}

ParserError
ParserError::AtRange(ParserError const& Prev, char const* Start, char const* End)
{
   return ParserError(Prev.CallStack, Start, End, Prev.hint());
}

ParserError
ParserError::AtRange(std::exception const& Prev, char const* Start, char const* End)
{
   return ParserError(std::list<std::string>(1, Prev.what()), Start, End);
}

ParserError
ParserError::Finalize(ParserError const& Prev, std::string const& Why, char const* beg, char const* end)
{
   return ParserError(Prev.CallStack, Why, Prev.Pos, Prev.End, beg, end, Prev.hint());
}

ParserError
ParserError::Finalize(std::exception const& Prev, std::string const& Why, char const* beg, char const* end)
{
   return ParserError(std::list<std::string>(1, Prev.what()), Why, NULL, NULL, beg, end);
}

void
ParserError::AssembleMessage()
{
   Msg = ColorError("Parser error: ");
   for (std::list<std::string>::const_iterator I = CallStack.begin(); I != CallStack.end(); ++I)
   {
      if (I == CallStack.begin())
      {
	 Msg = Msg + (*I);
      }
      else
      {
	 Msg = Msg + '\n' + ColorNote("note: ") + (*I);
      }
   }
   if (!Hint.empty())
   {
      Msg += '\n' + ColorHint("hint: ") + Hint;
   }
}

} // namespace Parser
