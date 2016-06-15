// ENDHEADER

#include "parser.h"
#include "common/terminal.h"

namespace Parser
{

//
// Color support
//

bool Color = terminal::is_cout_terminal();

void SetShowColors(ShowColors c)
{
   if (c == ColorNever)
      Color = false;
   else if (c == ColorAuto)
      Color = terminal::is_cout_terminal();
   else if (c == ColorAlways)
      Color = true;
   else
      Color = false;
}

std::string ConditionalColorText(std::string s, terminal::color c)
{
   return Color ? terminal::color_text(s, c) : s;
}

std::string ConditionalColorText(std::string s, terminal::color c1, terminal::color c2)
{
   return Color ? terminal::color_text(s, c1, c2) : s;
}

std::string ColorHighlight(std::string s)
{
   return ConditionalColorText(s, terminal::color::Bold);
}

std::string ColorHighlight(char c)
{
   return ConditionalColorText(std::string(1,c), terminal::color::Bold);
}

std::string ColorError(std::string s)
{
   return ConditionalColorText(s, terminal::color::Red, terminal::color::Bold);
}

std::string ColorWarning(std::string s)
{
   return ConditionalColorText(s, terminal::color::Magenta, terminal::color::Bold);
}

std::string ColorNote(std::string s)
{
   return ConditionalColorText(s, terminal::color::Cyan, terminal::color::Bold);
}

std::string ColorHint(std::string s)
{
   return ConditionalColorText(s, terminal::color::Blue, terminal::color::Bold);
}

std::string ColorPrompt(std::string s)
{
   return ConditionalColorText(s, terminal::color::Green, terminal::color::Bold);
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
   : CallStack(Prev.CallStack), Hint(Prev.Hint), Pos(Prev.Pos), End(Prev.End)
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
   : CallStack(CallStack_), Hint(Hint_), Pos(Position), End(NULL)
{
   this->AssembleMessage();
}

ParserError::ParserError(std::list<std::string> const& CallStack_, char const* Position,
			 char const* End_, std::string const& Hint_)
   : CallStack(CallStack_), Hint(Hint_), Pos(Position), End(End_)
{
   this->AssembleMessage();
}

ParserError::ParserError(std::list<std::string> const& CallStack_, 
			 std::string const& Why, char const* Position, char const* End_,
			 char const* beg, char const* end,
			 std::string const& Hint_)
   : CallStack(CallStack_), Hint(Hint_), Pos(NULL), End(NULL)
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
