// -*- C++ -*- $Id$

std::string indent_typeid(std::string In)
{
   std::ostringstream Out;
   std::string::const_iterator I = In.begin();
   int Indent = 0;
   while (I != In.end())
   {
      if (*I == '<')
      {
	 ++Indent;
	 Out << '\n';
	 for (int i = 0; i < Indent*3; ++i) Out << ' ';
	 Out << '<';
      }
      else if (*I == '>')
      {
	 --Indent;
	 Out << '\n';
	 for (int i = 0; i < Indent*3; ++i) Out << ' ';
	 out << ">";
      }
      else if (*I == ',')
      {
	 Out << '\n';
	 for (int i = 0; i < (Indent-1)*3; ++i) Out << ' ';
	 Out << " , ";
      }
   }
   return Out.str();
}
