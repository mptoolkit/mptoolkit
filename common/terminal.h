// -*- C++ -*- $Id$

/*
   terminal.h

   Some simple terminal control functions.

   Created 2006-05-19 Ian McCulloch
*/

#if !defined(TERMINAL_H_DSJFHJFH43Y439Y437YYH9YREW)
#define TERMINAL_H_DSJFHJFH43Y439Y437YYH9YREW

#include <utility>

namespace terminal
{

// return the size of the output terminal, as (rows, cols)
std::pair<int, int> size();

// return the number of rows of the output terminal
int rows();

// return the number of columns in the output terminal
int columns();

} // namespace terminal


#endif
