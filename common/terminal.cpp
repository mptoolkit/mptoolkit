// -*- C++ -*- $Id$

#include "terminal.h"
#include <stdlib.h>

#if HAVE_CONFIG_H
#include "config.h"
#endif

#if HAVE_TERMIOS_H
#include <termios.h>
#endif

#if HAVE_SYS_IOCTL_H
#include <sys/ioctl.h>
#endif

namespace terminal
{

int rows()
{
   // primary option: environment string LINES
   char const* rows_string = getenv("LINES");
   if (rows_string)
   {
      int n_rows = atoi(rows_string);
      if (n_rows > 0)
         return n_rows;
   }
   // fall back to ioctl
#ifdef TIOCGWINSZ
   winsize ws;
   if (!ioctl(1, TIOCGWINSZ, &ws) && ws.ws_row)
      return ws.ws_row;
#endif
   // last resort
   return 25;
}

int columns()
{
   // primary option: environment string COLUMNS
   char const* col_string = getenv("COLUMNS");
   if (col_string)
   {
      int n_cols = atoi(col_string);
      if (n_cols > 0)
         return n_cols;
   }
   // fall back to ioctl
#ifdef TIOCGWINSZ
   winsize ws;
   if (!ioctl(1, TIOCGWINSZ, &ws) && ws.ws_col)
      return ws.ws_col;
#endif
   // last resort
   return 80;
}

std::pair<int, int> tsize()
{
   return std::make_pair(rows(), columns());
}

} // namespace terminal

