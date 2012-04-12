// -*- C++ -*- $Id$

#include "quantumnumbers/all_symmetries.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/mpoperatorlist.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include "common/environment.h"
#include <boost/algorithm/string.hpp>

int main(int argc, char** argv)
{
   if (argc <= 1)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: mp-attr wavefunction [attribute [= value]] ...\n";
      return 1;
   }

   typedef std::map<std::string, std::string> AttribToSetType;
   typedef AttribToSetType::const_iterator AttribSetIter;
   std::map<std::string, std::string> AttributesToSet;

   typedef std::vector<std::string>::const_iterator AttribPrintIter;
   std::vector<std::string> AttributesToPrint;

   int arg = 2;
   while (arg < argc)
   {
      std::string a = argv[arg++];
      std::string Attrib;
      std::string::const_iterator Delim = std::find(a.begin(), a.end(), '=');
      if (Delim == a.end() && (arg >= argc || argv[arg][0] != '='))
      {
         Attrib = a;
         boost::trim(Attrib);
      }
      else
      {
         // attribute = value
         Attrib = std::string(static_cast<std::string const&>(a).begin(), Delim);
         boost::trim(Attrib);
         if (Delim == a.end())
         {
            // in this case, the first character of the next argument must be '='
            a = argv[arg++];
            Delim = a.begin();
         }
         ++Delim;
         if (Delim == a.end() && arg < argc)
         {
            a = argv[arg++];
            Delim = a.begin();
         }
         std::string Value(Delim, static_cast<std::string const&>(a).end());
         boost::trim(Value);
         AttributesToSet[Attrib] = Value;
      }
      AttributesToPrint.push_back(Attrib);
   }

   bool Readonly = AttributesToSet.empty();
   pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(argv[1], 
                                                          getenv_or_default("MP_CACHESIZE", 8*1024*1024),
                                                          Readonly);

   if (!AttributesToSet.empty())
   {
      pvalue_lock<MPWavefunction> PsiLock(Psi.lock());
      for (AttribSetIter I = AttributesToSet.begin(); I != AttributesToSet.end(); ++I)
         PsiLock->Attributes()[I->first] = I->second;
   }
   
   if (AttributesToPrint.size() == 0)
   {
      std::cout << Psi->Attributes();
   }

   // if there is only one attribute to print, show it without the Attrib=
   if (AttributesToSet.empty() && AttributesToPrint.size() == 1)
   {
      std::cout << Psi->Attributes()[AttributesToPrint[0]] << '\n';
   }
   else
   {
      for (AttribPrintIter I = AttributesToPrint.begin(); I != AttributesToPrint.end(); ++I)
      {
         std::cout << (*I) << "=" << Psi->Attributes()[*I] << '\n';
      }
   }

   if (!Readonly)
   {
      pheap::ShutdownPersistent(Psi);
   }
}
