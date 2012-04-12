// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "pheap/pvalueptr.h"

typedef std::pair<long, long> mypair;
typedef pvalue_ptr<mypair> ppair;

int main()
{
   pheap::PHeapLog.SetStream(std::cerr);  
   pheap::PHeapLog.SetThreshold(100);
   
   ppair N = pheap::OpenPersistent("testdata", 16*1024*1024);

   TRACE(N->first)(N->second);

   ppair M = new mypair(*N);
   pheap::ShutdownPersistent(M);
}

