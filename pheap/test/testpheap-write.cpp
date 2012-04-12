// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "pheap/pvalueptr.h"

typedef std::pair<long, long> mypair;
typedef pvalue_ptr<mypair> ppair;

int main()
{
   pheap::PHeapLog.SetStream(std::cerr);  
   pheap::PHeapLog.SetThreshold(100);
   
   pheap::Initialize("testdata", 1, 16384, 16*1024*1024);

   ppair N = new mypair(100, 200);

   pheap::ShutdownPersistent(N);
}

