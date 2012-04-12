// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "pheap/pvalueptr.h"

typedef pvalue_ptr<std::pair<int, int> > ppair;

int main()
{
   pheap::PHeapLog.SetStream(std::cerr);  
   pheap::PHeapLog.SetThreshold(100);
   
   pheap::Initialize("testdata", 1, 1024*1024, 16*1024*1024);

   ppair N = new std::pair<int, int>(100, 200);

   pheap::ShutdownPersistent(N);
}

