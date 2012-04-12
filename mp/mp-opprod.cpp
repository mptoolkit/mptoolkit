// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include "mp/copyright.h"

int main(int argc, char** argv)
{
   if (argc < 5 || argc > 6)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: mp-opprod <lattice> <operator1> <operator2> <out-operator> [quantum-number]\n";
      return 1;
   }

   pvalue_ptr<OperatorList> System = pheap::OpenPersistent(argv[1], 655360);

   MPOperator Op1 = (*System)[argv[2]];
   MPOperator Op2 = (*System)[argv[3]];
   std::string OutOperator = argv[4];

   QuantumNumbers::QuantumNumber Q;
   if (argc == 6)
   {
      Q = QuantumNumbers::QuantumNumber(Op1.GetSymmetryList(), argv[5]);
   }
   else
   {
      QuantumNumbers::QuantumNumberList QList = transform_targets(Op1.TransformsAs(), Op2.TransformsAs());
      if (QList.size() != 1)
      {
         std::cerr << "quantum number of operator-operator product is ambiguous,\n"
                   << "  choices are ";
         std::copy(QList.begin(), QList.end(), std::ostream_iterator<QuantumNumber>(std::cerr, " "));
         std::cerr << '\n';
         exit(1);
      }
      Q = QList[0];
   }

   MPOperator Out = prod(Op1, Op2, Q);
   (*System.mutate())[OutOperator] = Out;

   std::cout << "Created operator: " << OutOperator << std::endl;

   pheap::ShutdownPersistent(System);
}
