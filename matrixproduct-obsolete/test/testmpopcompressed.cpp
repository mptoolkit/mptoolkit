// -*- C++ -*- $Id$

#include "models/hubbard-so4.h"
#include "matrixproduct/mpopcompressed.h"
#include "common/trace.h"

template <typename T>
struct PrintRunLengthCompressed : public boost::static_visitor<>
{
   PrintRunLengthCompressed(std::ostream& out_, int NestLevel_) : out(out_), NestLevel(NestLevel_) {}

   void DoNest() const
   {
      out << std::string(NestLevel*4, ' ');
   }

   void operator()(T const& x) const
   {
      this->DoNest();
      out << "Value:\n";
      this->DoNest();
      out << x << '\n';
   };

   void operator()(run_length_array<T> const& x) const
   {
      this->DoNest();
      out << "Array: size " << x.size() << '\n';
      int i = 0;
      for (typename run_length_array<T>::const_iterator I = x.begin(); I != x.end(); ++I)
      {
         this->DoNest();
         out << "Array element " << i++ << '\n';
         I->apply_visitor(PrintRunLengthCompressed<T>(out, NestLevel+1));
      }
   }

   void operator()(run_length_repeat<T> const& x) const
   {
      this->DoNest();
      out << "Repeat: size " << x.size() << '\n';
      x.nested().apply_visitor(PrintRunLengthCompressed<T>(out, NestLevel+1));
   }

   std::ostream& out;
   int NestLevel;
};

template <typename T>
std::ostream& operator<<(std::ostream& out, run_length_compressed<T> const& x)
{
   x.apply_visitor(PrintRunLengthCompressed<T>(out, 0));
   return out;
}

int main()
{
   SymmetryList Symmetry("Q:SU(2),S:SU(2)");
   QuantumNumbers::QNConstructor<QuantumNumbers::SU2,QuantumNumbers::SU2> QN(Symmetry);

   SiteBlock A = CreateSO4HubbardSiteA();
   SiteBlock B = CreateSO4HubbardSiteB();

   Lattice L = join(Lattice(A), Lattice(B));
   L = repeat(L, 6);

   MPOpCompressed C0 = CreateMPOperator(L, "C",6);
#if 0
   //TRACE(C0.size())(C0);

   MPOpCompressed CH1 = CreateMPOperator(L, "C", 7);

   ProductBasis<BasisList, BasisList> PB(C0.front().Basis1(), CH1.front().Basis1(), QN(0,0));
   SimpleOperator C(PB.Basis(), QN(0,0));
   C(0,0) = 1;

   TRACE(C0)(CH1);

   MPOpCompressed H = do_prod(C, PB, C0, CH1);
   TRACE(C)(H);

   SumBasis<BasisList> sb(C0.front().Basis1(), CH1.front().Basis1());
   SimpleOperator S(C0.front().Basis1(), sb.Basis());
   S(0,0) = 1.0;
   S(0,1) = 1.0;

   MPOpCompressed SS = do_sum(S, sb, C0, CH1);
   TRACE(S)(SS);

   SimpleOperator STrunc(S.Basis2(), make_vacuum_basis(S.GetSymmetryList()));
   for (unsigned i = 0; i < S.Basis2().size(); ++i)
   {
      STrunc(i,0) = 1.0;
   }

   TRACE(STrunc);

   MPOpCompressed SSTrunc = inject_right(SS, STrunc);
   TRACE(STrunc)(SSTrunc);

   MPOpCompressed SSSTrunc = inject_left(STrunc, SSTrunc);
   TRACE(STrunc)(SSSTrunc);
#endif

   SumBasis<BasisList> Cxsb(C0.front().Basis1(), C0.front().Basis1());
   SimpleOperator Cx(C0.front().Basis1(), Cxsb.Basis());
   Cx(0,0) = 1.0;
   Cx(0,1) = 1.0;

   MPOpCompressed SumTest = do_sum(Cx, Cxsb, C0, C0);
   SimpleOperator Cxr(Cx.Basis2(), C0.back().Basis2());
   Cxr(0,0) = 1.0;
   Cxr(1,0) = 1.0;
   Cx = Cx * Cxr;
   SumTest = inject_right(SumTest, Cx);
   TRACE(C0)(Cx)(SumTest);
   
}
