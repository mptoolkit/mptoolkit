// ENDHEADER

#include "infinitewavefunctionleft.h"
#include "linearalgebra/arpack_wrapper.h"
#include "mps/packunpack.h"

std::tuple<std::complex<double>, int>
overlap_arpack(InfiniteWavefunctionLeft const& x, ProductMPO const& StringOp,
	       InfiniteWavefunctionLeft const& y,
	       QuantumNumbers::QuantumNumber const& Sector, int Iter, double Tol, int Verbose)
{

   int Length = statistics::lcm(x.size(), y.size(), StringOp.size());

   ProductMPO Str = StringOp * ProductMPO::make_identity(StringOp.LocalBasis2List(), Sector);

   LinearWavefunction xPsi = get_left_canonical(x).first;
   LinearWavefunction yPsi = get_left_canonical(y).first;

   StateComponent Init = MakeRandomStateComponent(Str.Basis1(), x.Basis1(), y.Basis1());

   int ncv = Iter;
   double tolsave = Tol;
   int ncvsave = ncv;
   ApplyToPackedStateComponent<LeftMultiplyOperator> Func(LeftMultiplyOperator(xPsi, x.qshift(), Str, yPsi, y.qshift(), Length),
							  Str.Basis1(), x.Basis1(), y.Basis1());
   int n = Func.pack_size();

   int NumEigen = 1;

   LinearAlgebra::Vector<std::complex<double>> RightEigen =
      LinearAlgebra::DiagonalizeARPACK(Func, n, NumEigen, Tol, nullptr, ncv, true, Verbose);

   return std::make_tuple(RightEigen[0], Length);
}

