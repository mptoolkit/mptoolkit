// -*- C++ -*-

#include "eigensolver.h"
#include "mp-algorithms/lanczos.h"

struct MPSMultiply
{
   MPSMultiply(StateComponent const& E_, OperatorComponent const& H_, StateComponent const& F_)
      : E(E_), H(H_), F(F_)
   {
   }

   StateComponent operator()(StateComponent const& Psi) const
   {
      StateComponent R = operator_prod_inner(H, E, Psi, herm(F));
      return R;
   }

   StateComponent const& E;
   OperatorComponent const& H;
   StateComponent const& F;
};

LocalEigensolver::LocalEigensolver()
   : FidelityScale(0.1), MaxTol(1e-4), MinTol(1-10), MinIter(2), MaxIter(20), Verbose(0)
{
}

void
LocalEigensolver::SetInitialFidelity(int UnitCellSize, double f)
{
   FidelityAv_ = statistics::moving_exponential<double>(exp(log(0.25)/UnitCellSize));
   FidelityAv_.push(f);
}

double
LocalEigensolver::Solve(StateComponent& C,
			StateComponent const& LeftBlockHam,
			OperatorComponent const& H,
			StateComponent const& RightBlockHam)
{
   DEBUG_CHECK_EQUAL(C.Basis1(), LeftBlockHam.Basis2());
   DEBUG_CHECK_EQUAL(C.Basis2(), RightBlockHam.Basis1());
   DEBUG_CHECK_EQUAL(LeftBlockHam.LocalBasis(), H.Basis1());
   DEBUG_CHECK_EQUAL(RightBlockHam.LocalBasis(), H.Basis2());
   DEBUG_CHECK_EQUAL(C.LocalBasis(), H.LocalBasis2());

   StateComponent ROld = C;

   if (EvolveDelta == 0.0)
   {
      LastTol_ = std::min(std::sqrt(this->AverageFidelity()) * FidelityScale, MaxTol);
      LastTol_ = std::max(LastTol_, MinTol);
      //LastTol_ = std::min(this->AverageFidelity() * FidelityScale, MaxTol);
      LastIter_ = MaxIter;
      if (Verbose > 2)
      {
	 std::cerr << "Starting eigensolver.  Initial guess vector has dimensions "
		   << C.Basis1().total_dimension() << " x " << C.LocalBasis().size()
		   << " x " << C.Basis2().total_dimension() << '\n';
      }
      LastEnergy_ = Lanczos(C, MPSMultiply(LeftBlockHam, H, RightBlockHam),
			    LastIter_, LastTol_, MinIter, Verbose-1);
      
   }
   else
   {
      C = operator_prod_inner(H, LeftBlockHam, ROld, herm(RightBlockHam));
      LastEnergy_ = inner_prod(ROld, C).real();
      C = ROld - EvolveDelta * C; // imaginary time evolution step
      C *= 1.0 / norm_frob(C);    // normalize
      LastIter_ = 1;
      LastTol_ = 0.0;
   }

   LastFidelity_ = std::max(1.0 - norm_frob(inner_prod(ROld, C)), 0.0);
   FidelityAv_.push(LastFidelity_);
   return LastEnergy_;
}
