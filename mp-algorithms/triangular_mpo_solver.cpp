// -*- C++ -*- $Id$

#include "triangular_mpo_solver.h"
#include "mps/momentum_operations.h"
#include "mps/operator_actions.h"
#include "mp-algorithms/arnoldi.h"
#include "mp-algorithms/gmres.h"

// For identifying an eigenvalue 1 of the transfer matrix, we need
// an epsilon tolerance.  1E-12 proved to be a bit too small in some
// cases.
double const EigenUnityEpsilon = 1E-10;

// A note on orthogonalizing vectors in non-orthogonal Hilbert spaces.
// Suppose <l| and |r> are left and right eigenvectors, satisfying <l|r>=1.
// If we want to orthogonalize some bra <x| against <l|, we want do do:
// <x'| = <x| - conj(c)*<l|
// where c = <x|r>
// Proof: <x'|r> = <x|r> - c*<l|r> = 0
//
// In the other direction, if we want to orthogonalize |y> against |r>, then
// |y'> = |y> - c*|r>
// where c = <l|y>

// binomial function.  This can be implemented with integer arithmetic.
// gmp contains a version with arbitrary size integers, but we're never going
// to overflow a long (or even an int!) with an MPO contraction.
long Binomial(int n, int k)
{
   if (k > n/2)
      k = n-k;     // take advantage of symmetry
   long r = 1;
   long n_k = n-k;
   for (long i = 1; i <= k; ++i)
   {
      DEBUG_CHECK((r*(n_k-i))%i == 0);
      r = (r * (n_k+i)) / i;
   }
   return r;
}

// Functor to evaluate (1-T)(x) where T is the generalized transfer matrix.
// x is defined in the Basis1.  Optionally, we also optionally orthogonalize against a
// Unit matrix, which is an eigenvalue 1 of T
//      -Psi*-
//        |
//  T =  Op*
//        |
//      -Psi--
struct OneMinusTransferLeft
{
   OneMinusTransferLeft(FiniteMPO const& Op, LinearWavefunction const& Psi, QuantumNumber const& QShift, 
			MatrixOperator const& LeftUnit,
			MatrixOperator const& RightUnit, bool Orthogonalize)
      : Op_(Op), Psi_(Psi), 
	QShift_(QShift), LeftUnit_(LeftUnit), 
	RightUnit_(RightUnit), Orthogonalize_(Orthogonalize) { }

   MatrixOperator operator()(MatrixOperator const& x) const
   {
      MatrixOperator r = x-delta_shift(inject_left(x, Op_, Psi_), QShift_);
      if (Orthogonalize_ && r.TransformsAs() == RightUnit_.TransformsAs())
	 {
	    DEBUG_TRACE(inner_prod(r, RightUnit_))("this should be small");
	    DEBUG_TRACE(inner_prod(LeftUnit_, r));
	    r -= conj(inner_prod(r, RightUnit_)) * LeftUnit_; // orthogonalize to the identity
	    DEBUG_TRACE(inner_prod(r, RightUnit_))("this should be zero");
	 }
      return r;
   }

   FiniteMPO const& Op_;
   LinearWavefunction const& Psi_;
   QuantumNumber const& QShift_;
   MatrixOperator const& LeftUnit_;
   MatrixOperator const& RightUnit_;
   bool Orthogonalize_;
};

template <typename Func>
MatrixOperator
LinearSolve(Func F, MatrixOperator Rhs, int Verbose = 0)
{
   MatrixOperator Guess = MakeRandomMatrixOperator(Rhs.Basis1(), Rhs.Basis2(), Rhs.TransformsAs());
   //Rhs;
   int m = 30;
   int max_iter = 10000;
   double tol = 1e-12;
   GmRes(Guess, F, Rhs, m, max_iter, tol, LinearAlgebra::Identity<MatrixOperator>(), Verbose);
   return Guess;
}


// Calculate the (complex) eigenvalue that is closest to 1.0
// using Arnoldi.
template <typename T>
std::complex<double>
FindClosestUnitEigenvalue(MatrixOperator& M, T Func, int Verbose)
{
   int Iterations = 20;
   double Tol = 1E-16;
   std::complex<double> EtaL;
   EtaL = LinearSolvers::Arnoldi(M, Func, Iterations, Tol, LinearSolvers::LargestAlgebraicReal, Verbose);
   while (Iterations == 20)
   {
      Tol = 1E-14;
      EtaL = LinearSolvers::Arnoldi(M, Func, Iterations, Tol, LinearSolvers::LargestAlgebraicReal, Verbose);
   }
   return EtaL;
}

KComplexPolyType
DecomposeParallelParts(KMatrixPolyType& C, std::complex<double> Factor, 
		       MatrixOperator const& UnitMatrixLeft, 
		       MatrixOperator const& UnitMatrixRight)
{
   KComplexPolyType EParallel;
   // diagonal element is the identity, up to a unitary factor
   DEBUG_TRACE("Unit diagonal element")(Factor);
   
   // decompose C into components parallel and perpendicular to the identity
   // The only part we have to care about is a component with the same momentum as our unit operator
   for (KMatrixPolyType::iterator Ck = C.begin(); Ck != C.end(); ++Ck) // sum over momenta
   {
      ComplexPolyType CParallel;
      std::complex<double> K = Ck->first;  // the momentum (complex phase)
      for (MatrixPolyType::iterator I = Ck->second.begin(); I != Ck->second.end(); ++I)
      {
	 if (I->second.is_null())
	 {
	    TRACE(Ck->second);
	    PANIC("oops");
	 }

	 //if (!is_scalar(I->second.TransformsAs()))
	 //continue;

	 std::complex<double> Overlap = inner_prod(I->second, UnitMatrixRight);
	 DEBUG_TRACE(Overlap);
	 I->second -= conj(Overlap)*UnitMatrixLeft;
	 DEBUG_TRACE(inner_prod(I->second, UnitMatrixRight))("should be zero");
	 DEBUG_TRACE(inner_prod(UnitMatrixLeft, UnitMatrixRight));
	 if (norm_frob(Overlap) > 1E-16)
	    CParallel[I->first] = Overlap;
      }

      // Is this the same momentum as our unit operator?
      if (norm_frob(K - Factor) < 1E-12)
      {
	 DEBUG_TRACE("Component at equal momenta")(K);
	 // same momenta, these components diverge		  
	 for (int m = CParallel.degree(); m >= 0; --m)
	 {
	    if (CParallel.has_term(m))
	    {
	       EParallel[Factor][m+1] = conj(Factor) * CParallel[m]; // include momentum
	       for (int k = m+2; k <= CParallel.degree()+1; ++k)
	       {
		  if (EParallel[Factor].has_term(k))
		  {
		     EParallel[Factor][m+1] -= double(Binomial(k,m)) 
			* EParallel[Factor][k];
		  }
	       }
	       EParallel[Factor][m+1] *= 1.0 / (1.0 + m);
	    }
	 }
      }
      else
      {
	 for (int m = CParallel.degree(); m >= 0; --m)
	 {
	    // different momenta, we get a contribution both at K and Factor
	    if (CParallel.has_term(m))
	    {
	       std::complex<double> Term = CParallel[m] / (K - Factor);
	       EParallel[K][m] += Term;
	       EParallel[Factor][m] -= Term;
	    }
	 }
      }
   }
   return EParallel;
}

KMatrixPolyType
DecomposePerpendicularParts(KMatrixPolyType& C,
			    FiniteMPO const& Diag, 
			    MatrixOperator const& UnitMatrixLeft, 
			    MatrixOperator const& UnitMatrixRight, 
			    LinearWavefunction const& Psi,
			    QuantumNumber const& QShift,
			    bool HasEigenvalue1,
			    int Verbose)
{
   // Components perpendicular to the identity satisfy equation (24)
   KMatrixPolyType E;
   for (KMatrixPolyType::const_iterator I = C.begin(); I != C.end(); ++I) // sum over momenta
   {
      std::complex<double> K = I->first;  // the momentum (complex phase)

      if (Verbose)
	 std::cerr << "Momentum " << K << std::endl;

      DEBUG_TRACE("Momentum")(K)(HasEigenvalue1);
      for (int m = I->second.degree(); m >= 0; --m)
      {
	 if (Verbose)
	    std::cerr << "Degree " << m << std::endl;

	 DEBUG_TRACE("degree")(m);
	 MatrixOperator Rhs = conj(K) * I->second[m];
	 for (int k = m+1; k <= I->second.degree(); ++k)
	 {
	    // avoid accessing E[K][k] if it doesn't exist, to avoid adding a null term
	    if (E.has_element(K) && E[K].has_term(k))
	       Rhs -= double(Binomial(k,m)) * E[K][k];
	 }
            
	 // orthogonalize Rhs against the identity again, which is a kind of
	 // DGKS correction
	 if (HasEigenvalue1 && Rhs.TransformsAs() == UnitMatrixRight.TransformsAs())
	 {
	    DEBUG_TRACE(inner_prod(Rhs, UnitMatrixRight))("should be small");
	    Rhs -= conj(inner_prod(Rhs, UnitMatrixRight)) * UnitMatrixLeft;
	    DEBUG_TRACE(inner_prod(Rhs, UnitMatrixRight))("should be zero");
	    DEBUG_TRACE(inner_prod(Rhs, UnitMatrixLeft));
	 }


	 double RhsNorm2 = norm_frob_sq(Rhs);
	 RhsNorm2 = RhsNorm2 / (Rhs.Basis1().total_degree()*Rhs.Basis2().total_degree());
	 DEBUG_TRACE(RhsNorm2);
	 DEBUG_TRACE(HasEigenvalue1)(UnitMatrixLeft)(UnitMatrixRight)(K)(Diag)(Rhs);
	 if (RhsNorm2 > 1E-22)
	 {
	    E[K][m] = LinearSolve(OneMinusTransferLeft(K*Diag, Psi, QShift, 
						       UnitMatrixLeft, UnitMatrixRight, HasEigenvalue1), 
				  Rhs, Verbose);

	    DEBUG_CHECK(!E[K][m].is_null());
	    // do another orthogonalization
	    if (HasEigenvalue1 && E[K][m].TransformsAs() == UnitMatrixRight.TransformsAs())
	    {
	       DEBUG_TRACE(inner_prod(E[K][m], UnitMatrixRight))("should be small");
	       E[K][m] -= conj(inner_prod(E[K][m], UnitMatrixRight)) * UnitMatrixLeft;
	       DEBUG_TRACE(inner_prod(E[K][m], UnitMatrixRight))("should be zero");
	    }
	 }
      }
   }
   return E;
}

// Solve the components for the case where the diagonal operator is zero
KMatrixPolyType 
SolveZeroDiagonal(KMatrixPolyType const& C)
{
   // See equation (20)
   // E(n) = C(n-1)
   // which expands to, for the degree d component of the polynomial with momentum k,
   // E^k_d = e^{-ik} C^k_d - \sum_{j=d+1}^p (j d) E^k_j
   KMatrixPolyType E;
   for (KMatrixPolyType::const_iterator I = C.begin(); I != C.end(); ++I) // sum over momenta
   {
      std::complex<double> K = I->first;  // the momentum (complex phase)
      int MaxDegree = I->second.degree();
      for (int i = MaxDegree; i >= 0; --i)
      {
	 if (I->second.has_term(i))
	 {
	    DEBUG_CHECK(!I->second[i].is_null());
	    E[K][i] = conj(K) * I->second[i];
	 }
	 for (int j = i+1; j <= MaxDegree; ++j)
	 {
	    DEBUG_CHECK(!E[K][j].is_null());
	    E[K][i] -= double(Binomial(j,i)) * E.lookup_or_default(K)[j];
	 }
      }
   }
   return E;
}

// Solve an MPO in the left-handed sense, as x_L * Op = lambda * x_L
// We currently assume there is only one eigenvalue 1 of the transfer operator

KMatrixPolyType
SolveMPO_Left(LinearWavefunction const& Psi, QuantumNumber const& QShift,
              TriangularMPO const& Op, MatrixOperator const& LeftIdentity,
              MatrixOperator const& RightIdentity, int Verbose)
{
   CHECK_EQUAL(RightIdentity.Basis1(), Psi.Basis1());
   CHECK_EQUAL(RightIdentity.Basis2(), Psi.Basis1());
   CHECK_EQUAL(LeftIdentity.Basis1(), Psi.Basis1());
   CHECK_EQUAL(LeftIdentity.Basis2(), Psi.Basis1());

   int Dim = Op.Basis1().size();       // dimension of the MPO
   std::vector<KMatrixPolyType> EMatK(Dim);  // the vector of E matrices

   if (Verbose)
      std::cerr << "SolveMPO_Left: dimension is " << Dim << std::endl;

   // solve recursively column 0 onwards

   // Make sure the (0,0) part is identity
   OperatorClassification CheckIdent = classify(Op(0,0));
   CHECK(CheckIdent.is_identity())("(0,0) component of the MPO must be identity!")(CheckIdent)(Op(0,0));

   // Initialize the first E matrix.  These are operators acting in the Basis1()
   EMatK[0][1.0] = MatrixPolyType(LeftIdentity);

   for (int Col = 1; Col < Dim; ++Col)
   {
      if (Verbose)
      {
         std::cerr << "Solving column " << Col << '\n';
      }

      // Generate the next C matrices, C(n) = sum_{j<Col} Op(j,Col) E_j(n)
      KMatrixPolyType C;

      std::vector<std::vector<int> > Mask = mask_column(Op, Col);
      C = inject_left_mask(EMatK, Psi, QShift, Op.data(), Psi, Mask)[Col];

      // Now do the classification, based on the properties of the diagonal operator
      FiniteMPO Diag = Op(Col, Col);
      OperatorClassification Classification = classify(Diag);
      if (Classification.is_null())
      {
         DEBUG_TRACE("Zero diagonal element")(Col)(Diag);
	 if (Verbose)
	    std::cerr << "Zero diagonal matrix element at column " << Col << std::endl;
         EMatK[Col] = SolveZeroDiagonal(C);
      }
      else
      {
	 if (Verbose)
	    std::cerr << "Non-zero diagonal matrix element at column " << Col << std::endl;

	 // Non-zero diagonal element.
	 // In this case we have to consider a possible component with eigenvalue magnitude 1.
	 // We call this component (if it exists) the unit matrix.

         KComplexPolyType EParallel;  // components parallel to the identity at some momentum, may be zero

	 // Multiplication factor and the left and right eigenvalues.
	 std::complex<double> Factor = Classification.factor();
	 MatrixOperator UnitMatrixLeft = LeftIdentity;
	 MatrixOperator UnitMatrixRight = RightIdentity; //delta_shift(RightIdentity, QShift);

	 // If the diagonal operator is unitary, it might have an eigenvalue of magnitude 1
	 // Or if we are not in the scalar sector, then there might be an additional eigenvalue 1
	 // due to symmetry.  We don't attempt to handle the case where we have more than one
	 // eigenvalue 1 in the same symmetry sector.
	 if (Classification.is_unitary() && (!is_scalar(Diag.Basis2()[0]) || !Classification.is_complex_identity()))
	 {
	    if (Verbose)
	       std::cerr << "Non-trivial unitary at column " << Col 
			 << ", finding left eigenvalue..." << std::endl;
	    
	    // Find the largest eigenvalue
	    // We need initial guess vectors in the correct symmetry sector
	    UnitMatrixLeft = MakeRandomMatrixOperator(LeftIdentity.Basis1(), LeftIdentity.Basis2(), Diag.Basis2()[0]);

	    std::complex<double> EtaL = FindClosestUnitEigenvalue(UnitMatrixLeft, 
								  InjectLeftQShift(Diag, Psi, QShift), Verbose);

	    if (Verbose)
	       std::cerr << "Eigenvalue of unitary operator is " << EtaL << std::endl;

	    if (std::abs(norm_frob(EtaL) - 1.0) < EigenUnityEpsilon)
	    {
	       if (Verbose)
		  std::cerr << "Found an eigenvalue 1, so we need the right eigenvector..." << std::endl;

	       UnitMatrixRight = UnitMatrixLeft;
	       // we have an eigenvalue of magnitude 1.  Find the right eigenvalue too
	       std::complex<double> EtaR = FindClosestUnitEigenvalue(UnitMatrixRight, 
								     InjectRightQShift(Diag, Psi, QShift), Verbose);
	       if (Verbose)
		  std::cerr << "Right eigenvalue is " << EtaR << std::endl;

	       CHECK(norm_frob(EtaL-EtaR) < 1E-12)("Left and right eigenvalues do not agree!")(EtaL)(EtaR);
	       // we already determined that the norm is sufficiently close to 1, but 
	       // fine-tune normalization
	       Factor = EtaL / norm_frob(EtaL);

	       // normalize the left/right eigenvector pair
	       //	       UnitMatrixLeft *= 1.0 / (inner_prod(UnitMatrixRight, UnitMatrixLeft));
	       UnitMatrixRight *= 1.0 / (inner_prod(UnitMatrixLeft, UnitMatrixRight));
	       DEBUG_TRACE(inner_prod(UnitMatrixLeft, UnitMatrixRight));
	    }
	    else if (std::abs(norm_frob(EtaL) - 1.0) < EigenUnityEpsilon*100)
	    {
	       // If we get here, then we have an eigenvalue that is close to 1, but
	       // misses the test
	       std::cerr << "SolveMPO_Left: warning: Eigenvalue misses tolerance but is near 1"
			 << ", epsilon=" << std::abs(norm_frob(EtaL)-1.0) << '\n';
	    }
	 }

	 // If we have an eigenvalue equal to 1, then decompose C into parallel and perpendicular parts
	 bool HasEigenvalue1 = false;
         if (std::abs(norm_frob(Factor) - 1.0) < 1E-12)
         {
	    HasEigenvalue1 = true;
	    //DEBUG_TRACE(UnitMatrixLeft)(UnitMatrixRight);
	    if (Verbose)
	       std::cout << "Decomposing parts parallel to the unit matrix" << std::endl;
	    EParallel = DecomposeParallelParts(C, Factor, UnitMatrixLeft, UnitMatrixRight);
	 }

         // Now the remaining components, which is anything that is not proportional
	 // to an eigenvector of magnitude 1
      
	 if (Verbose)
	    std::cout << "Decomposing parts perpendicular to the unit matrix" << std::endl;
         KMatrixPolyType E = DecomposePerpendicularParts(C, Diag, UnitMatrixLeft, UnitMatrixRight, 
							 Psi, QShift, HasEigenvalue1, Verbose);

         // Reinsert the components parallel to the unit matrix (if any)
         for (KComplexPolyType::const_iterator I = EParallel.begin(); I != EParallel.end(); ++I)
         {
	    for (ComplexPolyType::const_iterator J = I->second.begin(); J != I->second.end(); ++J)
	    {
	       // Conj here because this comes from an overlap(x, RightUnitMatrix)
	       E[I->first][J->first] += conj(J->second) * UnitMatrixLeft;
	       DEBUG_TRACE(inner_prod(E[I->first][J->first], UnitMatrixRight));
	    }
         }

	 // Finally, set the E matrix element at this column
	 //DEBUG_TRACE(E[1.0]);
         EMatK[Col] = E;
      }

   }

   // The return value is the last column of our E matrices
   return EMatK.back();
}
