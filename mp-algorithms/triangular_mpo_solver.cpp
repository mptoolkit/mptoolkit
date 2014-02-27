
#include "triangular_mpo_solver.h"


template <typename T>
std::complex<double>
FindLargestEigenvalue(MatrixOperator& M, T Func, bool Verbose)
{
   int Iterations = 20;
   double Tol = 1E-14;
   std::complex<double> EtaL;
   EtaL = LinearSolvers::Arnoldi(M, Func, Iterations, Tol, LinearSolvers::LargestMagnitude, Verbose);
   while (Iterations == 20)
   {
      Tol = 1E-14;
      EtaL = LinearSolvers::Arnoldi(M, Func, Iterations, Tol, LinearSolvers::LargestMagnitude, Verbose);
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
	 if (!is_scalar(I->second.TransformsAs()))
	    continue;

	 std::complex<double> Overlap = inner_prod(I->second, UnitMatrixRight);
	 I->second -= Overlap*UnitMatrixLeft;
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
			    MPOperator const& Diag, 
			    MatrixOperator const& UnitMatrixLeft, 
			    MatrixOperator const& UnitMatrixRight, 
			    LinearWavefunction const& Psi,
			    QuantumNumber const& QShift,
			    bool Magnitude1)
{
   // Components perpendicular to the identity satisfy equation (24)
   KMatrixPolyType E;
   for (KMatrixPolyType::const_iterator I = C.begin(); I != C.end(); ++I) // sum over momenta
   {
      std::complex<double> K = I->first;  // the momentum (complex phase)
      DEBUG_TRACE("Momentum")(K);
      for (int m = I->second.degree(); m >= 0; --m)
      {
	 DEBUG_TRACE("degree")(m);
	 MatrixOperator Rhs = conj(K) * I->second[m];
	 for (int k = m+1; k <= I->second.degree(); ++k)
	 {
	    Rhs -= double(Binomial(k,m)) * E[K][k];
	 }
            
	 // orthogonalize Rhs against the identity again, which is a kind of
	 // DGKS correction
	 if (is_scalar(Rhs.TransformsAs()))
	    Rhs -= inner_prod(UnitMatrixRight, Rhs) * UnitMatrixLeft;

	 double RhsNorm2 = norm_frob_sq(Rhs);
	 RhsNorm2 = RhsNorm2 / (Rhs.Basis1().total_degree()*Rhs.Basis2().total_degree());
	 //TRACE(RhsNorm2);
	 if (RhsNorm2 > 1E-22)
	 {
	    E[K][m] = LinearSolve(OneMinusTransfer(K*Diag, Psi, QShift, UnitMatrixRight, UnitMatrixLeft, Magnitude1), 
				  Rhs);
	    // do another orthogonalization
	    if (Magnitude1 && is_scalar(E[K][m].TransformsAs()))
	       E[K][m] -= inner_prod(UnitMatrixRight, E[K][m]) * UnitMatrixLeft;
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
	    E[K][i] = conj(K) * I->second[i];
	 for (int j = i+1; j <= MaxDegree; ++j)
	 {
	    E[K][i] -= double(Binomial(j,i)) * E[K][j];
	 }
      }
   }
   return E;
}

// Solve an MPO in the left-handed sense, as x_L * Op = lambda * x_L
// We currently assume there is only one eigenvalue 1 of the transfer operator
KMatrixPolyType
SolveMPO_Left(LinearWavefunction const& Psi, QuantumNumber const& QShift,
              TriangularOperator const& Op, MatrixOperator const& Rho,
              MatrixOperator const& Identity, bool Verbose = false)
{
   int Dim = Op.Basis1().size();       // dimension of the MPO
   std::vector<KMatrixPolyType> EMatK(Dim);  // the vector of E matrices

   // Initialize the first E matrix
   EMatK[Dim-1][1.0] = MatrixPolyType(Identity);

   // solve recursively
   int Col = Dim-2;
   while (Col >= 0)
   {
      if (Verbose)
      {
         std::cerr << "Solving column " << Col << '\n';
      }

      // Generate the next C matrices, C(n) = sum_{j>Col} Op(j,Col) E_j(n)
      KMatrixPolyType C;

#if 0
      MPOperator M = extract_lower_column(Op, Col);
      C = apply_right(EMatK, Psi, M, Psi)[0];
#else
      std::vector<std::vector<int> > Mask;
      mask_lower_column(Op, Col, Mask);
      C = apply_right_mask(EMatK, Psi, QShift, Op.data(), Psi, Mask)[Col];
#endif

      // Now do the classification, based on the properties of the diagonal operator
      MPOperator Diag = Op(Col, Col);
      MPOperatorClassification Classification = classify(Diag);
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
	 MatrixOperator UnitMatrixLeft = Identity;
	 MatrixOperator UnitMatrixRight = Rho;

	 // If the diagonal operator is unitary, it might have an eigenvalue of magnitude 1
	 if (Classification.is_unitary() && !Classification.is_complex_identity())
	 {
	    if (Verbose)
	       std::cerr << "Non-trivial unitary at column " << Col << std::endl;
	    
	    // Find the largest eigenvalue
	    // We need initial guess vectors in the correct symmetry sector
	    UnitMatrixLeft = MakeRandomMatrixOperator(Identity.Basis1(), Identity.Basis2(), Diag.Basis2()[0]);

	    std::complex<double> EtaL = FindLargestEigenvalue(UnitMatrixLeft, 
							      ApplyLeftQShift(Diag, Psi, QShift), Verbose);

	    if (Verbose)
	       std::cerr << "Eigenvalue of unitary operator is " << EtaL << std::endl;

	    if (std::abs(norm_frob(EtaL) - 1.0) < 1E-12)
	    {

	       UnitMatrixRight = UnitMatrixLeft;
	       // we have an eigenvalue of magnitude 1.  Find the right eigenvalue too
	       std::complex<double> EtaR = FindLargestEigenvalue(UnitMatrixRight, 
								 ApplyRightQShift(Diag, Psi, QShift), Verbose);
	       CHECK(norm_frob(EtaL-EtaR) < 1E-12)("Left and right eigenvalues do not agree!")(EtaL)(EtaR);
	       Factor = EtaL / norm_frob(EtaL);
	    }
	 }

	 // If we have an eigenvalue with magnitude 1, then decompose C into parallel and perpendicular parts
	 bool Magnitude1 = false;
         if (std::abs(norm_frob(Factor) - 1.0) < 1E-12)
         {
	    Magnitude1 = true;
	    EParallel = DecomposeParallelParts(C, Factor, UnitMatrixLeft, UnitMatrixRight);
	 }

         // Now the remaining components, which is anything that is not proportional
	 // to an eigenvector of magnitude 1
      
         KMatrixPolyType E = DecomposePerpendicularParts(C, Diag, UnitMatrixLeft, UnitMatrixRight, 
							 Psi, QShift, Magnitude1);

         // Reinsert the components parallel to the unit matrix (if any)
         for (KComplexPolyType::const_iterator I = EParallel.begin(); I != EParallel.end(); ++I)
         {
	    for (ComplexPolyType::const_iterator J = I->second.begin(); J != I->second.end(); ++J)
	    {
	       E[I->first][J->first] += J->second * Identity;
	    }
         }

	 // Finally, set the E matrix element at this column
	 DEBUG_TRACE(E[1.0]);
         EMatK[Col] = E;
      }

      --Col;
   }

   // The return value is column 0 of our E matrices
   return EMatK[0];
}
