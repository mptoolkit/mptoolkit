// -*- C++ -*- $Id$
//
// CanonicalWavefunction: class to represent an MPS that is in (left or right) canonical form.
//
// The canonical form as all matrices in left (right) canonical form, and we also store
// the Lambda matrix (which is diagonal) at each partition.
//
// a LeftCanonicalWavefunction has all MPS matrices in left-canonical form and the corresponding
// lambda is on the right (Basis2) at each site.
//
// A RightCanonicalWavefunction has all MPS matrices in right-canonical form and the corresponding
// lambda matrix is on the left (Basis1) at each site.
//
// A CanonicalWavefunction is generally going to be read-only.  Non-const
// iterators and access is defined, but generally we would want to convert to
// a LinearWavefunction or a CenterWavefunction and use that instead.
//
// 

#if !defined(MPTOOLKIT_MPS_CANONICALWAVEFUNCTION_H)
#define MPTOOLKIT_MPS_CANONICALWAVEFUNCTION_H

#include "infinitewavefunction.h"
#include "linearwavefunction.h"
#include "linearalgebra/diagonalmatrix.h"

// The lambda matrix is real-valued and diagonal
typedef Tensor::IrredTensor<LinearAlgebra::DiagonalMatrix<double>,
			    VectorBasis, VectorBasis> RealDiagonalOperator;

// Common functions for Left and Right-orthogonalized versions
class CanonicalWavefunctionBase
{
   public:
      // storage for the MPS matrices
      typedef StateComponent                                        mps_type;
      typedef pvalue_handle<mps_type>                               mps_handle_type;
      typedef std::vector<mps_handle_type>                          mps_container_type;
      typedef mps_container_type::iterator                          base_mps_iterator;
      typedef mps_container_type::const_iterator                    const_base_mps_iterator;
      typedef pvalue_handle_iterator<base_mps_iterator>             mps_iterator;
      typedef const_pvalue_handle_iterator<const_base_mps_iterator> const_mps_iterator;

      // storage for the lambda matrix
      typedef RealDiagonalOperator lambda_type;
      typedef pvalue_handle<lambda_type> lambda_handle_type;
      typedef std::vector<lambda_handle_type>                          lambda_container_type;
      typedef lambda_container_type::iterator                          base_lambda_iterator;
      typedef lambda_container_type::const_iterator                    const_base_lambda_iterator;
      typedef pvalue_handle_iterator<base_lambda_iterator>             lambda_iterator;
      typedef const_pvalue_handle_iterator<const_base_lambda_iterator> const_lambda_iterator;

      int size() const { return MpsData.size(); }

      // Returns true if the state transforms irreducibly.  This is true
      // iff the left basis contains only a single site.
      bool is_irreducible() const;

      // returns true if the state can be considered to be a pure state; which 
      // implies that both the right hand basis is one-dimensional.
      // We can interpret a reducible left-hand basis also as a pure state, by summing
      // the states, so we do not require that the left-hand basis is one-dimensional.
      bool is_pure() const;

      // precondition: is_irreducible
      QuantumNumbers::QuantumNumber TransformsAs() const;

      // returns the left-most basis.
      VectorBasis Basis1() const;
      // returns the right-most basis.
      VectorBasis Basis2() const;

      mps_iterator begin() { return mps_iterator(Data.begin()); }
      mps_iterator end() { return mps_iterator(Data.end()); }

      const_mps_iterator begin() const { return const_mps_iterator(Data.begin()); }
      const_mps_iterator end() const { return const_mps_iterator(Data.end()); }

      lambda_iterator lambda_begin() { return lambda_iterator(Lambda.begin()); }
      lambda_iterator lambda_end() { return lambda_iterator(Lambda.end()); }

      const_lambda_iterator lambda_begin() const { return const_lambda_iterator(Lambda.begin()); }
      const_lambda_iterator lambda_end() const { return const_lambda_iterator(Lambda.end()); }

      // return the i'th MPS matrix.  Because they are stored by handle, we can't
      // return a reference, but the tensors are reference counted anyway so a copy is cheap
      mps_type operator[](int i) const;

      // returns the lambda matrix at site i (acting on [i].Basis2() for a LeftCanonicalWavefunction
      // or [i].Basis1() for a RightCanonicalWavefunction)
      lambda_type lambda(int i) const; 

   protected:
      // don't allow construction except via derived classes
      CanonicalWavefunctionBase();
      CanonicalWavefunctionBase(CanonicalWavefunctionBase const& Psi);

   private:
      void push_back(mps_type const& x);
      void push_back_lambda(lambda_type const& x);

      mps_container_type Data;
      lambda_container_type Lambda;
};

class LeftCanonicalWavefunction : public CanonicalWavefunctionBase
{
   public:
      // construction and orthogonalization from a LinearWavefunction
      explicit LeftCanonicalWavefunction(LinearWavefunction const& Psi);

      LeftCanonicalWavefunction(LeftCanonicalWavefunction const& Psi);

      LeftCanonicalWavefunction& operator=(LeftCanonicalWavefunction const& Psi);

      void check_structure() const;
      void debug_check_structure() const;

};

class RightCanonicalWavefunction : public CanonicalWavefunctionBase
{
   public:
      // construction and orthogonalization from a LinearWavefunction
      explicit RightCanonicalWavefunction(LinearWavefunction const& Psi);

      RightCanonicalWavefunction(RightCanonicalWavefunction const& Psi);

      RightCanonicalWavefunction& operator=(RightCanonicalWavefunction const& Psi);

      void check_structure() const;
      void debug_check_structure() const;
};

#endif
