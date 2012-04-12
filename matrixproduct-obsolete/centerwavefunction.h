// -*- C++ -*- $Id$
//
// CenterWavefunction: Wrapper around a LinearWavefunction,
// that adds a center matrix. 
//

#if !defined(MPWAVEFUNCTION_H_FUIYT49786Y709)
#define MPWAVEFUNCTION_H_FUIYT49786Y709

#include "linearwavefunction.h"
#include "splitoperator.h"

class CenterWavefunction
{
   public:
      typedef MPStateComponent ComponentType;
      typedef ComponentType::OperatorType OperatorType;
      typedef OperatorType::basis1_type MatrixBasisType;

      CenterWavefunction();

      // Construction from a LinearWavefunction
      CenterWavefunction(LinearWavefunction const& Psi);
      CenterWavefunction(MatrixOperator const& C, LinearWavefunction const& Psi);

      SymmetryList GetSymmetryList() const;

      // if this is irreducible, then returns the quantum number
      QuantumNumber TransformsAs() const;

      int LeftSize() const;
      int RightSize() const;

      bool is_null() const { return this->LookupLeft(0).Basis1().is_null(); }

      int size() const { return this->LeftSize() + this->RightSize(); }

      OperatorType const& Center() const { return CenterOp; }
      OperatorType& Center() { return CenterOp; }

      // Left and right components.  These can be NULL if we are at the edge of the system
      ComponentType const& Left() const { return LeftTop; }
      ComponentType const& Right() const { return RightTop; }

      ComponentType& Left() { return LeftTop; }
      ComponentType& Right() { return RightTop; }

      // Returns the boundary basis.  Normally, this will be one dimensional
      // and the RightBasisVacuum will be zero quantum number.  If the Center
      // matrix is rotated to the boundary, then this gives the center matrix basis
      MatrixBasisType LeftVacuumBasis() const;
      MatrixBasisType RightVacuumBasis() const;

      MatrixBasisType Basis1() const { return LeftVacuumBasis(); }
      MatrixBasisType Basis2() const { return RightVacuumBasis(); }

      CenterWavefunction& operator+=(CenterWavefunction const& x);
      CenterWavefunction& operator-=(CenterWavefunction const& x);

      // returns the component at n sites from the left-hand edge.
      // LookupLeft(LeftSize()-1) returns Left()
      ComponentType LookupLeft(int n) const 
      { RANGE_CHECK(n, 0, int(LeftStack.size())); 
      return std::size_t(n) == LeftStack.size() ? LeftTop : *LeftStack[n].load(); }

      void SetLeft(int n, ComponentType const& c)
      { 
         RANGE_CHECK(n, 0, int(LeftStack.size())); 
         if (std::size_t(n) == LeftStack.size())
            LeftTop = c;
         else LeftStack[n] = pvalue_ptr<ComponentType>(new ComponentType(c)); 
      }

#if 0
      ComponentType& LookupLeft(int n) 
      { RANGE_CHECK(n, 0, int(LeftStack.size())); 
      return std::size_t(n) == LeftStack.size() ? LeftTop : LeftStack[n]; }
#endif

      // returns the component at n sites from the right-hand edge.
      // LookupRight(RightSize()-1) returns Right()
      ComponentType LookupRight(int n) const 
      { RANGE_CHECK(n, 0, int(RightStack.size())); 
      return std::size_t(n) == RightStack.size() ? RightTop : *RightStack[n].load(); }

      void SetRight(int n, ComponentType const& c)
      { 
         RANGE_CHECK(n, 0, int(RightStack.size())); 
         if (std::size_t(n) == RightStack.size())
            RightTop = c;
         else RightStack[n] = pvalue_ptr<ComponentType>(new ComponentType(c)); 
      }

#if 0
      ComponentType& LookupRight(int n)
      { RANGE_CHECK(n, 0, int(RightStack.size())); 
      return std::size_t(n) == RightStack.size() ? RightTop : RightStack[n]; }
#endif

      ComponentType Lookup(int n) const
      {
	 if (n < this->LeftSize()) return LookupLeft(n);
	 else return LookupRight(this->size()-n-1);
      }

      ComponentType LookupLinear(int n) const
      {
	 if (n < this->LeftSize()-1) return this->LookupLeft(n);
         else if (n == this->LeftSize()-1) return prod(this->Left(), this->Center());
	 else return this->LookupRight(this->size()-n-1);
      }

      // Rotates the Center matrix to the left,
      // by PushRight(prod(Left, Center)); PopLeft();
      // and then calculating the (equivalent of) the
      // singular value decomposition on the new Right matrix
      // to get the new Center matrix.
      void RotateLeft();

      // the reverse of RotateLeft
      void RotateRight();

      // Rotates to the left, and expands the new left basis.
      // Equivalent to PushRight(prod(Left, Center)); PopLeft(); Center = ExpandBasis1(Right);
      void RotateLeftExpand();
      void RotateRightExpand();

      void RotateLeftTruncate(int NumStates);
      void RotateRightTruncate(int NumStates);

      // Rotates left until this->LeftSize() == 0
      void RotateToNormalForm();

      void PushLeft(ComponentType const& L);
      void PushRight(ComponentType const& R);

      void PopLeft();
      void PopRight();

      CenterWavefunction& normalize();
      CenterWavefunction& Normalize();

      LinearWavefunction AsLinearWavefunction() const;

      AttributeList const& Attributes() const { return Attr; }
      AttributeList& Attributes() { return Attr; }

      AttributeList& AttributesMutable() const { return Attr; }

   private:
      std::vector<pvalue_handle<ComponentType> > LeftStack, RightStack;
      ComponentType LeftTop, RightTop;
      OperatorType CenterOp;
      mutable AttributeList Attr;

   friend PStream::opstream& operator<<(PStream::opstream& out, CenterWavefunction const& Psi);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, CenterWavefunction& Psi);
};

inline
CenterWavefunction&
operator*=(CenterWavefunction& Psi, double x)
{
   Psi.Center() *= x;
   return Psi;
}

inline
CenterWavefunction&
operator*=(CenterWavefunction& Psi, std::complex<double> x)
{
   Psi.Center() *= x;
   return Psi;
}

inline
CenterWavefunction operator*(double x, CenterWavefunction const& Psi)
{
   CenterWavefunction Result(Psi);
   Result *= x;
   return Result;
}

inline
CenterWavefunction operator*(std::complex<double> x, CenterWavefunction const& Psi)
{
   CenterWavefunction Result(Psi);
   Result *= x;
   return Result;
}

inline
CenterWavefunction operator*(CenterWavefunction const& Psi, double x)
{
   CenterWavefunction Result(Psi);
   Result *= x;
   return Result;
}

inline
CenterWavefunction operator*(CenterWavefunction const& Psi, std::complex<double> x)
{
   CenterWavefunction Result(Psi);
   Result *= x;
   return Result;
}

inline
PStream::opstream& operator<<(PStream::opstream& out, CenterWavefunction const& MP)
{
   return out << MP.LeftStack << MP.RightStack << MP.LeftTop << MP.RightTop << MP.CenterOp << MP.Attr;
}

inline
PStream::ipstream& operator>>(PStream::ipstream& in, CenterWavefunction& MP)
{
   return in >> MP.LeftStack >> MP.RightStack >> MP.LeftTop >> MP.RightTop >> MP.CenterOp >> MP.Attr;
}

inline
std::ostream& operator<<(std::ostream& out, CenterWavefunction const& M)
{
   out << "Left matrices:\n";
   for (int i = 0; i < M.LeftSize(); ++i)
   {
      out << M.LookupLeft(i) << '\n';
   }
   out << "Center matrix:\n" << M.Center() << '\n';
   out << "Right matrices:\n";
   for (int i = M.RightSize(); i > 0; --i)
   {
      out << M.LookupRight(i-1) << '\n';
   }
   return out;
}

//
// norm_2
//

inline
double norm_2_sq(CenterWavefunction const& A)
{
   return trace(prod(A.Center(), adjoint(A.Center()), A.Center().TransformsAs())).real();
}

inline
double norm_2(CenterWavefunction const& A)
{
   return std::sqrt(norm_2_sq(A));
}

inline
double norm_frob_sq(CenterWavefunction const& A)
{
   return trace(prod(A.Center(), adjoint(A.Center()), A.Center().TransformsAs())).real();
}

inline
double norm_frob(CenterWavefunction const& A)
{
   return std::sqrt(norm_frob_sq(A));
}

inline
std::complex<double> overlap(CenterWavefunction const& A, CenterWavefunction const& B)
{
   return overlap(A.AsLinearWavefunction(), B.AsLinearWavefunction());
}

inline
std::complex<double>
expectation(CenterWavefunction const& Psi1, 
            MPOperator const& M, 
            CenterWavefunction const& Psi2)
{
   return expectation(Psi1.AsLinearWavefunction(), M, Psi2.AsLinearWavefunction());
}

inline
std::complex<double>
expectation(CenterWavefunction const& Psi1, 
            SplitOperator const& M, 
            CenterWavefunction const& Psi2)
{
   return expectation(Psi1.AsLinearWavefunction(), M.AsMPOperator(), Psi2.AsLinearWavefunction());
}

#endif
