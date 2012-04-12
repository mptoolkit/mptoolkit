// -*- C++ -*- $Id$
//
// pack & unpack functions for MatrixOperator.  These convert to/from
// an inner-product conserving vector representation.
// Preserving the inner product means we need an SU(2) coefficient.
//

#if !defined(PACKUNPACK_H_ASCHUIY243Y478Y648E89)
#define PACKUNPACK_H_ASCHUIY243Y478Y648E89

#include "state_component.h"  // for MatrixOperator definition
#include "linearalgebra/vector.h"
#include "linearalgebra/matrix.h"
#include <complex>

class PackMatrixOperator
{
   public:
      typedef std::complex<double> value_type;

      PackMatrixOperator(VectorBasis const& Basis1, 
                         VectorBasis const& Basis2, 
                         QuantumNumbers::QuantumNumber const& q);

      PackMatrixOperator(MatrixOperator const& m);

      value_type* pack(MatrixOperator const& m, value_type* Iter) const;

      MatrixOperator unpack(value_type const* Iter) const;
   
      std::size_t size() const { return Size_; }
   
   private:
      void Initialize(VectorBasis const& Basis1, 
                      VectorBasis const& Basis2, 
                      QuantumNumbers::QuantumNumber const& q_);

      // OffsetRecType contains the information for an individual dense block
      struct OffsetRecType
      {
         int r, c;          // row and column of the MatrixOperator
         int Offset;        // offset into the packed array of this section
         unsigned Size;     // size of this block
         OffsetRecType() {}
         OffsetRecType(int r_, int c_, int Offset_, int Size_)
            : r(r_), c(c_), Offset(Offset_), Size(Size_) {}
      };

      typedef  std::vector<OffsetRecType> OffsetArrayType;

      VectorBasis B1_, B2_;
      QuantumNumbers::QuantumNumber q_;
      int Size_;  // linear size
      std::vector<OffsetRecType> OffsetArray_;   // array of blocks
      LinearAlgebra::Matrix<int> OffsetMatrix_;  // offset of each valid block
};

// construct the full matrix representation of some superoperator
// given by the functor F : MatrixOperator -> MatrixOperator
template <typename F>
LinearAlgebra::Matrix<std::complex<double> >
ConstructSuperOperator(F f, MatrixOperator const& Init);

template <typename F>
LinearAlgebra::Matrix<std::complex<double> >
ConstructSuperOperator(F f, MatrixOperator const& Init)
{
   PackMatrixOperator Pack(Init);

   std::size_t Size = Pack.size();

   TRACE(Size)(Size*Size);

   LinearAlgebra::Matrix<std::complex<double> > Out(Size, Size);

   std::vector<std::complex<double> > L(Size), R(Size);
   for (unsigned i = 0; i < Size; ++i)
   {
      if (i > 0)
         R[i-1] = 0.0;

      R[i] = 1.0;

      MatrixOperator M = Pack.unpack(&R[0]);
      M = f(M);
      Pack.pack(M, &L[0]);
      for (unsigned j = 0; j < Size; ++j)
      {
         Out(j,i) = L[j];
      }
   }

   return Out;
}

template <typename F>
LinearAlgebra::SparseMatrix<std::complex<double> >
ConstructSuperOperatorSparse(F f, MatrixOperator const& Init)
{
   PackMatrixOperator Pack(Init);

   std::size_t Size = Pack.size();

   TRACE(Size)(Size*Size);

   LinearAlgebra::SparseMatrix<std::complex<double> > Out(Size, Size);

   std::vector<std::complex<double> > L(Size), R(Size);
   for (unsigned i = 0; i < Size; ++i)
   {
      if (i > 0)
         R[i-1] = 0.0;

      R[i] = 1.0;

      MatrixOperator M = Pack.unpack(&R[0]);
      M = f(M);
      Pack.pack(M, &L[0]);
      for (unsigned j = 0; j < Size; ++j)
      {
	 if (norm_frob(Out(j,i)) != 0)
	    Out(j,i) = L[j];
      }
   }

   return Out;
}

// write the operator to a file, column major format
template <typename F>
std::streambuf*
ConstructSuperOperatorBinaryFile(F f, MatrixOperator const& Init, std::streambuf* Out)
{
   PackMatrixOperator Pack(Init);

   std::size_t Size = Pack.size();

   TRACE(Size)(Size*Size);


   std::vector<std::complex<double> > L(Size), R(Size);
   for (unsigned i = 0; i < Size; ++i)
   {
      if (i > 0)
         R[i-1] = 0.0;

      R[i] = 1.0;

      MatrixOperator M = Pack.unpack(&R[0]);
      M = f(M);
      Pack.pack(M, &L[0]);
      for (unsigned j = 0; j < Size; ++j)
      {
	 double r = L[j].real();
	 double c = L[j].real();
	 Out->sputn(reinterpret_cast<char const*>(&r), sizeof(r));
	 Out->sputn(reinterpret_cast<char const*>(&c), sizeof(c));
      }
   }

   return Out;
}

#endif
