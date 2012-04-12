// -*- C++ -*- $Id$

#include "dataops.h"
#include "matrix_temp.h"

namespace LinearAlgebra
{

template <class Scalar>
void 
MatrixRef<Scalar, ConcreteClass>::resize(size_type NewRows, size_type NewCols, bool PreserveValues)
{
   // Early return in the case that the data block is not shared and we don't need a resize.
   if (NewRows == this->size1() && NewCols == this->size2() && this->stride2() == 1 && !Block.is_shared()) return;

   // this->cow(); // FIXME: pretty sure this is not needed???

   block_type NewData(NewRows*NewCols, NoHeader(), MatrixDimensions(NewRows, NewCols));
   if (PreserveValues)
   {
      PANIC("MatrixRef::resize() with PreserveValues=true is not yet implemented!");
      //      ops::fast_copy(Block.get(), Block.get()+std::min(Size, NewSize), NewData.get());
      //      ops::fast_fill(NewData.get() + std::min(Size, NewSize), NewData.get() + NewSize, Scalar());
   }
   Block = NewData;
}

#if defined(USE_PSTREAM)

// TODO: these are particularly brain-dead

template <int Format, typename T>
PStream::opstreambuf<Format>& operator<<(PStream::opstreambuf<Format>& out, Matrix<T> const& M)
{
   typename PStream::opstreambuf<Format>::size_type Rows = M.size1(), Cols = M.size2();
   out << Rows << Cols;

   typedef typename Matrix<T>::const_iterator1 OuterIter;
   typedef typename OuterIter::const_iterator InnerIter;

   for (OuterIter i = M.begin1(); i != M.end1(); ++i)
   {
      for (InnerIter j = i.begin(); j != i.end(); ++j)
      {
	 out << *j;
      }
   }
   return out;
}

template <int Format, typename T>
PStream::ipstreambuf<Format>& operator>>(PStream::ipstreambuf<Format>& in, Matrix<T>& M)
{
   typename PStream::ipstreambuf<Format>::size_type Rows, Cols;
   in >> Rows >> Cols;
   M.resize(Rows, Cols, false);

   typedef typename Matrix<T>::iterator1 OuterIter;
   typedef typename OuterIter::iterator InnerIter;

   for (OuterIter i = M.begin1(); i != M.end1(); ++i)
   {
      for (InnerIter j = i.begin(); j != i.end(); ++j)
      {
	 in >> *j;
      }
   }
   return in;
}
#endif

template <class Scalar>
Matrix<Scalar>::Matrix(size_type Rows, size_type Cols)
  : MatrixRef<Scalar, ConcreteClass>(Rows, Cols)
{
}

template <class Scalar>
Matrix<Scalar>::Matrix(size_type Rows, size_type Cols, Scalar const& fill_)
  : MatrixRef<Scalar, ConcreteClass>(Rows, Cols, fill_)
{
}

template <class Scalar>
template <typename Iter>
Matrix<Scalar>::Matrix(size_type Rows, size_type Cols, Iter first, Iter last)
  : MatrixRef<Scalar, ConcreteClass>(Rows, Cols)
{
   // FIXME: this really should call an overload of assign()
   PRECONDITION(last-first == Rows*Cols);
   ops::fast_copy(first, last, this->begin());
}

template <class Scalar>
template <typename S, typename D>
Matrix<Scalar>::Matrix(GenericMatrix<S, D> const& V)
  : MatrixRef<Scalar, ConcreteClass>(V.size1(), V.size2())
{
   this->assign(V);
}

template <class Scalar>
template <class S2, class D2>
Matrix<Scalar>&
Matrix<Scalar>::operator=(GenericMatrix<S2 const, D2> const& V)
{
   TempMatrix<S2> Temp(V);
   this->resize(V.size1(), V.size2(), false);
   this->assign(V.as_derived());
   return *this;
}

template <class Scalar>
Matrix<Scalar>&
Matrix<Scalar>::operator=(Matrix<Scalar> const& V)
{
   this->block() = V.block().copy();
   return *this;
}

} // namespace LinearAlgebra
