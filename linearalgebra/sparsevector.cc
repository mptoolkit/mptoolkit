// -*- C++ -*- $Id$

#include "dataops.h"
#include <cmath>
#include <boost/type_traits.hpp>

namespace LinearAlgebra
{

template <class T>
SparseVector<T>::SparseVector() 
  : Index_begin(NULL), Index_end(NULL), Index_end_of_storage(NULL),
    Data_begin(NULL), Data_end(NULL), Data_end_of_storage(NULL), ZeroValue()
{
}

template <class T>
inline
SparseVector<T>::~SparseVector() 
{
   this->destroy();
}

template <class T>
void SparseVector<T>::destroy()
{
   for (T* x = Data_begin; x != Data_end; ++x)
   {
      x->~T();
   }
   delete[] static_cast<char*>(static_cast<void*>(Data_begin));
   delete[] Index_begin;
}

template <class T>
SparseVector<T>::SparseVector(const SparseVector<T>& V) : ZeroValue()
{
   int VSize = V.size();

   if (VSize != 0)
   {
      Data_begin = static_cast<T*>(static_cast<void*>(new char*[sizeof(T) * VSize]));
      //      Data_begin = new T[VSize];
      Data_end = Data_end_of_storage = Data_begin + VSize;

      Index_begin = new int[VSize];
      Index_end = Index_end_of_storage = Index_begin + VSize;

      std::uninitialized_copy(V.Data_begin, V.Data_end, Data_begin);
      //      std::copy(V.Data_begin, V.Data_end, Data_begin);
      std::copy(V.Index_begin, V.Index_end, Index_begin);
   }
   else
   {
      Data_begin = Data_end = Data_end_of_storage = NULL;
      Index_begin = Index_end = Index_end_of_storage = NULL;
   }
}

template <class T>
SparseVector<T>& SparseVector<T>::operator=(const SparseVector<T>& V)
{ 
   // we can't handle self assignment automatically
   if (&V == this) return *this;

   this->destroy();

   int VSize = V.size();

   if (VSize != 0)
   {
      Data_begin = static_cast<T*>(static_cast<void*>(new char*[sizeof(T) * VSize]));
      Data_end = Data_end_of_storage = Data_begin + VSize;

      Index_begin = new int[VSize];
      Index_end = Index_end_of_storage = Index_begin + VSize;

      std::uninitialized_copy(V.Data_begin, V.Data_end, Data_begin);
      std::copy(V.Index_begin, V.Index_end, Index_begin);
   }
   else
   {
      Data_begin = Data_end = Data_end_of_storage = NULL;
      Index_begin = Index_end = Index_end_of_storage = NULL;
   }

   return *this; 
} 

template <class T>
int SparseVector<T>::raw_index(int i) const
{
   const_index_iterator I = std::find(index_begin(), index_end(), i);
   return I == index_end() ? -1 : I - index_begin();
}

template <class T>
bool SparseVector<T>::has_element(int i) const 
{ 
   return std::find(index_begin(), index_end(), i) != index_end();
}

// template <class T>
// typename dot_traits<T>::value_type
// magnitude(const SparseVector<T>& S) 
// { 
//    return sqrt(magnitude_squared(S)); 
// }

template <class T>
void
SparseVector<T>::clear()
{
   this->destroy();
   Data_begin = Data_end = Data_end_of_storage = NULL;
   Index_begin = Index_end = Index_end_of_storage = NULL;
}

template <class T>
void
SparseVector<T>::negate()
{
   for (T* D = Data_begin; D != Data_end; ++D) *D = -*D;
}

template <class T>
SparseVector<T>& 
SparseVector<T>::operator+=(const SparseVector<T>& V)
{
   reserve(size() + V.size()); // reserve the worst case memory

   const_iterator End = V.end();
   for (const_iterator I = V.begin(); I != End; ++I)
      add(I.index(), I.data());

   return *this;
}

template <class T>
SparseVector<T>& 
SparseVector<T>::operator-=(const SparseVector<T>& V)
{
   reserve(size() + V.size()); // reserve the worst case memory

   const_iterator End = V.end();
   for (const_iterator I = V.begin(); I != End; ++I)
      subtract(I.index(), I.data());

   return *this;
}

template <class T>
template <class Scalar>
SparseVector<T>& 
SparseVector<T>::add_vector(const SparseVector<T>& V, Scalar S)
{
   reserve(size() + V.size()); // reserve the worst case memory

   const_iterator End = V.end();
   for (const_iterator I = V.begin(); I != End; ++I)
      add(I.index(), S * I.data());
   
   return *this;
}

template <class T>
template <class Scalar>
SparseVector<T>& 
SparseVector<T>::operator*=(Scalar V)
{
   for (T* D = Data_begin; D != Data_end; ++D) *D *= V;
   return *this;
}

template <class T>
void
SparseVector<T>::reallocate_storage(int NewSize)
{
   this->destroy();

   if (NewSize > 0)
   {
      Data_begin = static_cast<T*>(static_cast<void*>(new char*[sizeof(T) * NewSize]));
      Data_end = Data_end_of_storage = Data_begin + NewSize;
      std::uninitialized_fill(Data_begin, Data_end, T());

      Index_begin = new int[NewSize];
      Index_end = Index_end_of_storage = Index_begin + NewSize;
   }
   else
   {
      Data_begin = Data_end = Data_end_of_storage = NULL;
      Index_begin = Index_end = Index_end_of_storage = NULL;
   }
}

template <class T>
void
SparseVector<T>::increase_capacity()
{
   increase_capacity(std::max<int>(int(std::ceil(capacity() * 1.5)), MinCapacity));
}

template <class T>
void
SparseVector<T>::increase_capacity(int NewCapacity)
{
   int OldCapacity = capacity();
   if (OldCapacity >= NewCapacity) return; // don't shrink the allocation

   int Size = this->size();

   //   T* new_Data_begin = new T[NewCapacity];
   T* new_Data_begin = static_cast<T*>(static_cast<void*>(new char[sizeof(T) * NewCapacity]));
   int* new_Index_begin = new int[NewCapacity];

   //   std::copy(Data_begin, Data_end, new_Data_begin);   
   std::uninitialized_copy(Data_begin, Data_end, new_Data_begin);   
   //   if (!boost::has_trivial_constructor<T>::value)
   //      std::uninitialized_fill(new_Data_begin+Size, new_Data_begin+NewCapacity, T());

   std::copy(Index_begin, Index_end, new_Index_begin);   

   this->destroy();

   Data_begin = new_Data_begin;
   Data_end_of_storage = Data_begin + NewCapacity;
   Data_end = Data_begin + Size;

   Index_begin = new_Index_begin;
   Index_end_of_storage = Index_begin + NewCapacity;
   Index_end = Index_begin + Size;
}

template <class T>
template <typename U>
int
SparseVector<T>::increase_capacity_and_append_new(int i, U const& S)
{
   int NewCapacity = std::max<int>(int(std::ceil(capacity() * 1.5)), MinCapacity);

   int Size = this->size();

   //   T* new_Data_begin = new T[NewCapacity];
   T* new_Data_begin = static_cast<T*>(static_cast<void*>(new char[sizeof(T) * NewCapacity]));
   int* new_Index_begin = new int[NewCapacity];

   //   std::copy(Data_begin, Data_end, new_Data_begin);   
   std::uninitialized_copy(Data_begin, Data_end, new_Data_begin);   
   //   if (!boost::has_trivial_constructor<T>::value)
   //      std::uninitialized_fill(new_Data_begin+Size, new_Data_begin+NewCapacity, T());

   std::copy(Index_begin, Index_end, new_Index_begin);

   // Do the assignment before delete[] Data_begin, in case S is a reference
   // to the existing array
   new_Index_begin[Size] = i;
   new (new_Data_begin+Size) T(S);
   //   new_Data_begin[Size] = S;
   ++Size;

   this->destroy();

   Data_begin = new_Data_begin;
   Data_end_of_storage = Data_begin + NewCapacity;
   Data_end = Data_begin + Size;

   Index_begin = new_Index_begin;
   Index_end_of_storage = Index_begin + NewCapacity;
   Index_end = Index_begin + Size;

   return Size-1;
}

} // namespace LinearAlgebra

