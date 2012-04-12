/* -*- C++ -*- $Id$

*** DEPRECIATED ***

  Defines MPOpComponent, which represents an element of
  a matrix product operator.

  TODO:
    Fix remaining functions that still access the members of MPOpComponent directly,
    then make the data members private.

    It would be better style to remove data() completely.

    I prefer the naming 'local basis' over the original 'site basis'.  I intend to eventually
    add functions with 'site' -> 'local' everywhere and eventually remove the 'site' variants.

*/

#if !defined(MPOPCOMPONENT_H_SDHCKJHKJLRHGIURHYULUHR)
#define MPOPCOMPONENT_H_SDHCKJHKJLRHGIURHYULUHR

#include "tensor/tensorproduct.h"

using namespace Tensor;

typedef IrredTensor<std::complex<double> > SimpleOperator;

typedef IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, 
                            VectorBasis, 
                            VectorBasis> MatrixOperator;



typedef IrredTensor<SimpleOperator, BasisList, BasisList> CompoundOperator;

template <typename Component>
struct MPOperatorComponent;

template <typename T>
PStream::opstream& operator<<(PStream::opstream& out, MPOperatorComponent<T> const& Op);

template <typename T>
PStream::ipstream& operator>>(PStream::ipstream& in, MPOperatorComponent<T>& Op);

template <typename Component>
struct MPOperatorComponent
{
   typedef Component OperatorType;
   typedef typename OperatorType::basis1_type basis1_type;
   typedef typename OperatorType::basis2_type basis2_type;

   // TODO: this definition should be conditional on basis1_type == basis2_type
   typedef basis1_type BasisType;

   typedef IrredTensor<OperatorType, BasisList, BasisList> mapped_type;
 
   typedef std::map<QuantumNumber, mapped_type> DataType;

   typedef typename DataType::value_type value_type;
   typedef typename DataType::iterator iterator;
   typedef typename DataType::const_iterator const_iterator;

   MPOperatorComponent() {}

   MPOperatorComponent(BasisList const& SBasis, basis1_type const& V1, basis2_type const& V2);

   SymmetryList const& GetSymmetryList() const { return SBasis_.GetSymmetryList(); }

   BasisList const& SiteBasis() const { return SBasis_; }
   BasisList const& LocalBasis() const { return SBasis_; }

   basis1_type const& Basis1() const { return Basis1_; }
   basis2_type const& Basis2() const { return Basis2_; }

   mapped_type& operator[](QuantumNumber const& q);

   mapped_type operator[](QuantumNumber const& q) const;

   iterator begin() { return Data_.begin(); }
   iterator end() { return Data_.end(); }
 
   const_iterator begin() const { return Data_.begin(); }
   const_iterator end() const { return Data_.end(); }

   iterator find(QuantumNumber const& q) { return Data_.find(q); }
   const_iterator find(QuantumNumber const& q) const { return Data_.find(q); }

   bool is_null() const { return SBasis_.size() == 0; }

   // sets the matrix elements this(a,b)(i,j) = x(a,b)
   void set_operator(int i, int j, SimpleOperator const& x);

   MPOperatorComponent& operator*=(double x);

   MPOperatorComponent& operator*=(std::complex<double> x);

   static MPOperatorComponent ConstructFullBasis1(BasisList const& SiteBasis, 
                                                  basis2_type const& Basis2);

   static MPOperatorComponent ConstructFullBasis2(basis1_type const& Basis1, 
                                                  BasisList const& SiteBasis);

   DataType& data() { return Data_; }
   DataType const& data() const { return Data_; }

   void check_structure() const;
   void debug_check_structure() const;

   BasisList SBasis_;
   basis1_type Basis1_;
   basis2_type Basis2_;

   DataType Data_;

   friend PStream::opstream& operator<< <>(PStream::opstream& out, MPOperatorComponent const& Op);
   friend PStream::ipstream& operator>> <>(PStream::ipstream& in, MPOperatorComponent& Op);
};

// backwards-compatibility typedef
typedef MPOperatorComponent<SimpleOperator> MPOpComponent;

typedef MPOperatorComponent<SimpleOperator> MPSimpleOpComponent;
typedef MPOperatorComponent<MatrixOperator> MPMatrixOpComponent;

template <typename T>
inline
void MPOperatorComponent<T>::debug_check_structure() const
{
#if !defined(NDEBUG)
   this->check_structure();
#endif
}

inline
MPOpComponent ConstructIdentity(BasisList const& B)
{
   QuantumNumber Ident(B.GetSymmetryList());
   BasisList Vac = make_vacuum_basis(B.GetSymmetryList());
   SimpleOperator IdentOperator = SimpleOperator::make_identity(Vac);

   MPOpComponent Result(B, Vac, Vac);
   for (std::size_t i = 0; i < B.size(); ++i)
   {
      set_element(Result[Ident].data(), i,i, IdentOperator);
   }
   return Result;
}

bool IsProportionalIdentity(MPOpComponent const& Op);

std::complex<double>
IdentityScale(MPOpComponent const& Op);

template <typename T>
std::ostream& operator<<(std::ostream& out, MPOperatorComponent<T> const& Op);

template <typename T>
inline
MPOperatorComponent<T>::MPOperatorComponent(BasisList const& SBasis, 
                                            basis1_type const& V1, 
                                            basis2_type const& V2)
  : SBasis_(SBasis), Basis1_(V1), Basis2_(V2)
{
}

template <typename T>
inline
typename MPOperatorComponent<T>::mapped_type&
MPOperatorComponent<T>::operator[](QuantumNumber const& q)
{
   if (Data_.find(q) == Data_.end())
   {
      Data_[q] = mapped_type(this->SiteBasis(), q);
   }
   return Data_[q];
}

template <typename T>
inline
typename MPOperatorComponent<T>::mapped_type
MPOperatorComponent<T>::operator[](QuantumNumber const& q) const
{
   typename DataType::const_iterator I = Data_.find(q);
   if (I == Data_.end())
   {
      return mapped_type(this->SiteBasis(), q);
   }
   return I->second;
}

namespace LinearAlgebra
{

template <>
struct interface<MPOpComponent>
{
   typedef void type;
};

// hermitian conjugation

template <typename T>
struct Herm<MPOperatorComponent<T> >
{
   typedef HermitianProxy<MPOperatorComponent<T> > result_type;
   typedef MPOperatorComponent<T> const& argument_type;

   result_type operator()(argument_type x) const
   {
      return result_type(x);
   }
};

// scalar_prod

// does Result' = sum_s A[s] * herm(B[s])
template <typename T>
struct ScalarProd<MPOperatorComponent<T>, HermitianProxy<MPOperatorComponent<T> > >
{
   typedef typename MPOperatorComponent<T>::OperatorType result_type;
   typedef MPOperatorComponent<T> const& first_argument_type;
   typedef HermitianProxy<MPOperatorComponent<T> > const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y) const;
};

template <typename T>
struct ScalarProd<HermitianProxy<MPOperatorComponent<T> >, MPOperatorComponent<T> >
{
   typedef typename MPOperatorComponent<T>::OperatorType result_type;
   typedef HermitianProxy<MPOperatorComponent<T> > const& first_argument_type;
   typedef MPOperatorComponent<T> const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y) const;
};

} // namespace LinearAlgebra

#include "mpopcomponent.cc"

#endif
