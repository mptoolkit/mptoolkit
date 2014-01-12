// -*- C++ -*- $Id$

#include "generic_mpo.h"

std::vector<BasisList> 
GenericMPO::LocalBasis1List() const
{
   std::vector<BasisList> Result;
   Result.reserve(Data_.size());
   for (unsigned i = 0; i < Data_.size(); ++i)
   {
      Result.push_back(Data_[i].LocalBasis1());
   }
   return Result;
}

std::vector<BasisList> 
GenericMPO::LocalBasis2List() const
{
   std::vector<BasisList> Result;
   Result.reserve(Data_.size());
   for (unsigned i = 0; i < Data_.size(); ++i)
   {
      Result.push_back(Data_[i].LocalBasis2());
   }
   return Result;
}

std::ostream&
operator<<(std::ostream& out, GenericMPO const& op)
{
   out << "Operator has a unit cell of " << op.size() << " sites.\n";
   for (unsigned i = 0; i < op.size(); ++i)
   {
      out << "Site " << i << ": " << op[i] << '\n';
   }
   return out;
}

bool
GenericMPO::is_null() const
{
   for (unsigned i = 0; i < Data_.size(); ++i)
   {
      if (!Data_[i].is_null())
         return false;
   }
   return true;
}

PStream::opstream& operator<<(PStream::opstream& out, GenericMPO const& op)
{
   return out << op.Data_;
}

PStream::ipstream& operator>>(PStream::ipstream& in, GenericMPO& op)
{
   return in >> op.Data_;
}

#if 0
GenericMPO&
operator*=(GenericMPO& x, double a)
{
   x.front() *= a;
   return x;
}

GenericMPO&
operator*=(GenericMPO& x, std::complex<double> a)
{
   x.front() *= a;
   return x;
}

GenericMPO operator*(double a, GenericMPO const& x)
{
   GenericMPO Result(x);
   Result *= a;
   return Result;
}

GenericMPO operator*(GenericMPO const& x, double a)
{
   GenericMPO Result(x);
   Result *= a;
   return Result;
}

GenericMPO operator*(std::complex<double> a, GenericMPO const& x)
{
   GenericMPO Result(x);
   Result *= a;
   return Result;
}

GenericMPO operator*(GenericMPO const& x, std::complex<double> a)
{
   GenericMPO Result(x);
   Result *= a;
   return Result;
}
#endif

void zero_unused_elements(GenericMPO& Op)
{
   bool Done = false;
   while (!Done)
   {
      Done = true;
      std::set<int> NextKeep;
      GenericMPO::iterator I = Op.begin();

      std::set<int> RowsToKeep;
      for (unsigned i = 0; i < I->size1(); ++i)
         for (unsigned j = 0; j < I->size2(); ++j)
            if (I->iterate_at(i,j))
               {
                  RowsToKeep.insert(j);
               }
      ++I;

      while (I != Op.end())
      {
         OperatorComponent Current = *I;
         for (unsigned i = 0; i < Current.size1(); ++i)
         {
            for (unsigned j = 0; j < Current.size2(); ++j)
            {
               if (!Current.iterate_at(i,j))
                  continue;

               if (RowsToKeep.count(i))
                  NextKeep.insert(j);
               else
               {
                  zero_element(Current.data(), i,j);
                  Done = false;
               }

            }
         }
         *I = Current;

         RowsToKeep = NextKeep;
         NextKeep.clear();

         ++I;
      }

      // now work backwards
      --I;
      std::set<int> ColumnsToKeep;
      for (unsigned i = 0; i < I->size1(); ++i)
         for (unsigned j = 0; j < I->size2(); ++j)
            if (I->iterate_at(i,j))
               ColumnsToKeep.insert(i);

      while (I != Op.begin())
      {
         --I;

         OperatorComponent Current = *I;
         for (unsigned i = 0; i < Current.size1(); ++i)
         {
            for (unsigned j = 0; j < Current.size2(); ++j)
            {
               if (!Current.iterate_at(i,j))
                  continue;

               if (ColumnsToKeep.count(j))
                  NextKeep.insert(i);
               else
               {
                  zero_element(Current.data(), i,j);
                  Done = false;
               }

            }
         }
         *I = Current;

         ColumnsToKeep = NextKeep;
         NextKeep.clear();
      }

   } // while (!Done)
}

bool update_mask(OperatorComponent const& x, OperatorComponent const& y, 
                 std::vector<int> const& M1, std::vector<int>& Mask, std::vector<int> const& M2)
{
   typedef OperatorComponent::data_type OpMatrixType;
   typedef OpMatrixType::outer_value_type VecOfVecType;
   //   typedef VecOfVecType::value_type InnerType;
   typedef LinearAlgebra::MapVector<SimpleRedOperator> InnerType;

   bool Updated = false;
   for (const_iterator<OperatorComponent::data_type>::type I = iterate(x.data()); I; ++I)
   {
      if (!M1[I.index()])
         continue;

      for (const_inner_iterator<OperatorComponent::data_type>::type J = iterate(I); J; ++J)
      {
         //         if (M1[J.index1()] && Mask[J.index2()])
         if (Mask[J.index2()])
         {
            const_iterator<InnerType>::type K = iterate(y.data().vec()[J.index2()]);
            while (K && !M2[K.index()])
               ++K;

            if (!K)
            {
               Mask[J.index2()] = 0;
               Updated = true;
            }
         }
      }
   }

   return Updated;
}

bool update_mask(OperatorComponent const& x, std::vector<int>& M1, std::vector<int>& M2)
{
   typedef OperatorComponent::data_type OpMatrixType;
   typedef OpMatrixType::outer_value_type VecOfVecType;
   //   typedef VecOfVecType::value_type InnerType;
   typedef LinearAlgebra::MapVector<SimpleRedOperator> InnerType;

   bool Updated = false;
   for (const_iterator<OperatorComponent::data_type>::type I = iterate(x.data()); I; ++I)
   {
      if (!M1[I.index()])
         continue;

      bool Found = false;
      for (const_inner_iterator<OperatorComponent::data_type>::type J = iterate(I); J; ++J)
      {
         if (M2[J.index2()])
         {
            M2[J.index2()] = 2;  // tag as accessible
            Found = true;
         }
      }
      if (!Found)
      {
         M1[I.index()] = 0;
         Updated = true;
      }
   }

   for (unsigned j = 0; j < M2.size(); ++j)
   {
      if (M2[j] == 2)
      {
         M2[j] = 1;
      }
      else if (M2[j] == 1)
      {
         M2[j] = 0;  // not accessible
         Updated = true;
      }
   }

   return Updated;
}

// remove unused components from the shared basis between x and y
bool cull_boundary(OperatorComponent& x, OperatorComponent& y)
{
   std::set<int> ToKeep;
   for (const_iterator<OperatorComponent::data_type>::type I = iterate(x.data()); I; ++I)
   {
      for (const_inner_iterator<OperatorComponent::data_type>::type J = iterate(I); J; ++J)
      {
         if (nnz(y.data().vec()[J.index2()]) > 0)
            ToKeep.insert(J.index2());
      }
   }
   if (ToKeep.size() != x.Basis2().size())
   {
      x = project_columns(x, ToKeep);
      y = project_rows(y, ToKeep);
      return true;
   }
   return false;
}

void cull_unused_elements(GenericMPO& Op)
{
   // We need at least one bond to optimize
   if (Op.size() < 2)
      return;

   bool Done = false;
   while (!Done)
   {
      Done = true;
      GenericMPO::iterator I = Op.begin();
      GenericMPO::iterator J = I; ++J;

      while (J != Op.end())
      {
         if (cull_boundary(*I, *J))
            Done = false;

	 I=J;
	 ++J;
      }

      // now work backwards and cull columns
      while (I != Op.begin())
      {
	 J=I;
	 --I;
       
         if (cull_boundary(*I, *J))
            Done = false;
      }
   } // while (!Done)
}

void mask_unused_elements(GenericMPO const& Op, std::vector<std::vector<int> >& Mask)
{
   if (Op.size() < 1)
      return;

   GenericMPO::const_iterator I = Op.begin();
   std::vector<std::vector<int> >::iterator M1 = Mask.begin();
   std::vector<std::vector<int> >::iterator M2 = Mask.begin()+1;

   bool Done = false;
   while (!Done)
   {
      Done = true;
      while (I != Op.end())
      {
         if (update_mask(*I, *M1, *M2))
            Done = false;

         ++I;
         ++M1;
         ++M2;
      }

      if (Done)
         break;

      while (I != Op.begin())
      {
         --I;
         --M1;
         --M2;

         if (update_mask(*I, *M1, *M2))
            Done = false;
      }
   }
}

void initialize_mask(GenericMPO const& Op, std::vector<std::vector<int> >& Mask)
{
   std::vector<std::vector<int> >(Op.size()+1).swap(Mask);
   for (unsigned i = 0; i < Op.size(); ++i)
   {
      Mask[i] = std::vector<int>(Op[i].Basis1().size(), true);
   }
   Mask.back() = std::vector<int>(Op.Basis2().size(), true);
}

SimpleOperator make_projector_onto(BasisList const& Basis, std::set<int> const& Onto)
{
   BasisList ProjectedBasis(Basis.GetSymmetryList());
   for (std::set<int>::const_iterator I = Onto.begin(); I != Onto.end(); ++I)
   {
      ProjectedBasis.push_back(Basis[*I]);
   }

   SimpleOperator Result(ProjectedBasis, Basis);
   int j = 0;
   for (std::set<int>::const_iterator I = Onto.begin(); I != Onto.end(); ++I)
   {
      Result(j++, *I) = 1.0;
   }
   return Result;
}

// classification

GenericMPOClassification::GenericMPOClassification()
   : Factor_(0.0), Product_(false), Unitary_(false),
     Identity_(false), PropUnitary_(false), PropIdentity_(false), Null_(false)
{
}

bool GenericMPOClassification::is_null() const
{
   return Null_;
}

bool GenericMPOClassification::is_product() const
{
   return Product_;
}

bool GenericMPOClassification::is_unitary() const
{
   return Unitary_;
}

bool GenericMPOClassification::is_prop_unitary() const
{
   return PropUnitary_;
}

bool GenericMPOClassification::is_prop_identity() const
{
   return PropIdentity_;
}

bool GenericMPOClassification::is_identity() const
{
   return Identity_;
}

bool GenericMPOClassification::is_unclassified() const
{
   return !Product_ && !Null_;
}

std::complex<double> GenericMPOClassification::factor() const
{
   return Factor_;
}

std::ostream& operator<<(std::ostream& out, GenericMPOClassification const& Class)
{
   out << "null: " << Class.is_null() << '\n';
   out << "product: " << Class.is_product() << '\n';
   out << "unitary: " << Class.is_unitary() << '\n';
   out << "prop_unitary: " << Class.is_prop_unitary() << '\n';
   out << "prop_identity: " << Class.is_prop_identity() << '\n';
   out << "complex_identity: " << Class.is_complex_identity() << '\n';
   out << "identity: " << Class.is_identity() << '\n';
   out << "factor: " << Class.factor() << '\n';
   return out;
}

std::complex<double> PropIdent(SimpleOperator const& X)
{
   DEBUG_PRECONDITION_EQUAL(X.Basis1(), X.Basis2());
   SimpleOperator Ident = SimpleOperator::make_identity(X.Basis1());
   std::complex<double> x = inner_prod(Ident, X) / norm_frob_sq(Ident);
   //   TRACE(x);
   //TRACE(X-x*Ident);
   if (norm_frob(X-x*Ident) > std::numeric_limits<double>::epsilon()*10)
      x = 0.0;
   return x;
}

GenericMPOClassification classify(GenericMPO const& Op)
{
   GenericMPOClassification Result;

   // Early return if the operator is null
   if (Op.is_null())
   {
      Result.Null_ = true;
      return Result;
   }

   bool IsPropIdentity = true;  // true if the operator is proportional to identity
   bool IsPropUnitary = true;   // true if the operator is proportional to a unitary operator
   bool IsUnitUnitary = true;   // true if the operator is unitary (ie. proportional, with factor 1.0)
   std::complex<double> Factor  = 1.0;

   for (unsigned i = 0; i < Op.size(); ++i)
   {
      // firstly, check to see if it is 1x1
      if (Op[i].Basis1().size() != 1 || Op[i].Basis2().size() != 1)
         return Result;  // default constructed return is unclassified

      if (IsPropUnitary)
      {
         SimpleRedOperator X = Op[i](0,0);

         if (IsPropIdentity)
         {
            if (X.Basis1() != X.Basis2() || !is_pure_scalar(X))
               IsPropIdentity = false;
            else
            {
               std::complex<double> x = PropIdent(X.scalar());
               if (x == 0.0)
                  IsPropIdentity = false;
               else
                  Factor *= x;
            }
         }

         if (!IsPropIdentity)
         {
            // is it unitary?
            std::complex<double> x = PropIdent(scalar_prod(X, herm(X)));
            std::complex<double> y = PropIdent(scalar_prod(herm(X), X));

            if (x == 0.0 || y == 0.0)
            {
               IsPropUnitary = false;
            }
            else
            {
               if (norm_frob(x - 1.0) + norm_frob(y - 1.0) > std::numeric_limits<double>::epsilon()*200)
                  IsUnitUnitary = false;

            }
         }

      }
   }

   Result.Product_ = true;
   if (IsPropUnitary)
   {
      Result.PropUnitary_ = true;

      if (IsPropIdentity)
      {
         Result.PropIdentity_ = true;
         //TRACE(Factor);
         Result.Identity_ = norm_frob(Factor - std::complex<double>(1.0, 0)) < std::numeric_limits<double>::epsilon()*100;

         // if we claim to be an identity operator, we might as well make it exact
         if (Result.Identity_)
            Factor = 1.0;

         Result.Unitary_ = norm_frob(norm_frob(Factor) - 1.0) < std::numeric_limits<double>::epsilon()*100;
         Result.Factor_ = Factor;
      }
      else
      {
         Result.Unitary_ = IsUnitUnitary;
      }
   }

   return Result;
}

std::vector<BasisList>
ExtractLocalBasis1(GenericMPO const& Op)
{
   std::vector<BasisList> Result(Op.size());
   for (unsigned i = 0; i < Op.size(); ++i)
   {
      Result[i] = Op[i].LocalBasis1();
   }
   return Result;
}
