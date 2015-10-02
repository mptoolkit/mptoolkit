// -*- C++ -*- $Id$

#include "wigner_eckart.h"

namespace Tensor
{

inline
int 
WignerEckartBasis<VectorBasis>::Map(int NonAbelianSubspace, int Projection) const
{
   std::map<int, int>::const_iterator I = Mapping[NonAbelianSubspace].find(Projection);
   if (I == Mapping[NonAbelianSubspace].end())
      return -1;
   // else
   return I->second;
}

template <typename T, typename B1T, typename B2T>
IrredTensor<T, B1T, B2T, DefaultStructure>
wigner_eckart(IrredTensor<T, B1T, B2T, DefaultStructure> const& x, 
              Projection const& p, 
              WignerEckartBasis<B1T> const& b1, 
              WignerEckartBasis<B1T> const& b2)
{
   DEBUG_CHECK_EQUAL(b1.AbelianBasis().GetSymmetryList(),
                     b2.AbelianBasis().GetSymmetryList());
   DEBUG_CHECK_EQUAL(b1.NonAbelianBasis().GetSymmetryList(),
                     b2.NonAbelianBasis().GetSymmetryList());
   DEBUG_CHECK_EQUAL(x.GetSymmetryList(), b1.NonAbelianBasis().GetSymmetryList());

   typedef typename const_iterator<IrredTensor<T, B1T, B2T, DefaultStructure> >::type iter_type;
   typedef typename const_inner_iterator<IrredTensor<T, B1T, B2T, DefaultStructure> >::type inner_iter_type;

   SymmetryList AbelianSList = b1.AbelianBasis().GetSymmetryList();

   QuantumNumber q = x.TransformsAs();

   // get the quantum number of how our abelian operator transforms
   QuantumNumber qAbelian = map_projection_to_quantum(p, AbelianSList);

   IrredTensor<T, B1T, B2T, DefaultStructure> Result(b1.AbelianBasis(), b2.AbelianBasis(), qAbelian);

   for (iter_type I = iterate(x); I; ++I)
   {
      QuantumNumber q1 = x.Basis1()[I.index()];
      QuantumNumbers::ProjectionList p1 = enumerate_projections(q1);
      for (unsigned p1i = 0; p1i < p1.size(); ++p1i)
      {
         int a1 = b1.Map(I.index(), p1i);
         if (a1 == -1)  // stop now if this projection is not needed
            continue;

         QuantumNumber q1Abelian = map_projection_to_quantum(p1[p1i], AbelianSList);

         for (inner_iter_type J = iterate(I); J; ++J)
         {
            QuantumNumber q2 = x.Basis2()[J.index2()];
            QuantumNumbers::ProjectionList p2 = enumerate_projections(q2);
            for (unsigned p2i = 0; p2i < p2.size(); ++p2i)
            {
               int a2 = b2.Map(J.index2(), p2i);
               if (a2 == -1)
                  continue;

               QuantumNumber q2Abelian = map_projection_to_quantum(p2[p2i], AbelianSList);

               // if this CG coefficient is identically zero, then bug out
               if (!is_transform_target(q2Abelian, qAbelian, q1Abelian))
                  continue;

               Result.data()(a1, a2) = clebsch_gordan(q2, q, q1, p2[p2i], p, p1[p1i]) * (*J);
            }
         }
      }
   }
   return Result;
}

// A DiagonalStructure will in general map into a non-diagonal
// structure here, unless the projection is scalar
template <typename T, typename B1T, typename B2T>
IrredTensor<T, B1T, B2T, DefaultStructure>
wigner_eckart(IrredTensor<T, B1T, B2T, DiagonalStructure> const& x, 
              Projection const& p, 
              WignerEckartBasis<B1T> const& b1, 
              WignerEckartBasis<B1T> const& b2)
{
   DEBUG_CHECK_EQUAL(b1.AbelianBasis().GetSymmetryList(),
                     b2.AbelianBasis().GetSymmetryList());
   DEBUG_CHECK_EQUAL(b1.NonAbelianBasis().GetSymmetryList(),
                     b2.NonAbelianBasis().GetSymmetryList());
   DEBUG_CHECK_EQUAL(x.GetSymmetryList(), b1.NonAbelianBasis().GetSymmetryList());

   typedef typename const_iterator<IrredTensor<T, B1T, B2T, DiagonalStructure> >::type iter_type;
   typedef typename const_inner_iterator<IrredTensor<T, B1T, B2T, DiagonalStructure> >::type inner_iter_type;

   SymmetryList AbelianSList = b1.AbelianBasis().GetSymmetryList();

   QuantumNumber q = x.TransformsAs();

   // get the quantum number of how our abelian operator transforms
   QuantumNumber qAbelian = map_projection_to_quantum(p, AbelianSList);

   IrredTensor<T, B1T, B2T, DefaultStructure> Result(b1.AbelianBasis(), b2.AbelianBasis(), qAbelian);

   for (iter_type I = iterate(x); I; ++I)
   {
      for (inner_iter_type J = iterate(I); J; ++J)
      {

	 QuantumNumber q1 = x.Basis1()[J.index1()];
	 QuantumNumbers::ProjectionList p1 = enumerate_projections(q1);
	 for (unsigned p1i = 0; p1i < p1.size(); ++p1i)
	 {
	    int a1 = b1.Map(J.index1(), p1i);
	    if (a1 == -1)  // stop now if this projection is not needed
	       continue;

	    QuantumNumber q1Abelian = map_projection_to_quantum(p1[p1i], AbelianSList);

            QuantumNumber q2 = x.Basis2()[J.index2()];
            QuantumNumbers::ProjectionList p2 = enumerate_projections(q2);
            for (unsigned p2i = 0; p2i < p2.size(); ++p2i)
            {
               int a2 = b2.Map(J.index2(), p2i);
               if (a2 == -1)
                  continue;

               QuantumNumber q2Abelian = map_projection_to_quantum(p2[p2i], AbelianSList);

               // if this CG coefficient is identically zero, then bug out
               if (!is_transform_target(q2Abelian, qAbelian, q1Abelian))
                  continue;

               Result.data()(a1, a2) = clebsch_gordan(q2, q, q1, p2[p2i], p, p1[p1i]) * (*J);
            }
         }
      }
   }
   return Result;
}


} // namespace Tensor
