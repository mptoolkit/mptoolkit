// -*- C++ -*- $Id$

#include "vectorbasissum.h"

std::pair<DiagonalProjection, DiagonalProjection>
basis_sum(VectorBasis const& B1, VectorBasis const& B2)
{
   DEBUG_CHECK_EQUAL(B1.GetSymmetryList(), B2.GetSymmetryList());
   // Get the dimension of the final subspaces
   typedef std::map<QuantumNumbers::QuantumNumber, int> SubDimMapType;
   SubDimMapType SubDim;
   for (unsigned i = 0; i < B1.size(); ++i)
      SubDim[B1[i]] += B1.dim(i);
   for (unsigned i = 0; i < B2.size(); ++i)
      SubDim[B2[i]] += B2.dim(i);

   // construct the new basis
   VectorBasis Result(B1.GetSymmetryList());
   std::vector<int> B1Map(B1.size());
   std::vector<int> B2Map(B2.size());
   for (SubDimMapType::const_iterator I = SubDim.begin(); 
	I != SubDim.end(); ++I)
   {
      Result.push_back(I->first, I->second);
   }

   // Fill the projectors
   DiagonalProjection B1Proj(B1, Result);
   std::vector<int> SubIndex(Result.size(), 0);
   for (unsigned i = 0; i < B1.size(); ++i)
   {
      // Get the index j with B1[i] == Result[j]
      unsigned j = 0;
      while (Result[j] != B1[i])
	 ++j;

      B1Proj(i,j) = LinearAlgebra::Matrix<double>(B1.dim(i), Result.dim(j), 0.0);

      for (int k = 0; k < B1.dim(i); ++k)
	 B1Proj(i,j)(k,SubIndex[j]+k) = 1.0;
      SubIndex[j] += B1.dim(i);
   }

   DiagonalProjection B2Proj(B2, Result);
   for (unsigned i = 0; i < B2.size(); ++i)
   {
      // Get the index j with B2[i] == Result[j]
      unsigned j = 0;
      while (Result[j] != B2[i])
	 ++j;

      B2Proj(i,j) = LinearAlgebra::Matrix<double>(B2.dim(i), Result.dim(j), 0.0);

      for (int k = 0; k < B2.dim(i); ++k)
	 B2Proj(i,j)(k,SubIndex[j]+k) = 1.0;
      SubIndex[j] += B2.dim(i);
   }

   return std::pair<DiagonalProjection, DiagonalProjection>(B1Proj, B2Proj);
}

