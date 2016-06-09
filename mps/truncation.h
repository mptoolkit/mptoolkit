// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mps/truncation.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Reseach publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#if !defined(TRUNCATION_H_SDCHJYUTYFHLUH53247895732489)
#define TRUNCATION_H_SDCHJYUTYFHLUH53247895732489

#include <limits>
#include <iostream>

// Structure to hold information on a particular density matrix eigenstate,
// or the square of a singular value
struct EigenInfo
{
   double Eigenvalue;
   int Degree;              // degree of the representation
   int Subspace, Index;     // subspace/index into the basis

   double Weight() const
   {
      return Degree*Eigenvalue;
   }

   double Entropy(double TotalWeight) const
   {
      if (TotalWeight < std::numeric_limits<double>::epsilon())
         return 0;
      double x = Eigenvalue/TotalWeight;
      return x > 0 ? -x * Degree * log(x) : 0;
   }

   EigenInfo(double Eigen_, int Degree_, int QuantumNumber_, int LinearIndex_)
      : Eigenvalue(Eigen_), Degree(Degree_), Subspace(QuantumNumber_), Index(LinearIndex_) {}

   bool operator==(EigenInfo const& Other) const
   {
      return Subspace == Other.Subspace && Index == Other.Index;
   }

   bool operator!=(EigenInfo const& Other) const
   {
      return Subspace != Other.Subspace || Index != Other.Index;
   }

   bool operator<(EigenInfo const& Other) const
   {
      return Subspace < Other.Subspace || (Subspace == Other.Subspace && Index < Other.Index);
   }

   // Functor to return true if the EigenInfo has a given quantum number
   struct QuantumNumberEqualTo
   {
      typedef bool result_type;
   
      QuantumNumberEqualTo(VectorBasis const& B_, QuantumNumbers::QuantumNumber const& q_)
	 : B(B_), q(q_) {}

      bool operator()(EigenInfo const& e) const
      {
	 return B[e.Subspace] == q;
      }
   
      VectorBasis const& B;
      QuantumNumbers::QuantumNumber q;
   };

};

//
// There are two ways we can sort the eigenvalues, by total weight or
// by eigenvalue (ignoring degeneracy).  Sorting by total weight
// is slightly better for finite systems, but possibly doesnt
// work well in iDMRG where it favours higher spin representations
// that increase the degeneracy of the groundstate unnecessarily.
//

inline
bool EigenCompareWeight(EigenInfo const& x, EigenInfo const& y)
{
   return y.Eigenvalue*y.Degree < x.Eigenvalue*x.Degree;
}

inline
bool EigenCompare(EigenInfo const& x, EigenInfo const& y)
{
   return y.Eigenvalue < x.Eigenvalue;
}

// the information we use to determine which states to keep is getting bigger, make it a structure
// this is the information that we set to determine the truncation of the basis
struct StatesInfo
{
   static int const DefaultMaxStates = 100000;  // we need a large default for the max states to keep

   int MinStates;
   int MaxStates;
   double TruncationCutoff;
   double EigenvalueCutoff;
   bool TruncateRelative;
   // if TruncateRelative is true, cutoff values are relative to total weight and range 0..1
   // otherwise the cutoff values are in absolute units

   StatesInfo()
      : MinStates(0), 
        MaxStates(DefaultMaxStates), 
        TruncationCutoff(0), 
        EigenvalueCutoff(0),
        TruncateRelative(false){}
};

inline
std::ostream& operator<<(std::ostream& out, StatesInfo const& s)
{
   return out << "MinStates: " << s.MinStates 
              << " MaxStates: " << s.MaxStates
              << " TruncationCutoff: " << s.TruncationCutoff 
              << " EigenvalueCutoff: " << s.EigenvalueCutoff;
}

inline
PStream::opstream& operator<<(PStream::opstream& out, StatesInfo const& s)
{
   out << s.MinStates << s.MaxStates << s.TruncationCutoff << s.EigenvalueCutoff;
   return out;
}

inline
PStream::ipstream& operator>>(PStream::ipstream& in, StatesInfo& s)
{
   in >> s.MinStates >> s.MaxStates >> s.TruncationCutoff >> s.EigenvalueCutoff;
   return in;
}

// TruncationInfo is used to return statistics on the truncation that is performed
struct TruncationInfo
{
   double TotalWeight_;
   double KeptWeight_;
   int TotalStates_;
   int KeptStates_;
   double TotalEntropy_;
   double KeptEntropy_;
   double SmallestKeptEigenvalue_;
   double LargestDiscardedEigenvalue_;

   // return the total weight of all states considered
   double TotalWeight() const { return TotalWeight_; }

   // return the weight of the states that were kept
   double KeptWeight() const { return KeptWeight_; }

   // return the total number of states that were considered
   int TotalStates() const { return TotalStates_; }

   // return the number of states that were kept
   int KeptStates() const { return KeptStates_; }

   // return the total entropy (natural units) of all of the states
   double TotalEntropy() const { return TotalEntropy_; }

   // return the entropy (natural units) of the kept states
   double KeptEntropy() const { return KeptEntropy_; }

   // return the normalized truncation error
   double TruncationError() const { return TotalWeight_ < std::numeric_limits<double>::epsilon() ? 0.0 :
      1.0 - KeptWeight_ / TotalWeight_; }

   // return the smallest eigenvalue of the kept states
   double SmallestKeptEigenvalue() const { return SmallestKeptEigenvalue_; }

   // return the largest eigenvalue of the discarded states
   double LargestDiscardedEigenvalue() const { return LargestDiscardedEigenvalue_; }
};

// functor to get the cumulative sum of eigenvalues
template <class FwdIter>
double DensitySigma(FwdIter first, FwdIter last)
{
   double Acc = 0;
   while (first != last)
   {
      Acc += first->Eigenvalue * first->Degree;
      ++first;
   }
   return Acc;
}

// functor to get the entropy of a set of eigenvalues
template <typename FwdIter>
double DensityEntropyBase2(FwdIter first, FwdIter last, double EigenSum)
{
   if (EigenSum < std::numeric_limits<double>::epsilon())
      return 0;
   double E = 0;
   while (first != last)
   {
      double x = first->Eigenvalue/EigenSum;
      if (x > 0)
	 E -= x * log2(x) * first->Degree;
      else if (x < -1E-10)
      {
         WARNING("Entropy: eigenvalue is negative and big!")(x)(first->Eigenvalue)(EigenSum);
      }
      ++first;
   }
   return E;
}

template <typename FwdIter>
double DensityEntropyBaseE(FwdIter first, FwdIter last, double EigenSum)
{
   if (EigenSum < std::numeric_limits<double>::epsilon())
      return 0;
   double E = 0;
   while (first != last)
   {
      double x = first->Eigenvalue/EigenSum;
      if (x > 0)
	 E -= x * log(x) * first->Degree;
      else if (x < -1E-10)
      {
         WARNING("Entropy: eigenvalue is negative and big!")(x)(first->Eigenvalue)(EigenSum);
      }
      ++first;
   }
   return E;
}

template <typename FwdIter>
double DensityEntropy(FwdIter first, FwdIter last, double EigenSum, bool Base2 = false)
{
   return Base2 ? 
      DensityEntropyBase2(first, last, EigenSum)
      : 
      DensityEntropyBaseE(first, last, EigenSum);
}

// functor to get the normalized truncation error of a set of eigenvalues
template <typename FwdIter>
double DensityTruncation(FwdIter first, FwdIter last, double EigenSum)
{
   double E = 0;
   while (first != last)
   {
      double x = first->Eigenvalue * first->Degree;
      if (x > 0)
	 E += x;

      ++first;
   }
   return 1.0 - E/EigenSum;
}

// Truncate to the given truncation error, normalized to norm 1
template <typename FwdIter>
FwdIter 
TruncateFixTruncationErrorRelative(FwdIter first, FwdIter last, 
                                   StatesInfo const& States, TruncationInfo& Info)
{
   Info.TotalWeight_ = DensitySigma(first, last);
   Info.TotalStates_ = std::distance(first, last);
   Info.TotalEntropy_ = DensityEntropy(first, last, Info.TotalWeight_);
   Info.KeptStates_ = 0;
   Info.KeptWeight_ = 0;
   Info.KeptEntropy_ = 0;
   Info.SmallestKeptEigenvalue_ = 0;
   double const RequiredSigma = Info.TotalWeight_ * (1.0 - States.TruncationCutoff);
   double const RequiredWeight = Info.TotalWeight_ * States.EigenvalueCutoff;
   while (first != last && (Info.KeptStates_ < States.MinStates || 
                            (first->Weight() > RequiredWeight && 
			     Info.KeptWeight_ < RequiredSigma && 
			     Info.KeptStates_ < States.MaxStates)))
   {
      ++Info.KeptStates_;
      Info.KeptWeight_ += first->Weight();
      Info.KeptEntropy_ += first->Entropy(Info.TotalWeight_);
      Info.SmallestKeptEigenvalue_ = first->Weight();
      ++first;
   }
   if (first == last)
      Info.LargestDiscardedEigenvalue_ = 0;
   else
      Info.LargestDiscardedEigenvalue_ = first->Weight();
   return first;
}

// Truncate to the given truncation error, absolute normalization
template <typename FwdIter>
FwdIter 
TruncateFixTruncationErrorAbsolute(FwdIter first, FwdIter last, 
                                   StatesInfo const& States, TruncationInfo& Info)
{
   Info.TotalWeight_ = DensitySigma(first, last);
   Info.TotalStates_ = std::distance(first, last);
   Info.TotalEntropy_ = DensityEntropy(first, last, Info.TotalWeight_);
   Info.KeptStates_ = 0;
   Info.KeptWeight_ = 0;
   Info.KeptEntropy_ = 0;
   Info.SmallestKeptEigenvalue_ = 0;
   double const RequiredSigma = Info.TotalWeight_ - States.TruncationCutoff;
   double const RequiredWeight = States.EigenvalueCutoff;
   while (first != last && (Info.KeptStates_ < States.MinStates || 
                            (first->Weight() > RequiredWeight && 
			     Info.KeptWeight_ < RequiredSigma && 
			     Info.KeptStates_ < States.MaxStates)))
   {
      ++Info.KeptStates_;
      Info.KeptWeight_ += first->Weight();
      Info.KeptEntropy_ += first->Entropy(Info.TotalWeight_);
      Info.SmallestKeptEigenvalue_ = first->Weight();
      ++first;
   }
   if (first == last)
      Info.LargestDiscardedEigenvalue_ = 0;
   else
      Info.LargestDiscardedEigenvalue_ = first->Weight();
   return first;
}

template <typename FwdIter>
FwdIter 
TruncateFixTruncationError(FwdIter first, FwdIter last, 
                           StatesInfo const& States, TruncationInfo& Info)
{
   return States.TruncateRelative 
      ? TruncateFixTruncationErrorRelative(first, last, States, Info)
      : TruncateFixTruncationErrorAbsolute(first, last, States, Info);
}

#endif
