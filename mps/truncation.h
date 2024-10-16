// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mps/truncation.h
//
// Copyright (C) 2004-2024 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

// This header exists to provide some declarations of common classes to avoid
// recursive includes, since state_component.h uses functionality from
// density.h, but cannot include.

#if !defined(MPTOOLKIT_MPS_TRUNCATION_H)
#define MPTOOLKIT_MPS_TRUNCATION_H

#include <limits>
#include <iostream>
#include <quantumnumbers/quantumnumber.h>

// Structure to hold information on a particular density matrix eigenstate,
// or the square of a singular value
struct EigenInfo
{
   double Eigenvalue;
   int Subspace, Index;     // subspace/index into the basis
   QuantumNumbers::QuantumNumber Q;

   double Degree() const
   {
      return degree(Q);
   }

   double Weight() const
   {
      return degree(Q)*Eigenvalue;
   }

   double Entropy(double TotalWeight) const
   {
      if (TotalWeight < std::numeric_limits<double>::epsilon())
         return 0;
      double x = Eigenvalue/TotalWeight;
      return x > 0 ? -x * this->Degree() * log(x) : 0;
   }
   EigenInfo(double Eigen_, int QuantumNumber_, int LinearIndex_, QuantumNumbers::QuantumNumber Q_)
      : Eigenvalue(Eigen_), Subspace(QuantumNumber_), Index(LinearIndex_), Q(Q_) {}

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

      QuantumNumberEqualTo(QuantumNumbers::QuantumNumber const& q_)
         :q(q_) {}

      bool operator()(EigenInfo const& e) const
      {
         return e.Q == q;
      }

      QuantumNumbers::QuantumNumber q;
   };

};

//
// There are two ways we can sort the eigenvalues, by total weight or
// by eigenvalue (ignoring degeneracy).  Sorting by total weight
// is slightly better for finite systems, but possibly doesnt
// work well in iDMRG where it favours higher spin representations
// that increase the degeneracy of the groundstate unnecessarily.
// Default value is false.
//
// For finite DMRG, we can set EigenSortByWeight=true

extern bool EigenSortByWeight;

inline
bool EigenCompareWeight(EigenInfo const& x, EigenInfo const& y)
{
   return y.Eigenvalue*degree(y.Q) < x.Eigenvalue*degree(x.Q);
}

inline
bool EigenCompare(EigenInfo const& x, EigenInfo const& y)
{
   if (EigenSortByWeight)
      return y.Eigenvalue*degree(y.Q) < x.Eigenvalue*degree(x.Q);
   else
      return y.Eigenvalue < x.Eigenvalue;
}

// the information we use to determine which states to keep is getting bigger, make it a structure
// this is the information that we set to determine the truncation of the basis
struct StatesInfo
{
   static int const DefaultMaxStates = 1000000;  // we need a large default for the max states to keep

   int MinStates;
   int MaxStates;
   double TruncationCutoff;
   double EigenvalueCutoff;
   bool TruncateRelative;
   // if TruncateRelative is true, cutoff values are relative to total weight and range 0..1
   // otherwise the cutoff values are in absolute units

   StatesInfo()
      : MinStates(1),
        MaxStates(DefaultMaxStates),
        TruncationCutoff(0),
        EigenvalueCutoff(0),
        TruncateRelative(false){}

   explicit StatesInfo(int MaxStates_)
      : MinStates(1), MaxStates(MaxStates_), TruncationCutoff(0), EigenvalueCutoff(0), TruncateRelative(false)
   {}

   // helper function that returns true if the configuration of the StatesInfo is such that
   // we want to keep eigenvalues that are exactly zero
   bool KeepZeroEigenvalues() const { return TruncationCutoff < 0 && EigenvalueCutoff < 0; }
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
   int ExtraStates_;
   double TotalEntropy_;
   double KeptEntropy_;
   double SmallestKeptEigenvalue_;
   double LargestDiscardedEigenvalue_;

   TruncationInfo();

   // return the total weight of all states considered
   double TotalWeight() const { return TotalWeight_; }

   // return the weight of the states that were kept
   double KeptWeight() const { return KeptWeight_; }

   // return the total number of states that were considered
   int TotalStates() const { return TotalStates_; }

   // return the number of states that were kept
   int KeptStates() const { return KeptStates_; }

   // return the number of extra that were kept
   int ExtraStates() const { return ExtraStates_; }

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

inline
TruncationInfo::TruncationInfo()
   : TotalWeight_(0), KeptWeight_(0), TotalStates_(0), KeptStates_(0), ExtraStates_(0), TotalEntropy_(0), KeptEntropy_(0), SmallestKeptEigenvalue_(0), LargestDiscardedEigenvalue_(0)
{
}

// functor to get the cumulative sum of eigenvalues
template <class FwdIter>
double DensitySigma(FwdIter first, FwdIter last)
{
   double Acc = 0;
   while (first != last)
   {
      Acc += first->Eigenvalue * degree(first->Q);
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
         E -= x * log2(x) * degree(first->Q);
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
         E -= x * log(x) * degree(first->Q);
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
      double x = first->Eigenvalue * degree(first->Q);
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
   Info.ExtraStates_ = 0;
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
   Info.ExtraStates_ = 0;
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

// TruncateExtraStates - intended for finding optimal states to keep in an expanded environment.
// Keep at least NumStates states, with at least StatesPerSector states in each quantum number sector.
// StatesPerSector would normally be small, typically 0 or 1.
// If AllowZeroWeight is false, then never keep eigenvalues <= 0, even if we don't reach NumStates or
// StatesPerSector.
// Generally speaking, for pre-expansion it is harmless to keep eigenvalues that have a zero eigenvalue here,
// and it is possible that they can end up with non-zero weight in the wavefunction.
// Note that the implementation of the SVD returns random vectors for states that have a structurally zero singular value.
// On the other hand, for post-expansion states that have zero weight have no overlap with the
// system reduced density matrix, and therefore there is no point keeping them (i.e. they are not reachable by applying operators
// from the system MPO to the kept states).
template <typename FwdIter>
std::vector<EigenInfo>
TruncateExtraStates(FwdIter first, FwdIter last, int NumStates, int StatesPerSector, bool AllowZeroWeight)
{
   std::vector<EigenInfo> Result;
   Result.reserve(NumStates);

   std::map<QuantumNumbers::QuantumNumber, int> KeptStatesPerSector;
   // first part: keep the next NumStatesWithWeight in order from heighest weight,
   // as long as they have non-zero weight.
   auto f = first;
   while (NumStates > 0 && f != last && (AllowZeroWeight || f->Eigenvalue > 0.0))
   {
      Result.push_back(*f);
      --NumStates;
      ++KeptStatesPerSector[f->Q];
      ++f;
   }
   // at this point, we have reached NumStates states in Result, OR we've hit the end of the list and f == last, OR all remaining
   // eigenvalues are zero.

   // second part: ensire that we keep at least StatesPerSector states in each quantum number sector, but only keep
   // zero eigenvalues if StatesPerSectorAllowZeroWeight is true
   while (f != last)
   {
      if (KeptStatesPerSector[f->Q] < StatesPerSector && (AllowZeroWeight || f->Eigenvalue > 0.0))
      {
         Result.push_back(*f);
         ++KeptStatesPerSector[f->Q];
      }
      ++f;
   }

   return Result;
}

#endif
