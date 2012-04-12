// -*- C++ -*- $Id$

#if !defined(DMRG_H_HHVERTHYEHYYOIUDLE89P)
#define DMRG_H_HHVERTHYEHYYOIUDLE89P

#include "matrixproduct/mpstate.h"
#include "matrixproduct/splitoperator.h"
#include "matrixproduct/centerwavefunction.h"
#include "matrixproduct/operatorstack.h"
#include "siteoperator/sitebasis.h"
#include "common/conflist.h"
#include <fstream>
#include <boost/shared_ptr.hpp>

struct DMRG
{
   typedef CenterWavefunction WavefunctionType;
   typedef std::complex<double> complex;
   typedef OperatorStack<MPStateComponent> SuperblockOperator;
   typedef OperatorStack<MatrixOperator>   TransformOperator;

   DMRG() {}

   DMRG(CenterWavefunction const& Psi_, SplitOperator const& Ham_, bool Verbose = false);

   // Adds x to the 'orthogonal set', that we explicitly orthogonalize the
   // wavefunction against
   void AddOrthogonalState(CenterWavefunction x);

   void StartSweep(bool IncrementSweepNumber = true, double Broad = 0);
   
   void EndSweep();    // statistics for end of sweep

   void StartIteration();  // prepare statistics for start of iteration
   void EndIteration();    // statistics for end of iteration

   void CreateLogFiles(std::string const& BasePath, ConfList const& Conf);
   void RestoreLogFiles(std::string const& BasePath, ConfList const& Conf);

   int LeftSize() const { return Psi.LeftSize(); }
   int RightSize() const { return Psi.RightSize(); }

   // Calculates the groundstate, using Iterations number
   // of Lanczos steps
   double Solve(int Iterations);

   // Does a truncation and shifts the Center matrix to the right.
   // The new left basis is automatically expanded.
   void ShiftRightAndExpand();

   // Does a truncation and shifts the Center matrix to the left.
   // The new right basis is automatically expanded.
   void ShiftLeftAndExpand();

   // 'Expands' the left basis to cover all states.
   void ExpandLeft();

   // 'Expands' the right basis to cover all states.
   void ExpandRight();

   // returns the energy, calculated from a matrix-vector product
   double Energy() const;

   // returns ||psi - psi_prev||^2, where psi_prev is the previous
   // wavefunction (set to psi by StartSweep())
   double FidelityLoss() const;

   // Inserts sites into the lattice to the left and right of
   // the center matrix.  A delta shift is applied to the left matrices
   // such that the new wavefunction transforms as NewTarget.
   // The Hamiltonian matrix elements are recalculated from scratch.
   // The OrthogonalSet must be empty.  The NewTarget must be obtainable
   // from the current target by a positive delta shift.
   void InsertSitesDeltaShift(std::vector<SiteBasis> const& LeftSites,
                              std::vector<SiteBasis> const& RightSites,
                              QuantumNumber NewTarget,
                              LinearOperator const& NewHam);

   // This version inserts sites and makes the Center matrix a non-scalar tensor
   void InsertSites(std::vector<SiteBasis> const& LeftSites,
		    std::vector<SiteBasis> const& RightSites,
		    QuantumNumber NewTarget,
		    SplitOperator const& NewHam);

   CenterWavefunction& Wavefunction();
   CenterWavefunction const& Wavefunction() const;

   void PrepareConvergenceTest();
   bool IsConverged() const;

   TruncationInfo TruncateLeft(StatesInfo const& SInfo, double CFactor);
   TruncationInfo TruncateRight(StatesInfo const& SInfo, double CFactor);

   void DebugCheckBasis() const;

   SuperblockOperator HamMatrices;
   CenterWavefunction Psi;
   MatrixOperator PsiPrevC; // Center matrix of the old wavefunction projected onto new basis
   SplitOperator Ham;

   SplitOperator HamSquared;

   QuantumNumber Ident;

   std::vector<CenterWavefunction> Ortho;              // set of wavefunctions that we want to be
   std::vector<TransformOperator> PsiOrthoProjector;  // orthogonal to

   bool SaveDiscardedWavefunction; // if true, save the discarded wavefunctions on a right-moving sweep
   bool NormalizeWavefunction; // should we normalize the wavefunction after each truncation?
   boost::optional<double> LastOverlap;
   bool IsPsiConverged; 
   bool IsConvergedValid;
   bool TestConverged;            // true if we should calculate (H-E)^2 at end of sweep
   KeepListType KeepList;

   // some global statistics
   int TotalSweepNumber;          // this is incremented every sweep
   int TotalSweepRecNumber;       // this is incremented conditionally
   int TotalNumIterations;
   int TotalNumMultiplies;

   // some statistics, for current sweep
   int SweepNumIterations;
   int SweepSumStates;            // sum of IterationNumStates
   int SweepMaxStates;
   int SweepNumMultiplies;        // number of mat-vec multiplies this sweep
   double SweepEnergy;            // lowest energy seen this sweep
   double SweepTruncation;        // cumulative truncation error this sweep
   double SweepTruncationFrac;    // cumulative fractional-truncation error this sweep
   double SweepEntropy;           // maximum entropy this sweep
   double SweepStartTime;         // wall time at the start of the sweep
   double SweepTruncatedEnergy;   // sum of (E_0 - E_truncated) over the sweep
   double SweepEnergyError;       // standard error of the energy at each iteration
   double SweepLastMixFactor;     // the last used mix factor, for the .sweep log file

   // some statistics, for current iteration
   int IterationNumMultiplies;
   int IterationNumStates;
   double IterationEnergy;
   double IterationTruncation;
   double IterationTruncationFrac;
   double IterationEntropy;
   double IterationEnergyBeforeTrunc;
   std::vector<double> IterationEnergyVec;

   // Log files are not serialized, but initialized by CreateLogFiles or
   // RestoreLogFiles
   boost::shared_ptr<std::ofstream> EnergyLog, DiagLog, SweepLog, CpuLog, DensityLog;

   // Following items are not persistent, but loaded from the configuration

   // for convergence, we require 2 conditions:
   // Overlap < ConvergenceOverlapTruncationScale * SweepTruncation
   // and OverlapDifference < ConvergenceOverlapDifferenceOverlapScale * Overlap
   double ConvergenceOverlapTruncationScale;
   double ConvergenceOverlapDifferenceOverlapScale;
   double ConvergenceSweepTruncMin;
   double ConvergenceOverlapMin;
   bool MixUseEnvironment;
   bool UseDGKS;  // use DGKS correction in the lanczos for orthogonal states
   std::string Solver; // "davidson" or "lanczos", or "arnoldi"

   bool TwoStageTruncation;
   double ExtraFrac;
};

PStream::opstream& operator<<(PStream::opstream& out, DMRG const& d);
PStream::ipstream& operator>>(PStream::ipstream& out, DMRG& d);

#endif
