// -*- C++ -*-

#include "wavefunction/mpwavefunction.h"
#include "lattice/latticesite.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "mp/copyright.h"
#include <boost/program_options.hpp>
#include "common/environment.h"
#include "common/prog_options.h"
#include "interface/inittemp.h"
#include "tensor/tensor_eigen.h"
#include "lattice/unitcell-parser.h"
#include "wavefunction/operator_actions.h"
#include "lattice/infinite-parser.h"



namespace prog_opt = boost::program_options;

void PrintFormat(std::complex<double> const& Value, bool ShowRealPart, bool ShowImagPart, 
		 bool ShowMagnitude, bool ShowArgument,
		 bool ShowRadians)
{
   if (ShowRealPart)
   {
      std::cout << std::setw(20) << Value.real() << "    ";
   }
   if (ShowImagPart)
   {
      std::cout << std::setw(20) << Value.imag() << "    ";
   }
   if (ShowMagnitude)
   {
      std::cout << std::setw(20) << LinearAlgebra::norm_frob(Value) << "    ";
   }
   if (ShowArgument)
   {
      double Arg =  std::atan2(Value.imag(), Value.real());
      if (!ShowRadians)
	 Arg *= 180.0 / math_const::pi;
      std::cout << std::setw(20) << Arg << "    ";
   }
}

void ShowHeading(bool ShowRealPart, bool ShowImagPart, 
		 bool ShowMagnitude, bool ShowArgument, bool ShowRadians)
{
   if (ShowRealPart)
      std::cout << "#real                   ";
   if (ShowImagPart)
      std::cout << "#imag                   ";
   if (ShowMagnitude)
      std::cout << "#magnitude              ";
   if (ShowArgument)
      std::cout << "#argument" << (ShowRadians ? "(rad)" : "(deg)") << "          ";
}

// copied from mp-ispectrum.cpp
// inject_left for a FiniteMPO.  This can have support on multiple wavefunction unit cells

// TODO: convert this to act on InfiniteWavefunctionLeft and UnitCellMPO

// On input, m is an operator acting on Psi.Basis1()
// result' is an operator acting on Psi.Basis2()
MatrixOperator
contract_from_left(MatrixOperator const& m, 
		   InfiniteWavefunctionLeft const& Psi,
		   FiniteMPO const& Op)
{
   CHECK(Op.size() % Psi.size() == 0)(Op.size())(Psi.size());
   DEBUG_CHECK_EQUAL(m.Basis2(), Psi.Basis1());
   if (Op.is_null())
      return MatrixOperator();

   // we currently only support simple irreducible operators
   DEBUG_CHECK_EQUAL(m.Basis1(), Psi.Basis1());
   CHECK_EQUAL(Op.Basis1().size(), 1);
   CHECK_EQUAL(Op.Basis2().size(), 1);
   CHECK_EQUAL(Op.Basis1()[0], m.TransformsAs());
   StateComponent E(Op.Basis1(), m.Basis1(), m.Basis2());
   E[0] = m;
   E.debug_check_structure();
   InfiniteWavefunctionLeft::const_mps_iterator I = Psi.begin();
   FiniteMPO::const_iterator OpIter = Op.begin();
   while (OpIter != Op.end())
   {
      if (I == Psi.end())
      {
	 I = Psi.begin();
	 E = delta_shift(E, Psi.qshift());
      }
      E = contract_from_left(*OpIter, herm(*I), E, *I);
      ++I; ++OpIter;
   }
   return E[0];
}

MatrixOperator
contract_from_right(MatrixOperator const& m, 
		    LinearWavefunction const& Psi1,
		    QuantumNumbers::QuantumNumber const& QShift,
		    GenericMPO const& Op, 
		    LinearWavefunction const& Psi2)
{
   PRECONDITION_EQUAL(Psi1.size(), Psi2.size());
   if (Op.is_null())
      return MatrixOperator();

   // we currently only support simple irreducible operators
   CHECK_EQUAL(Op.Basis1().size(), 1);
   CHECK_EQUAL(Op.Basis2().size(), 1);
   StateComponent E(Op.Basis2(), m.Basis1(), m.Basis2());
   E[0] = m;
   MatrixOperator Result = m;
   LinearWavefunction::const_iterator I1 = Psi1.end();
   LinearWavefunction::const_iterator I2 = Psi2.end();
   GenericMPO::const_iterator OpIter = Op.end();
   while (OpIter != Op.begin())
   {
      if (I1 == Psi1.begin())
      {
         I1 = Psi1.end();
         I2 = Psi2.end();
         E = delta_shift(E, adjoint(QShift));
      }
      --I1; --I2; --OpIter;
      E = contract_from_right(herm(*OpIter), *I1, E, herm(*I2));
   }
   return E[0];
}

MatrixOperator
contract_from_right(MatrixOperator const& m, 
		    InfiniteWavefunctionLeft const& Psi,
		    GenericMPO const& Op)
{
   LinearWavefunction PsiLinear = get_left_canonical(Psi).first;
   return contract_from_right(m, PsiLinear, Psi.qshift(), Op, PsiLinear);
}

int RoundUp(int numToRound, int multiple)
{
    int remainder = std::abs(numToRound) % multiple;
    if (remainder == 0)
        return numToRound;

    if (numToRound < 0)
       return -(std::abs(numToRound) - remainder);
    else
        return numToRound + multiple - remainder;
}

int RoundDown(int numToRound, int multiple)
{
    int remainder = std::abs(numToRound) % multiple;
    if (remainder == 0)
        return numToRound;

    if (numToRound < 0)
       return -(std::abs(numToRound) + multiple - remainder);
    else
        return numToRound - remainder;
}

std::complex<double>
mean(std::vector<std::complex<double>> const& v)
{
   std::complex<double> sum = 0.0;
   for (auto x : v)
      sum += x;
   return (1.0 / v.size()) * sum;
}

double
variance(std::vector<std::complex<double>> const& v)
{
   std::complex<double> me = mean(v);
   double var = 0.0;
   for (auto x : v)
      var += norm_frob_sq(me-x);
   return var / v.size();
}

int main(int argc, char** argv)
{
   try
   {
      bool ShowRealPart = false, ShowImagPart = false, ShowMagnitude = false;
      bool ShowCartesian = false, ShowPolar = false, ShowArgument = false;
      bool ShowRadians = false;
      int UnitCellSize = 0;
      std::string StringOpStr;
      std::string PsiStr;
      std::string Op1Str, Op2Str;
      int Length = 100;
      bool IncludeFirst = false;
      int Verbose = 0;
      bool Quiet = false;
      bool UseTempFile = false;
      bool Average = false;
      bool Connected = false;
      bool Separate = false;
      std::string LatticeFile;
      bool IncludeOverlap = true;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("cart,c", prog_opt::bool_switch(&ShowCartesian),
	  "show the result in cartesian coordinates [equivalent to --real --imag]")
	 ("polar,p", prog_opt::bool_switch(&ShowPolar),
	  "show the result in polar coodinates [equivalent to --mag --arg]")
	 ("real,r", prog_opt::bool_switch(&ShowRealPart),
	  "display only the real part of the result")
	 ("imag,i", prog_opt::bool_switch(&ShowImagPart),
	  "display only the imaginary part of the result")
	 ("mag,m", prog_opt::bool_switch(&ShowMagnitude),
          "display the magnitude of the result")
         ("arg,a", prog_opt::bool_switch(&ShowArgument),
          "display the argument (angle) of the result")
	 ("radians", prog_opt::bool_switch(&ShowRadians),
	  "display the argument in radians instead of degrees")
	 ("lattice,l", prog_opt::value(&LatticeFile),
	  "Use this lattice file, instead of specifying the lattice for each operator")
         ("string,s", prog_opt::value(&StringOpStr),
          "calculate a string correlation with the given operator S, as operator1 "
	  "\u2297 S \u2297  S .... S \u2297 operator2")
         ("includefirst,d", prog_opt::bool_switch(&IncludeFirst),
          "when calculating a string correlation, apply the string operator also "
          "to the first site, as operator1*S")
	 ("unitcell,u", prog_opt::value(&UnitCellSize),
	  "Use this unit cell size [default lattice unit cell size]")
         ("length,n", prog_opt::value(&Length),
          "calculate correlator to this length [default 100]")
	 ("average", prog_opt::bool_switch(&Average),
	  "average the correlation over shifts by the unitcell")
	 ("connected", prog_opt::bool_switch(&Connected),
	  "calculate the connected correlation <AB> - <A><B> [generally not compatible with --string]")
	 ("separate", prog_opt::bool_switch(&Separate),
	  "print the 'connected' part <A><B> as a separate column [not yet implemented]")
         ("tempfile", prog_opt::bool_switch(&UseTempFile),
          "a temporary data file for workspace (path set by environment MP_BINPATH)")
         ("quiet", prog_opt::bool_switch(&Quiet),
          "don't show the preamble and column headings")
         ("verbose,v", prog_opt_ext::accum_value(&Verbose),
          "Extra debug output [can be used multiple times])")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi", prog_opt::value(&PsiStr), "psi")
         ("op1", prog_opt::value(&Op1Str), "op1")
         ("op2", prog_opt::value(&Op2Str), "op2")
         ;

      prog_opt::positional_options_description p;
      p.add("psi", 1);
      p.add("op1", 1);
      p.add("op2", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);    

      if (vm.count("help") > 0 || vm.count("op2") == 0)
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options] <psi> <operator1> <operator2>\n";
         std::cerr << desc << '\n';
         return 1;
      }

      // check consistency of some options
      if (Connected && vm.count("string"))
      {
	 std::cerr << "mp-icorrelation: fatal: cannot construct a connected string correlation function!\n";
	 return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      if (!Quiet)
	 print_preamble(std::cout, argc, argv);

      // If no output switches are used, default to showing everything
      if (!ShowRealPart && !ShowImagPart && !ShowMagnitude
	  && !ShowCartesian && !ShowPolar && !ShowArgument)
      {
	 ShowCartesian = true;
	 ShowPolar = true;
      }
      
      if (ShowCartesian)
      {
	 ShowRealPart = true;
	 ShowImagPart = true;
      }
      if (ShowPolar)
      {
	 ShowMagnitude = true;
	 ShowArgument = true;
      }
      if (ShowRadians)
	 ShowArgument = true;      

      // Load the wavefunction
      pvalue_ptr<MPWavefunction> PsiPtr;
      if (UseTempFile)
      {
	  mp_pheap::InitializeTempPHeap(Verbose);
	  PsiPtr = pheap::ImportHeap(PsiStr);
      }
      else
      {
         PsiPtr = pheap::OpenPersistent(PsiStr, mp_pheap::CacheSize(), true);
      }
      InfiniteWavefunctionLeft Psi = PsiPtr->get<InfiniteWavefunctionLeft>();

      // If a lattice was explicitly specified, then load it
      pvalue_ptr<InfiniteLattice> Lattice;
      if (!LatticeFile.empty())
	 Lattice = pheap::ImportHeap(LatticeFile);

      // Load the operators
      UnitCellMPO Op1;
      int Lattice1UnitCellSize;
      {
	 if (Lattice)
	 {
	    Op1 = ParseUnitCellOperator(Lattice->GetUnitCell(), 0, Op1Str);
	    Lattice1UnitCellSize = Lattice->GetUnitCell().size();
	 }
	 else
	 {
	    InfiniteLattice Lattice1;
	    boost::tie(Op1, Lattice1) = ParseUnitCellOperatorAndLattice(Op1Str);
	    Lattice1UnitCellSize = Lattice1.GetUnitCell().size();
	 }
      }

      UnitCellMPO Op2;
      int Lattice2UnitCellSize;
      {
	 if (Lattice)
	 {
	    Op2 = ParseUnitCellOperator(Lattice->GetUnitCell(), 0, Op2Str);
	    Lattice2UnitCellSize = Lattice->GetUnitCell().size();
	 }
	 else
	 {
	    InfiniteLattice Lattice2;
	    boost::tie(Op2, Lattice2) = ParseUnitCellOperatorAndLattice(Op2Str);
	    Lattice2UnitCellSize = Lattice2.GetUnitCell().size();
	 }
      }

      // consistency checks
      if (Op1.Commute() * Op2.Commute() != LatticeCommute::Bosonic)
      {
	 std::cerr << "mp-icorrelation: fatal: cannot mix bosonic and fermionic operators!\n";
	 return 1;
      }
      if (adjoint(Op1.TransformsAs()) != Op2.TransformsAs())
      {
	 std::cerr << "mp-icorrelation: fatal: product of Op1 and Op2 is not a scalar!\n";
	 return 1;
      }

      if (Connected)
      {
	 if (Op1.Commute() != LatticeCommute::Bosonic)
	 {
	    std::cerr << "mp-icorrelation: fatal: cannot construct a connected correlator of a non-bosonic operator!\n";
	    return 1;
	 }
	 if (is_scalar(Op1.TransformsAs()))
	 {
	    std::cerr << "mp-icorrelation: fatal: cannot construct a connected correlator of a non-scalar operator!\n";
	    return 1;
	 }
      }

      // the string operator
      // Note that we don't need to shift the quantum numbers of the string operator,
      // because it is always multiplied by some other operator (eg the Operator2 JW string)
      // that handles properly the quantum numbers.
      ProductMPO StringOp;
      int StringSize = 0;
      int LatticeStringUnitCellSize = 0;
      if (vm.count("string"))
      {
	 if (Lattice)
	 {
	    StringOp = ParseProductOperator(*Lattice, StringOpStr);
	 }
	 else
	 {
	    InfiniteLattice StringLattice;
	    boost::tie(StringOp, StringLattice) = ParseProductOperatorAndLattice(StringOpStr);
	    LatticeStringUnitCellSize = StringLattice.GetUnitCell().size();
	 }
	 StringSize = StringOp.size();
      }
      else
      {
	 // No specified string operator, use the identity operator.
	 LatticeStringUnitCellSize = std::max(Lattice1UnitCellSize, Lattice2UnitCellSize);
	 StringOp = ProductMPO::make_identity(Basis1FromSiteList(*Op1.GetSiteList()));
	 StringSize = LatticeStringUnitCellSize;
      }

      // Set the unit cell size, if it wasn't specified on the command line
      if (UnitCellSize == 0)
      {
	 UnitCellSize = StringSize;
      }

      // make sure that the other lattice sizes divide evenly the UnitCellSize
      if (UnitCellSize % Lattice1UnitCellSize != 0)
      {
	 std::cerr << "mp-icorrelation: fatal: operator1 lattice unit cell size must divide the prime unit cell size.\n";
	 return 1;
      }
      if (UnitCellSize % Lattice2UnitCellSize != 0)
      {
	 std::cerr << "mp-icorrelation: fatal: operator2 lattice unit cell size must divide the prime unit cell size.\n";
	 return 1;
      }
      if (UnitCellSize % StringSize != 0)
      {
	 std::cerr << "mp-icorrelation: fatal: string operator lattice unit cell size must divide the prime unit cell size.\n";
	 return 1;
      }

      // make the wavefunction a multiple of the unit cell size
      if (Psi.size() % UnitCellSize != 0)
      {
	 int sz = statistics::lcm(Psi.size(), UnitCellSize);
	 std::cerr << "mp-icorrelation: warning: extending wavefunction to size " << sz
		   << " to be commensurate to the unit cell size.\n";
	 Psi = repeat(Psi, sz / Psi.size());
      }


      // TODO: better idea is to 
      // require that the unit cell of Op1 and Op2 are the same, and divide UnitCellSize.
      // Then we dont have to care about the separate Latice{1|2}UnitCellSize

      // Fix the StringOp size to the unit cell size
      StringOp = repeat(StringOp, UnitCellSize / StringOp.size());

      // get a finite MPO version of the string operator
      UnitCellMPO StringOpPerUnitCell = UnitCellMPO(Op1.GetSiteList(), FiniteMPO(StringOp), 
						    LatticeCommute::Bosonic, 0);

      MatrixOperator LeftIdent = MatrixOperator::make_identity(Psi.Basis1());


      int const PsiUnitCellSize = Psi.size() / UnitCellSize;


      // For generality we keep track of the offsets (with respect to the wavefunction unit cell)
      // that we want to use for the left and right operators separately.
      std::vector<int> LeftOffsets;
      std::vector<int> RightOffsets;
      for (int i = 0; i < Psi.size()/UnitCellSize; ++i)
      {
	 LeftOffsets.push_back(i);
	 RightOffsets.push_back(i);
      }

      // the Op1 matrices, indexed by the LeftOffset
      // Also keep track of the offset, measured in UnitCellSize, from the
      // logical '0' point of the operator to the edge of the operator.
      std::map<int, MatrixOperator> LeftOperator;
      std::map<int, int> LeftOperatorOffsetUnits;

      // The JW string for Op2
      FiniteMPO JW = Op2.GetJWStringUnit();
      JW = repeat(JW, UnitCellSize/JW.size());

      FiniteMPO WavefunctionStringJW = StringOpPerUnitCell.MPO() * JW;
      WavefunctionStringJW = repeat(WavefunctionStringJW, Psi.size()/WavefunctionStringJW.size());

      // the density matrix at the unit cell boundary
      MatrixOperator Rho = Psi.lambda(Psi.size());
      Rho = Rho*Rho;

      // for connected correlators
      std::map<int, std::complex<double>> ExpectOp1;
      std::map<int, std::complex<double>> ExpectOp2;

      if (Verbose > 0)
	 std::cerr << "Constructing matrices for Operator1\n";
      for (int x : LeftOffsets)
      {
	 UnitCellMPO Op = translate(Op1, x*UnitCellSize);
	 Op.ExtendToCoverUnitCell(Psi.size());
	 // incorporate the string operator, if necessary
	 for (int i = (IncludeFirst ? x : x+1); i < (Op.size()+Op.offset())/UnitCellSize; ++i)
	 {
	    Op = Op * translate(StringOpPerUnitCell, i*UnitCellSize);
	 }
	 // incorporate Op2 JW string, which goes over the entire length of the operator

	 LeftOperator[x] = contract_from_left(MatrixOperator::make_identity(Psi.Basis1()),
					      Psi, Op.MPO() * repeat(JW, Op.size()/JW.size()));
	 LeftOperatorOffsetUnits[x] = (Op.size()+Op.offset())/UnitCellSize - x;

	 if (Connected)
	    ExpectOp1[x] = inner_prod(LeftOperator[x], Rho);
      }

      // the Op2 matrices, indexed by the RightOffset
      // here the OffsetUnits measures the distance from the left edge of the operator
      // to the '0' point, in units of the UnitCellSize
      std::map<int, MatrixOperator> RightOperator;
      std::map<int, int> RightOperatorOffsetUnits;

      if (Verbose > 0)
	 std::cerr << "Constructing matrices for Operator2\n";

      for (int y : RightOffsets)
      {
	 UnitCellMPO Op = translate(Op2, y*UnitCellSize);
	 Op.ExtendToCoverUnitCell(Psi.size());
	 int Offset = Op.offset() + y*UnitCellSize;
	 // incorporate the string operator, if necessary
	 for (int i = y-1; i >= (Op.offset()/UnitCellSize); --i)
	 {
	    Op = Op * translate(StringOpPerUnitCell, i*UnitCellSize);
	 }
	 RightOperator[y] = contract_from_right(Rho, Psi, Op.MPO());
	 RightOperator[y].delta_shift(adjoint(Psi.qshift()));
	 RightOperatorOffsetUnits[y] = y - Op.offset()/UnitCellSize;

	 if (Connected)
	    ExpectOp2[y] = trace(RightOperator[y]);
      }

      // show the heading
      if (!Quiet)
      {
	 if (Average)
	 {
	    std::cout << "#n    ";
	    ShowHeading(ShowRealPart, ShowImagPart, ShowMagnitude, ShowArgument, ShowRadians);
	    std::cout << "#nsamples  #stdev";
	 }
	 else
	 {
	    std::cout << "#x    #y    ";
	    ShowHeading(ShowRealPart, ShowImagPart, ShowMagnitude, ShowArgument, ShowRadians);
	 }
	 std::cout << std::endl;
      }

      // for calculating averages over offsets
      std::vector<std::complex<double>> Av;

      // n is the length of the correlation
      for (int n = 0; n <= Length; ++n)
      {
	 // x and y offsets
	 for (int x : LeftOffsets)
	 {
	    for (int y : RightOffsets)
	    {
	       // Is it possible to make a correlation of length n out of these offsets?

	       if ((x+n) % PsiUnitCellSize != y)
		  continue;

	       // Otherwise we have a valid correlation

	       // Avoid calculating overlapping terms if we don't need to
	       if (!IncludeOverlap && LeftOperatorOffsetUnits[x] + RightOperatorOffsetUnits[y] > n)
		  continue;

	       if (Verbose > 1)
		  std::cerr << "correlator at x=" << x << " y=" << y << '\n';

	       // the value of our correlator
	       std::complex<double> e;

	       // Is it an overlapping term?

	       if (LeftOperatorOffsetUnits[x] + RightOperatorOffsetUnits[y] > n)
	       {
		  if (Verbose > 0)
		     std::cerr << "Term at " << x << "," << (x+n) << " has overlapping operators\n";
		  // the term is overlapping, calculate it 'by hand'
		  UnitCellMPO MyOp1 = translate(Op1, x*UnitCellSize);
		  // incorporate the string operator
		  for (int Shift = IncludeFirst ? 0 : 1; Shift < n; ++Shift)
		  {
		     MyOp1 = MyOp1 * translate(StringOpPerUnitCell, (Shift+x)*UnitCellSize);
		  }

		  UnitCellMPO MyOp2 = translate(Op2, (x+n)*UnitCellSize);
		  
		  UnitCellMPO Op = dot(MyOp1, MyOp2);
		  
		  // Evaluate the MPO
		  Op.ExtendToCoverUnitCell(Psi.size());
		  e = expectation(Psi, Op.MPO());
	       }
	       else
	       {
		  // Do we need to add another wavefunction unit cell?
		  while (LeftOperatorOffsetUnits[x] + RightOperatorOffsetUnits[y] < n)
		  {
		     // add another unit cell to LeftOperator[x]
		     if (Verbose > 0)
			std::cerr << "Adding another cell to operator at offset " << x << '\n';
		     LeftOperator[x] = contract_from_left(delta_shift(LeftOperator[x], Psi.qshift()), 
							  Psi, WavefunctionStringJW);
		     LeftOperatorOffsetUnits[x] += WavefunctionStringJW.size()/UnitCellSize;
		  }
		  CHECK_EQUAL(LeftOperatorOffsetUnits[x] + RightOperatorOffsetUnits[y], n);
		  
		  e = inner_prod(LeftOperator[x], RightOperator[y]);
	       }

	       // output the expectation value
	       if (Connected)
		  e -= ExpectOp1[x] * ExpectOp2[x];

	       if (Average)
	       {
		  Av.push_back(e);
	       }
	       else
	       {
		  std::cout << std::left << std::setw(5) << x << ' '
			    << std::left << std::setw(5) << (x+n) << ' ';
		  PrintFormat(e, ShowRealPart, ShowImagPart, ShowMagnitude, ShowArgument, ShowRadians);
		  std::cout << '\n';
		  std::cout << std::flush;
	       }
	    }
	 }

	 if (Average && Av.size() > 0)
	 {
	    std::complex<double> e = mean(Av);
	    double v = variance(Av);
	    std::cout << std::left << std::setw(5) << n << ' ';
	    PrintFormat(e, 
			ShowRealPart, ShowImagPart, ShowMagnitude, ShowArgument, ShowRadians);
	    std::cout << std::setw(10) << Av.size() << ' ' 
		      << std::setw(20) << std::sqrt(v) << '\n';
	    Av.clear();
	    std::cout << std::flush;
	 }
      }

      pheap::Shutdown();

   }
   catch (std::exception& e)
   {
      std::cerr << "Exception: " << e.what() << '\n';
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!\n";
      return 1;
   }
}
