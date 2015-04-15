// -*- C++ -*- $Id$

#include "unitcell.h"

LatticeSite
flip_conj(LatticeSite const& A)
{
   LatticeSite Result;

   if (A.empty())
      return Result;

   SiteBasis ReflectedBasis = adjoint(A.Basis1());
   for (LatticeSite::const_iterator ai = A.begin(); ai != A.end(); ++ai)
   {
      Result[ai->first] = flip_conj(ai->second, ReflectedBasis);
   }
   return Result;
}

//
// UnitCell members
//

UnitCell::UnitCell()
{
}

UnitCell::UnitCell(UnitCell const& Other)
   : Data_(Other.Data_), OperatorMap_(Other.OperatorMap_)
{
   for (operator_map_type::iterator I = OperatorMap_.begin(); I != OperatorMap_.end(); ++I)
   {
      I->second.SetUnitCell(this);
   }
}

UnitCell::UnitCell(LatticeSite const& s)
   : Data_(1, s)
{
   OperatorMap_["I"] = UnitCellMPO(this, this->identity_mpo(), LatticeCommute::Bosonic, 0);
}

UnitCell::UnitCell(LatticeSite const& s, LatticeSite const& t)
   : Data_(1, s)
{
   Data_.push_back(t);
   OperatorMap_["I"] = UnitCellMPO(this, this->identity_mpo(), LatticeCommute::Bosonic, 0);
}

UnitCell::UnitCell(LatticeSite const& s, LatticeSite const& t, LatticeSite const& u)
   : Data_(1, s)
{
   Data_.push_back(t);
   Data_.push_back(u);
   OperatorMap_["I"] = UnitCellMPO(this, this->identity_mpo(), LatticeCommute::Bosonic, 0);
}

UnitCell::UnitCell(LatticeSite const& s, LatticeSite const& t, LatticeSite const& u, LatticeSite const& v)
   : Data_(1, s)
{
   Data_.push_back(t);
   Data_.push_back(u);
   Data_.push_back(v);
   OperatorMap_["I"] = UnitCellMPO(this, this->identity_mpo(), LatticeCommute::Bosonic, 0);
}

UnitCell::UnitCell(SymmetryList const& sl, LatticeSite const& s)
   : Data_(1, CoerceSL(sl,s))
{
   OperatorMap_["I"] = UnitCellMPO(this, this->identity_mpo(), LatticeCommute::Bosonic, 0);
}

UnitCell::UnitCell(SymmetryList const& sl, LatticeSite const& s, LatticeSite const& t)
   : Data_(1, CoerceSL(sl, s))
{
   Data_.push_back(CoerceSL(sl, t));
   OperatorMap_["I"] = UnitCellMPO(this, this->identity_mpo(), LatticeCommute::Bosonic, 0);
}

UnitCell::UnitCell(SymmetryList const& sl, LatticeSite const& s, LatticeSite const& t, LatticeSite const& u)
   : Data_(1, CoerceSL(sl, s))
{
   Data_.push_back(CoerceSL(sl, t));
   Data_.push_back(CoerceSL(sl, u));
   OperatorMap_["I"] = UnitCellMPO(this, this->identity_mpo(), LatticeCommute::Bosonic, 0);
}

UnitCell::UnitCell(SymmetryList const& sl, LatticeSite const& s, LatticeSite const& t, 
                 LatticeSite const& u, LatticeSite const& v)
   : Data_(1, CoerceSL(sl, s))
{
   Data_.push_back(CoerceSL(sl, t));
   Data_.push_back(CoerceSL(sl, u));
   Data_.push_back(CoerceSL(sl, v));
   OperatorMap_["I"] = UnitCellMPO(this, this->identity_mpo(), LatticeCommute::Bosonic, 0);
}

UnitCell::UnitCell(int RepeatCount, UnitCell const& l)
{
   Data_.reserve(l.size()*RepeatCount);
   for (int i = 0; i < RepeatCount; ++i)
   {
      Data_.insert(Data_.end(), l.begin(), l.end());
   }
   OperatorMap_["I"] = UnitCellMPO(this, this->identity_mpo(), LatticeCommute::Bosonic, 0);
}

UnitCell::UnitCell(UnitCell const& x1, UnitCell const& x2)
   : Data_(x1.Data_)
{
   Data_.insert(Data_.end(), x2.begin(), x2.end());
   OperatorMap_["I"] = UnitCellMPO(this, this->identity_mpo(), LatticeCommute::Bosonic, 0);
}

UnitCell::UnitCell(UnitCell const& x1, UnitCell const& x2, UnitCell const& x3)
   : Data_(x1.Data_)
{
   Data_.insert(Data_.end(), x2.begin(), x2.end());
   Data_.insert(Data_.end(), x3.begin(), x3.end());
   OperatorMap_["I"] = UnitCellMPO(this, this->identity_mpo(), LatticeCommute::Bosonic, 0);
}

UnitCell::UnitCell(UnitCell const& x1, UnitCell const& x2, UnitCell const& x3, UnitCell const& x4)
   : Data_(x1.Data_)
{
   Data_.insert(Data_.end(), x2.begin(), x2.end());
   Data_.insert(Data_.end(), x3.begin(), x3.end());
   Data_.insert(Data_.end(), x4.begin(), x4.end());
   OperatorMap_["I"] = UnitCellMPO(this, this->identity_mpo(), LatticeCommute::Bosonic, 0);
}

UnitCell::UnitCell(int Size, LatticeSite const& s)
   : Data_(Size, s)
{
   OperatorMap_["I"] = UnitCellMPO(this, this->identity_mpo(), LatticeCommute::Bosonic, 0);
}

UnitCell&
UnitCell::operator=(UnitCell const& Other)
{
   Data_ = Other.Data_;
   OperatorMap_ = Other.OperatorMap_;
   for (operator_map_type::iterator I = OperatorMap_.begin(); I != OperatorMap_.end(); ++I)
   {
      I->second.SetUnitCell(this);
   }
   return *this;
}

LatticeSite const& 
UnitCell::operator[](int n) const
{
   return Data_[n];
}

bool
UnitCell::operator_exists(std::string const& s) const
{
   if (OperatorMap_.find(s) != OperatorMap_.end())
      return true;
   if (Data_.size() == 1 && Data_.front().operator_exists(s))
      return true;
   return false;
}

UnitCellMPO
UnitCell::Operator(std::string const& Op) const
{
   operator_map_type::const_iterator I = OperatorMap_.find(Op);
   if (I != OperatorMap_.end())
      return I->second;

   if (Data_.size() != 1)
   {
      PANIC("Operator not found in unit cell")(Op);
   }

   return this->LocalOperator(Op, 0);
}

UnitCellMPO
UnitCell::OperatorAtCell(std::string const& Op, int n) const
{
   operator_map_type::const_iterator I = OperatorMap_.find(Op);
   if (I != OperatorMap_.end())
      return UnitCellMPO(this, I->second.MPO(), I->second.Commute(), I->second.offset()+n);

   if (Data_.size() != 1)
   {
      PANIC("Operator not found in unit cell")(Op);
   }

   return this->LocalOperator(Op, n, 0);
}


UnitCellMPO&
UnitCell::Operator(std::string const& Op)
{
   return OperatorMap_[Op];
}

bool
UnitCell::operator_exists(std::string const& Op, int n) const
{
   if (n < 0 || n >= int(Data_.size()))
      return false;
   return Data_[n].operator_exists(Op);
}

LatticeCommute
UnitCell::Commute(std::string const& Op, int n) const
{
   CHECK(0 <= n && n < int(Data_.size()))("Site index is out of range")(n)(Data_.size());
   return Data_[n][Op].Commute();
}

UnitCellMPO
UnitCell::LocalOperator(std::string const& Op, int Cell, int n) const
{
   CHECK(0 <= n && n < int(Data_.size()))("Site index is out of range")(n)(Data_.size());

   SiteOperator Operator = Data_[n][Op];
   std::string SignOperator = Operator.Commute().SignOperator();

   FiniteMPO Result(Data_.size());

   BasisList Vacuum = make_vacuum_basis(Operator.GetSymmetryList());
   BasisList Basis = make_single_basis(Operator.TransformsAs());

   // Assemble the JW-string.
   for (int i = 0; i < n; ++i)
   {
      if (!Data_[i].operator_exists(SignOperator))
      {
	 WARNING("JW-string operator doesn't exist at a lattice site, using the identity")(i)(SignOperator);
      }
      SimpleOperator Op = Data_[i].operator_exists(SignOperator) ? Data_[i][SignOperator] : Data_[i]["I"];
      Result[i] = OperatorComponent(Op.Basis1(), Op.Basis2(), Basis, Basis);
      Result[i](0,0) = Op;
   }
   Result[n] = OperatorComponent(Operator.Basis1(), Operator.Basis2(), Basis, Vacuum);
   Result[n](0,0) = Operator;
   for (int i = n+1; i < this->size(); ++i)
   {
      SimpleOperator I = Data_[i]["I"];
      Result[i] = OperatorComponent(I.Basis1(), I.Basis2(), Vacuum, Vacuum);
      Result[i](0,0) = I;
   }

   return UnitCellMPO(this, Result, Operator.Commute(), Cell*this->size());
}

UnitCellMPO
UnitCell::LocalOperator(std::string const& Op, int n) const
{
   return this->LocalOperator(Op, 0, n);
}

UnitCellMPO
UnitCell::OperatorFunction(std::string const& Op, 
			   std::vector<std::complex<double> > const& Params) const
{
   PANIC("Operator function is not defined");
   return UnitCellMPO(this, FiniteMPO(), LatticeCommute::None, 0);
}

UnitCellMPO
UnitCell::OperatorFunction(std::string const& Op, int n,
			   std::vector<std::complex<double> > const& Params) const
{
   PANIC("Operator function is not defined");
   return UnitCellMPO(this, FiniteMPO(), LatticeCommute::None, 0);
}

UnitCellMPO
UnitCell::swap_gate(int i, int j) const
{
   return this->swap_gate(0,i,0,j);
}

UnitCellMPO
UnitCell::swap_gate(int Cell_i, int i, int Cell_j, int j) const
{
   // normal-order the sites
   if (Cell_i > Cell_j || (Cell_i == Cell_j && i > j))
   {
      std::swap(Cell_i, Cell_j);
      std::swap(i,j);
   }

   if (Cell_i == Cell_j && i == j)
   {
      return UnitCellMPO(this, this->identity_mpo(), LatticeCommute::Bosonic, Cell_i);
   }

   BasisList Basis_i = this->operator[](i).Basis1();
   BasisList Basis_j = this->operator[](j).Basis1();

   ProductBasis<BasisList, BasisList> Basis_ij(Basis_i, Basis_j);
   ProductBasis<BasisList, BasisList> Basis_ji(Basis_j, Basis_i);

   // The actual gate operator
   SimpleOperator Op = ::swap_gate(Basis_i, Basis_j, Basis_ji, Basis_ij);

   // decompose it back into sites   
   OperatorComponent R1, R2;
   boost::tie(R1, R2) = decompose_local_tensor_prod(Op, Basis_ji, Basis_ij);

   // now turn this into a FiniteMPO
   FiniteMPO Result(this->size() * (Cell_j-Cell_i+1));
   for (int n = 0; n <i; ++n)
   {
      Result[n] = OperatorComponent::make_identity(this->LocalBasis(n % this->size()));
   }
   Result[i] = R1;
   for (int n = i + 1; n < (Cell_j-Cell_i)*this->size() + j; ++n)
   {
      Result[n] = OperatorComponent::make_identity(this->LocalBasis(n % this->size()), R1.Basis2());
   }
   Result[(Cell_j-Cell_i)*this->size() + j] = R2;
   for (int n = (Cell_j-Cell_i)*this->size() + j + 1; n < (Cell_j-Cell_i+1)*this->size(); ++n)
   {
      Result[n] = OperatorComponent::make_identity(this->LocalBasis(n % this->size()));
   }

   Result.debug_check_structure();
   return UnitCellMPO(this, Result, LatticeCommute::Bosonic, Cell_i);
}

FiniteMPO
UnitCell::identity_mpo(QuantumNumbers::QuantumNumber const& q) const
{
   FiniteMPO Result(this->size());
   BasisList b = make_single_basis(q);
   for (int i = 0; i < this->size(); ++i)
   {
      Result[i] = OperatorComponent(this->LocalBasis(i), b, b);
      Result[i](0,0) = this->operator[](i)["I"];
   }
   return Result;
}

FiniteMPO
UnitCell::identity_mpo() const
{
   return this->identity_mpo(QuantumNumbers::QuantumNumber(this->GetSymmetryList()));
}

FiniteMPO
UnitCell::string_mpo(std::string const& OpName, QuantumNumbers::QuantumNumber const& Trans) const
{
   FiniteMPO Result(this->size());

   BasisList Vacuum = make_single_basis(Trans);

   // Assemble the JW-string
   for (int i = 0; i < this->size(); ++i)
   {
      if (!this->operator_exists(OpName))
      {
	 WARNING("JW-string operator doesn't exist at a lattice site, using the identity")(i)(OpName);
      }
      SimpleOperator Op = Data_[i].operator_exists(OpName) ? Data_[i][OpName] : this->operator[](i)["I"];
      Result[i] = OperatorComponent(Op.Basis1(), Op.Basis2(), Vacuum, Vacuum);
      Result[i](0,0) = Op;
   }
   return Result;
}

FiniteMPO
UnitCell::string_mpo(LatticeCommute Com, QuantumNumbers::QuantumNumber const& Trans) const
{
   return this->string_mpo(Com.SignOperator(), Trans);
}

PStream::opstream& operator<<(PStream::opstream& out, UnitCell const& L)
{
   out << L.Data_;
   out << L.OperatorMap_;
   return out;
}

PStream::ipstream& operator>>(PStream::ipstream& in, UnitCell& L)
{
   in >>  L.Data_;
   in >>  L.OperatorMap_;
   // set the UnitCell member of each UnitCellOperator
   for (UnitCell::operator_map_type::iterator I = L.OperatorMap_.begin(); I != L.OperatorMap_.end(); ++I)
   {
      I->second.SetUnitCell(&L);
   }
   return in;
}

UnitCell repeat(UnitCell const& x, int RepeatCount)
{
   return UnitCell(RepeatCount, x);
}

UnitCell join(UnitCell const& x, UnitCell const& y)
{
   return UnitCell(x, y);
}

UnitCell join(UnitCell const& x, UnitCell const& y, UnitCell const& z)
{
   return UnitCell(x,y,z);
}

UnitCell join(UnitCell const& x, UnitCell const& y, UnitCell const& z, UnitCell const& w)
{
   return UnitCell(x,y,z,w);
}

UnitCell join(UnitCell const& x, UnitCell const& y, UnitCell const& z, UnitCell const& w,
	     UnitCell const& v)
{
   return join(UnitCell(x,y,z,w), v);
}

#if 0
FiniteMPO
UnitCell::Parse(std::string const& s)
{
}
#endif