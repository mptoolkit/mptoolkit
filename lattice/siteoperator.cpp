// -*- C++ -*- $Id$

#include "siteoperator.h"
#include <boost/lexical_cast.hpp>

// SiteBasis

SiteBasis::SiteBasis(SymmetryList const& SList)
  : Basis_(SList), Label_(new std::vector<std::string>())
{
}

SiteBasis::SiteBasis(std::string const& S)
  : Basis_(SymmetryList(S)), Label_(new std::vector<std::string>())
{
}

void
SiteBasis::push_back(std::string const& Label, QuantumNumber const& q)
{
   Basis_.push_back(Coerce(q, this->GetSymmetryList()));
   Label_.mutate()->push_back(Label);
}

int
SiteBasis::Lookup(std::string const& Label) const
{
   int Pos = std::find(Label_->begin(), 
		       Label_->end(), Label) - Label_->begin();
   if (Pos == int(Label_->size())) { PANIC("Invalid site basis element")(Label); }
   return Pos;
}

int
SiteBasis::LookupOrNeg(std::string const& Label) const
{
   int Pos = std::find(Label_->begin(), 
		       Label_->end(), Label) - Label_->begin();
   if (Pos == int(Label_->size())) return -1;
   return Pos;
}

std::ostream& operator<<(std::ostream& out, SiteBasis const& Basis)
{
   for (std::size_t i = 0; i < Basis.size(); ++i)
   {
      if (i != 0) out << ", ";
      out << "{ " << i << ", |" << Basis.Label(i) << ">, " 
          << Basis.qn(i) << " }";
   }
   return out;
}

void show_projections(std::ostream& out, SiteBasis const& Basis)
{
   for (std::size_t i = 0; i < Basis.size(); ++i)
   {
      if (i != 0) out << ", ";

      out << "{ ";
      std::vector<QuantumNumbers::Projection> Projections;
      enumerate_projections(Basis.qn(i), std::back_inserter(Projections));
      for (size_type j = 0; j < Projections.size(); ++j)
      {
        if (j != 0) out << ", ";
	out << '|' << Basis.Label(i) << ' ' << Basis.qn(i)
	    << ", " << Projections[j]
	    << ">";
      }
      out << " }";
   }
}

PStream::opstream& operator<<(PStream::opstream& out, SiteBasis const& B)
{
   return out << B.Basis_ << B.Label_;
}

PStream::ipstream& operator>>(PStream::ipstream& in, SiteBasis& B)
{
   return in >> B.Basis_ >> B.Label_;
}

// SiteProductBasis

SiteProductBasis::SiteProductBasis(SiteBasis const& B1, SiteBasis const& B2)
  : Basis_(B1.GetSymmetryList()), 
    Basis1_(B1), Basis2_(B2), ProductBasis_(B1, B2)
{
   for (std::size_t i = 0; i < ProductBasis_.size(); ++i)
   {
      ProductBasis<BasisList, BasisList>::source_type s = ProductBasis_.rmap(i);
      std::string Label = B1.Label(s.first) + " x " 
	 + B2.Label(s.second) + " [" 
	 + boost::lexical_cast<std::string>(ProductBasis_[i]) + ']';
      Basis_.push_back(Label, ProductBasis_[i]);
   }
}

std::ostream& operator<<(std::ostream& out, SiteProductBasis const& Basis)
{
   for (size_type i = 0; i < Basis.size(); ++i)
   {
      if (i != 0) out << ", ";
      out << "{ " << i << " = (" << Basis.PBasis().rmap(i).first
	  << "x" << Basis.PBasis().rmap(i).second << "), |" << Basis.Basis().Label(i) << ">, " 
          << Basis.PBasis()[i] << " }";
   }
   return out;
}


// SiteOperator

PStream::opstream& operator<<(PStream::opstream& out, SiteOperator const& Op)
{
   return out << Op.base() << Op.Basis_ << Op.Com_;
}

PStream::ipstream& operator>>(PStream::ipstream& in, SiteOperator& Op)
{
   in >> Op.base() >> Op.Basis_ >> Op.Com_;
   return in;
}

std::ostream& operator<<(std::ostream& out, SiteOperator const& Op)
{
   out << '[' << Op.TransformsAs() << "] ";
   bool first = true;
   for (std::size_t i = 0; i < Op.size1(); ++i)
   {
      for (std::size_t j = 0; j < Op.size2(); ++j)
      {
         if (!is_transform_target(Op.Basis()[j].second, Op.TransformsAs(), Op.Basis()[i].second)) 
            continue;

	 std::complex<double> x = Op(i,j);
	 if (x == std::complex<double>()) continue;

	 if (x.imag() != 0)
	 {
	    if (!first) out << " + ";
	    out << '(' << x.real();
	    if (x.imag() >= 0) out << '+';
	    out << x.imag() << "I) ";
	 }
	 else // x.imag() == 0
	 {
	    if (!first)
	    {
	       if (x.real() >= 0) out << " + ";
	       else out << " - ";
	    }
	    else if (x.real() < 0) out << '-';
	    if (std::abs(x.real()) != 1) out << std::abs(x.real()) << ' ';
	 }
	 first = false;

	 out << "|" << Op.Basis().Label(i) << "><" << Op.Basis().Label(j) << "|";
      }
   }
   return out;
}

void show_projections(std::ostream& out, SiteOperator const& Op)
{
   out << "transforms as " << Op.TransformsAs() << '\n';

   SiteBasis Basis(Op.Basis());

   std::vector<QuantumNumbers::Projection> Projections;
   enumerate_projections(Op.TransformsAs(), std::back_inserter(Projections));
   for (size_type km = 0; km < Projections.size(); ++km)
   {
      out << "Projection " << std::setw(10) << Projections[km] << " :\n";

      for (std::size_t i = 0; i < Basis.size(); ++i)
      {
	 QuantumNumber qi = Basis.qn(i);
	 for (std::size_t j = 0; j < Basis.size(); ++j)
	 {
            if (!is_transform_target(Op.Basis()[j].second, Op.TransformsAs(), Op.Basis()[i].second)) 
               continue;

	    std::complex<double> x(Op(i,j));
	    if (x == 0.0) continue;

	    //	    out << x << std::endl;

	    QuantumNumber qj = Basis.qn(j);

	    std::vector<QuantumNumbers::Projection> mi, mj;
	    enumerate_projections(qi, std::back_inserter(mi));
	    enumerate_projections(qj, std::back_inserter(mj));

	    for (std::size_t ii = 0; ii < mi.size(); ++ii)
	    {
	       for (std::size_t jj = 0; jj < mj.size(); ++jj)
	       {
		  std::complex<double> y = x * clebsch_gordan(qj,    Op.TransformsAs(),  qi,
							      mj[jj], Projections[km], mi[ii]);

		  if (numerics::norm_2(y) > 1E-10)
		  {
		     out << "   " << y
			 << " |" << Basis.Label(i) << ' ' << qi
			 << ", " << mi[ii]
			 << "> <"
			 << Basis.Label(j) << ' '
			 << qj
			 << ", " << mj[jj]
			 << "|\n";
		  }
	       }
	    }
	 }
      }
      //      out << '\n';
   }
}

#if 0
std::vector<std::pair<SiteOperator, SiteOperator> >
decompose_tensor_prod(SiteOperator const& S, SiteProductBasis const& SPBasis)
{
   SiteBasis BLeft = SPBasis.Basis1();
   SiteBasis BRight = SPBasis.Basis2();

   ProductBasis<BasisList, BasisList> ab = SPBasis.PBasis();
   ProductBasis<BasisList, BasisList> alpha(BLeft.Basis(), adjoint(BLeft.Basis()));
   ProductBasis<BasisList, BasisList> beta(BRight.Basis(), adjoint(BRight.Basis()));
   std::vector<std::pair<SiteOperator, SiteOperator> > Result;

   IrredTensor<std::complex<double> > Partial = partial_transpose(S.base(), ab, ab, alpha, beta);
   LinearAlgebra::Matrix<std::complex<double> > M(Partial.data());
   LinearAlgebra::Matrix<std::complex<double> > U, Vt;
   LinearAlgebra::Vector<double> D;
   LinearAlgebra::SingularValueDecomposition(M, U, D, Vt);

   for (std::size_t i = 0; i < size(D); ++i)
   {
      if (D[i] < 1E-10) continue;

      SiteOperator A;
      for (std::size_t k = 0; k < size1(U); ++k)
      {
         if (LinearAlgebra::norm_2(U(k,i)) > 1E-10)
         {
            if (A.is_null())
               A = SiteOperator(BLeft, alpha[k]);
            else
               CHECK_EQUAL(alpha[k], A.TransformsAs());

            int l,m;
            boost::tie(l,m) = alpha.rmap(k);
            A(l,m) = U(k,i) * D[i];
         }
      }

      SiteOperator B;
      for (std::size_t k = 0; k < size2(Vt); ++k)
      {
         if (LinearAlgebra::norm_2(Vt(i,k)) > 1E-10)
         {
            if (B.is_null())
               B = SiteOperator(BRight, beta[k]);
            else
               CHECK_EQUAL(beta[k], B.TransformsAs());

            int l,m;
            boost::tie(l,m) = beta.rmap(k);
            B(l,m) = Vt(i,k);
         }
      }
      Result.push_back(std::make_pair(A,B));
   }
   return Result;
}
#endif

SiteOperator flip_conj(SiteOperator const& s, SiteBasis const& ReflectedBasis)
{
   DEBUG_CHECK_EQUAL(ReflectedBasis, adjoint(s.Basis()));
   SiteOperator Result(ReflectedBasis, adjoint(s.TransformsAs()), s.Commute());
   Result.data() = s.data();
   return Result;
}