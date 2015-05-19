// -*- C++ -*- $Id$

#include "random_wavefunc.h"
#include "linearalgebra/matrix_utility.h"

bool
WavefunctionDesc::Flip(std::vector<BasisList> const& Basis, int Site, int NewState, 
                       QuantumNumber const& NewQuantumNumber)
{ 
   Projection Delta = difference(NewQuantumNumber, Height[Site]);
   //   TRACE(Delta);
   State[Site] = NewState;
   for (int i = Site; i >= 0; --i)
   {
      if (!is_possible(Height[i], Delta)) return false;
      //      TRACE(i)(Height[i])(change(Height[i], Delta));

      QuantumNumbers::QuantumNumberList N = transform_targets(Height[i], Basis[i][State[i]]);
      if (N.size() < 1) return false;

      Height[i] = change(Height[i], Delta);

      if (!is_transform_target(Height[i+1], Basis[i][State[i]], Height[i])) 
      {
         //DEBUG_TRACE(Height[i+1])(Basis[i][State[i]])(Height[i]);
         return false;
      }
   }
   return true;
}

WavefunctionDesc::WavefunctionDesc(std::vector<BasisList> const& L) 
   : State(L.size()), Height(L.size()+1)
{
   QuantumNumber Q(L[0].GetSymmetryList());
   Height[L.size()] = Q;
   for (int i = L.size()-1; i >= 0; --i)
   {
      BasisList Basis = L[i];
      State[i] = rand() % Basis.size();

      QuantumNumbers::QuantumNumberList QList = transform_targets(Q, Basis[State[i]]);
      Q = QList[rand() % QList.size()];
      Height[i] = Q;
   }
}

WavefunctionDesc::WavefunctionDesc(std::vector<BasisList> const& L,
                                   QuantumNumber Q) 
   : State(L.size()), Height(L.size()+1)
{
   Height[L.size()] = Q;
   for (int i = L.size()-1; i >= 0; --i)
   {
      BasisList Basis = L[i];
      State[i] = rand() % Basis.size();

      QuantumNumbers::QuantumNumberList QList = transform_targets(Q, Basis[State[i]]);
      Q = QList[rand() % QList.size()];
      Height[i] = Q;
   }
}

std::ostream& operator<<(std::ostream& out, WavefunctionDesc const& Config)
{
   int i = Config.State.size();
   out << "Height: " << Config.Height[i] << '\n';
   while (i > 0)
   {
      --i;
      out << "State: " << Config.State[i] << '\n';
      out << "Height: " << Config.Height[i] << '\n';
   }
   return out;
}

void WavefunctionDesc::CheckValid(std::vector<BasisList> const& Basis) const
{
   for (int i = Basis.size()-1; i >= 0; --i)
   {
      CHECK(is_projection(Basis[i][State[i]], difference(Height[i], Height[i+1])));
   }
}

WavefunctionDesc
CreateRandomConfiguration(std::vector<BasisList> const& Basis, 
                          QuantumNumber const& q, double Beta,
                          QuantumNumber RightBoundary)
{
   if (RightBoundary.is_null())
      RightBoundary = QuantumNumber(Basis[0].GetSymmetryList());
   double now = clock() / CLOCKS_PER_SEC;
   WavefunctionDesc Psi(Basis, RightBoundary);
   while (q != Psi.TransformsAs())
   {
      double c = weight(difference(q, Psi.TransformsAs()));
      int Site = rand() % Basis.size();
      int NewState = rand() % Basis[Site].size();
      QuantumNumbers::QuantumNumberList QList = transform_targets(Psi.Height[Site+1], Basis[Site][NewState]);
      //      for (std::size_t Qi = 0; Qi < QList.size(); ++Qi)
      {
         std::size_t Qi = std::size_t(LinearAlgebra::random<double>() * QList.size());

	 QuantumNumber NewQ(QList[Qi]);
	 WavefunctionDesc New = Psi;
	 //	 TRACE(Site)(c)(NewState)(Psi.Height[Site+1])(Psi.Height[Site])(Psi.State[Site])(NewQ)(Basis[Site].qn(Psi.State[Site]))(Basis[Site].qn(NewState))(QList);
	 if (New.Flip(Basis, Site, NewState, NewQ))
	 {
	    //	    TRACE(Psi.TransformsAs())(New.TransformsAs());
	    if (clock() / CLOCKS_PER_SEC > now+10)
	    {
	       now = clock() / CLOCKS_PER_SEC;
	       TRACE("slow convergence!")(Site)(q)(NewQ)(New.TransformsAs())(Psi.TransformsAs());
	    }
	    double Newc = weight(difference(q, New.TransformsAs()));
            //TRACE(Newc);
	    if (LinearAlgebra::random<double>() < exp(Beta * (c-Newc))) 
	    {
               //TRACE("Accepted");
	       Psi = New;
               //	       continue;
	    }
	 }
      }
   }
   return Psi;
}

LinearWavefunction CreateRandomWavefunction(std::vector<BasisList> const& Basis, 
					QuantumNumber const& q, double Beta,
                                        QuantumNumber const& RightBoundary)
{
   WavefunctionDesc Psi = CreateRandomConfiguration(Basis, q, Beta, RightBoundary);

   BasisList Vac = make_single_basis(RightBoundary);
   VectorBasis B2(Vac);
   QuantumNumber Ident(Basis[0].GetSymmetryList());  // the scalar quantum number

   LinearWavefunction Result;

   MatrixOperator Center(B2, B2, Ident);
   Center(0,0) = LinearAlgebra::Matrix<double>(1,1,1);

   for (int i = Basis.size()-1; i >= 0; --i)
   {
      //      TRACE(i)(Psi.Height[i]);
      VectorBasis B1(Basis[i].GetSymmetryList());
      B1.push_back(Psi.Height[i], 1);
      MatrixOperator Next(B1, B2, Basis[i][Psi.State[i]]);
      Next(0,0) = LinearAlgebra::Matrix<double>(1,1,1);
      StateComponent R(Basis[i], B1, B2);
      R[Psi.State[i]] = Next;
      R = prod(R, Center);
      Center = TruncateBasis1(R);
      Result.push_front(R);
      B2 = B1;
   }

   return Result;
}

LinearWavefunction CreateRandomWavefunction(std::vector<BasisList> const& Basis, 
					QuantumNumber const& q, double Beta)
{
   return CreateRandomWavefunction(Basis, q, Beta, QuantumNumber(Basis[0].GetSymmetryList()));
}

#if 0

LinearWavefunction CreateRandomWavefunction(Lattice const& L, QuantumNumber const& q, double Beta,
                                        QuantumNumber const& RightBoundary)
{
   std::vector<BasisList> Basis(L.size());
   //   Lattice::const_iterator Li = L.end();
   Lattice::const_iterator Li = L.begin();
   for (int i = 0; i < L.size(); ++i)
   {
      //      --Li;
      Basis[i] = Li->Basis1().Basis();
      ++Li;
   }
   return CreateRandomWavefunction(Basis, q, Beta, RightBoundary);
}

LinearWavefunction CreateRandomWavefunction(Lattice const& L, QuantumNumber const& q, double Beta)
{
   return CreateRandomWavefunction(L, q, Beta, QuantumNumber(L.GetSymmetryList()));
}

LinearWavefunction CreateRandomWavefunction(Lattice const& L, QuantumNumber const& q, 
                                        double Beta,
                                        QuantumNumber const& RightBoundary,
                                        int Count)
{
   LinearWavefunction Psi = CreateRandomWavefunction(L, q, Beta, RightBoundary);
   while (Count > 1)
   {
      std::cout << "Working..." << std::endl;
      LinearWavefunction P2 = CreateRandomWavefunction(L, q, Beta, RightBoundary);
      //      P2 *= 2.0 * (double(rand()) / RAND_MAX) - 1.0;
      //      Psi = Psi + P2;
      --Count;
   }
   //   Psi.normalize();
   return Psi;
}

LinearWavefunction CreateRandomWavefunction(Lattice const& L, QuantumNumber const& q, 
                                        double Beta,
                                        int Count)
{
   return CreateRandomWavefunction(L, q, Beta, QuantumNumber(L.GetSymmetryList()), Count);
}

#endif

std::complex<double>
Amplitude(LinearWavefunction const& Psi, WavefunctionDesc const& Config)
{
   MatrixOperator x = MatrixOperator::make_identity(Psi.Basis2());
   LinearWavefunction::const_iterator I = Psi.end();
   int s = Psi.size();
   while (I != Psi.begin())
   {
      --I; --s;
      x = prod((*I)[Config.State[s]], x, Config.Height[s]);
      //DEBUG_CHECK(norm_frob_sq(x) > 0)(x);
   }
   CHECK_EQUAL(s, 0);
   if (iterate_at(x.data(), 0, 0))
      return x(0,0)(0,0);
   else return 0.0;
}
