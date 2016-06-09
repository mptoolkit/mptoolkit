TriangularMPO TriangularOneSite(SimpleOperator const& x);

// A one-site operator with the given momentum, in angular units
TriangularMPO TriangularOneSite(SimpleOperator const& x, double Momentum);

TriangularMPO TriangularTwoSite(SimpleOperator const& x, SimpleOperator const& y,
                                     QuantumNumbers::QuantumNumber const& Trans);

TriangularMPO TriangularTwoSite(SimpleOperator const& x, SimpleOperator const& y);

TriangularMPO TriangularTwoSiteExponential(SimpleOperator const& x, SimpleOperator const& y, 
                                                std::complex<double> Factor, QuantumNumber const& Trans);

TriangularMPO TriangularTwoSiteExponential(SimpleOperator const& x, SimpleOperator const& y, 
                                                std::complex<double> Factor);

TriangularMPO TriangularThreeSite(SimpleOperator const& x, 
                                       SimpleOperator const& y, 
                                       SimpleOperator const& z);

// a two point interaction between sites n1 and n2 of a lattice.
// This probably should use SiteOperator's rather than SimpleOperator's, so that
// we can handle fermions, especially for the case n2<n1.
// this is assumed to be a bosonic operator; if n2 < n1 we swap the sites
// (which is incorrect for fermions)
TriangularMPO TwoPointOperator(std::vector<BasisList> const& Sites, 
                                    int n1, SimpleOperator const& x1,
                                    int n2, SimpleOperator const& x2);

// This is a restricted implementation; we must have
// PRECONDITION: 0 <= n1 < Sites.size()
// PRECONDITION: Sites.size() <= n2 < 2*Sites.size()
// That is, the operator1 must be in the first unit cell, the operator2 must be
// in the second unit cell.
TriangularMPO TwoPointExponentialOperator(std::vector<BasisList> const& Sites, 
                                          int n1, SimpleOperator const& x1,
                                          int n2, SimpleOperator const& x2,
                                          std::complex<double> Factor);


// a two-point string operator where the String term is inserted
// at sites n1+1, n1+2, ..., n2-1.
// Because this function can be used to implement fermionic operators,
// we demand normal ordering of sites; it is an error to call this function
// with n2<n1.
TriangularMPO TwoPointStringOperator(std::vector<BasisList> const& Sites, 
                                     int n1, SimpleOperator const& x1,
                                     SimpleOperator const& String,
                                     int n2, SimpleOperator const& x2);

TriangularMPO TwoPointStringOperator(std::vector<BasisList> const& Sites, 
                                     int n1, SimpleOperator const& x1,
                                     std::vector<SimpleOperator> const& String,
                                     int n2, SimpleOperator const& x2);

// A one-site operator on a lattice with a given momentum, in angular units per unit cell
TriangularMPO OnePointOperator(std::vector<BasisList> const& Sites, 
                                    int n, SimpleOperator const& x, double Momentum = 0);

// A one-site operator on a lattice with a given momentum, in angular units per unit cell
TriangularMPO OnePointStringOperator(std::vector<BasisList> const& Sites, 
					  std::vector<SimpleOperator> const& String,
					  int n, SimpleOperator const& x, double Momentum = 0);

// A one-site operator on a lattice with a given momentum, in angular units per unit cell
// Here we exclude the string term on the unit cell where the operator is defined
TriangularMPO OnePointStringOperatorExclude(std::vector<BasisList> const& Sites, 
					    std::vector<SimpleOperator> const& String,
					    int n, SimpleOperator const& x, double Momentum = 0);

TriangularMPO OneCellStringOperator(std::vector<BasisList> const& Sites, 
				    std::vector<SimpleOperator> const& String,
				    std::vector<SimpleOperator> const& CellOp,
				    double Momentum = 0);


TriangularMPO TriangularOneSite(SimpleOperator const& x)
{
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());
   BasisList b(x.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(x.TransformsAs());

   SimpleOperator I = SimpleOperator::make_identity(x.Basis1());

   OperatorComponent data(x.Basis1(), b, b);

   data.set_operator(0,0,I);
   data.set_operator(1,1,I);
   data.set_operator(0,1,x);

   return TriangularMPO(data);
}

TriangularMPO TriangularOneSite(SimpleOperator const& x, double Momentum)
{
   std::complex<double> r = std::exp(std::complex<double>(0.0, 1.0) * Momentum);

   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());
   BasisList b(x.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(x.TransformsAs());

   SimpleOperator I = SimpleOperator::make_identity(x.Basis1());

   OperatorComponent data(x.Basis1(), b, b);
   data.set_operator(0,0,r*I);
   data.set_operator(1,1,I);
   data.set_operator(0,1,x);

   return TriangularMPO(data);
}

TriangularMPO TriangularTwoSite(SimpleOperator const& x, SimpleOperator const& y, QuantumNumber const& Trans)
{
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());
   BasisList b(x.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(y.TransformsAs());
   b.push_back(Trans);

   SimpleOperator I = SimpleOperator::make_identity(x.Basis1());

   OperatorComponent data(x.Basis1(), b, b);
   data.set_operator(0,0,I);
   data.set_operator(2,2,I);
   data.set_operator(0,1,x);
   data.set_operator(1,2,y);

   return TriangularMPO(data);
}

TriangularMPO TriangularTwoSite(SimpleOperator const& x, SimpleOperator const& y)
{
   QuantumNumber q;
   if (num_transform_targets(x.TransformsAs(), y.TransformsAs()) == 1)
      q = transform_targets(x.TransformsAs(), y.TransformsAs())[0];
   else
      q = QuantumNumber(x.GetSymmetryList());
   return TriangularTwoSite(x, y, q);
}

TriangularMPO TriangularTwoSiteExponential(SimpleOperator const& x, SimpleOperator const& y, std::complex<double> Factor, QuantumNumber const& Trans)
{
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());
   BasisList b(x.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(y.TransformsAs());
   b.push_back(Trans);

   SimpleOperator I = SimpleOperator::make_identity(x.Basis1());

   OperatorComponent data(x.Basis1(), b, b);
   data.set_operator(0,0,I);
   data.set_operator(2,2,I);
   data.set_operator(1,1,Factor*I);
   data.set_operator(0,1,x);
   data.set_operator(1,2,y);

   return TriangularMPO(data);
}

TriangularMPO TriangularTwoSiteExponential(SimpleOperator const& x, SimpleOperator const& y, std::complex<double> Factor)
{
   QuantumNumber q;
   if (num_transform_targets(x.TransformsAs(), y.TransformsAs()) == 1)
      q = transform_targets(x.TransformsAs(), y.TransformsAs())[0];
   else
      q = QuantumNumber(x.GetSymmetryList());
   return TriangularTwoSiteExponential(x, y, Factor, q);
}



TriangularMPO TriangularThreeSite(SimpleOperator const& x, SimpleOperator const& y, SimpleOperator const& z, 
				   QuantumNumber const& yz_trans)
{
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());

   BasisList b(x.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(yz_trans);
   b.push_back(z.TransformsAs());
   b.push_back(Ident);

   SimpleOperator I = SimpleOperator::make_identity(x.Basis1());

   OperatorComponent data(x.Basis1(), b, b);
   data.set_operator(0,0,I);
   data.set_operator(3,3,I);
   data.set_operator(0,1,x);
   data.set_operator(1,2,y);
   data.set_operator(2,3,z);

   return TriangularMPO(data);
}

TriangularMPO TriangularStretchedTwoSite(SimpleOperator const& x, int NumNeighbors,
					  SimpleOperator const& y)
{
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());
   BasisList b(x.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(y.TransformsAs());
   for (int i = 0; i < NumNeighbors; ++i)
   {
      b.push_back(y.TransformsAs());
   }
   b.push_back(Ident);

   SimpleOperator I = SimpleOperator::make_identity(x.Basis1());

   OperatorComponent data(x.Basis1(), b, b);
   data.set_operator(0,0,I);
   data.set_operator(0,1,x);
   for (int i = 0; i < NumNeighbors; ++i)
   {
      data.set_operator(i+1, i+2, I);
   }
   data.set_operator(1+NumNeighbors,2+NumNeighbors,y);
   data.set_operator(2+NumNeighbors,2+NumNeighbors,I);

   return TriangularMPO(data);
}


TriangularMPO TriangularThreeSite(SimpleOperator const& x, SimpleOperator const& y, SimpleOperator const& z)
{
   QuantumNumbers::QuantumNumberList ql = transform_targets(y.TransformsAs(), z.TransformsAs());
   CHECK_EQUAL(ql.size(), 1)("ambiguous coupling of y and z")(y.TransformsAs())(z.TransformsAs());
   return TriangularThreeSite(x,y,z, ql[0]);
}

TriangularMPO TriangularFourSite(SimpleOperator const& w, SimpleOperator const& x, QuantumNumber const& wx_trans,
                                 SimpleOperator const& y, QuantumNumber const& wxy_trans,
                                 SimpleOperator const& z)
{
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());

   BasisList b(w.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(adjoint(w.TransformsAs()));
   b.push_back(adjoint(wx_trans));
   b.push_back(adjoint(wxy_trans));
   b.push_back(Ident);

   SimpleOperator I = SimpleOperator::make_identity(x.Basis1());

   OperatorComponent data(x.Basis1(), b, b);
   data.set_operator(0,0,I);
   data.set_operator(4,4,I);
   data.set_operator(0,1,w);
   data.set_operator(1,2,x);
   data.set_operator(2,3,y);
   data.set_operator(3,4,z);

   return TriangularMPO(data);
}

TriangularMPO TriangularFourSite(SimpleOperator const& w, SimpleOperator const& x, 
                                  SimpleOperator const& y, SimpleOperator const& z)
{
   QuantumNumbers::QuantumNumberList ql = transform_targets(w.TransformsAs(), x.TransformsAs());
   CHECK_EQUAL(ql.size(), 1)("ambiguous coupling of w and x")(w.TransformsAs())(x.TransformsAs());
   QuantumNumber wx = ql[0];

   ql = transform_targets(wx, y.TransformsAs());
   CHECK_EQUAL(ql.size(), 1)("ambiguous coupling of wx and y")(y.TransformsAs())
      (y.TransformsAs())(z.TransformsAs());
   QuantumNumber wxy = ql[0];

   return TriangularFourSite(w, x, wx, y, wxy, z);
}

#if 0
TriangularMPO ZigZagChain(SimpleOperator const& S, SimpleOperator const& T, double J1, double J2)
{
   QuantumNumbers::QuantumNumber Ident(S.GetSymmetryList());

   BasisList b(S.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(S.TransformsAs());
   b.push_back(S.TransformsAs());
   b.push_back(Ident);

   SimpleOperator I = SimpleOperator::make_identity(S.Basis1());

   OperatorComponent data(S.Basis1(), b, b);
   data.set_operator(0,0,I);
   data.set_operator(3,3,I);
   data.set_operator(1,0,S);
   data.set_operator(2,1,I);
   data.set_operator(3,2,J2*T);
   data.set_operator(3,1,J1*T);

   return TriangularMPO(data);
}

// temporary hacks for specific hamiltonians

TriangularMPO TriangularTwoSitePBC(SimpleOperator const& x, SimpleOperator const& y, QuantumNumber const& Trans)
{
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());
   BasisList b(x.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(y.TransformsAs());
   b.push_back(y.TransformsAs());
   b.push_back(Trans);

   SimpleOperator I = SimpleOperator::make_identity(x.Basis1());

   OperatorComponent data(x.Basis1(), b, b);
   data.set_operator(0,0,I);
   data.set_operator(2,2,I);
   data.set_operator(3,3,I);
   data.set_operator(1,0,y);
   data.set_operator(3,1,x);

   return TriangularMPO(data);
}

TriangularMPO TriangularTwoSitePBC(SimpleOperator const& x, SimpleOperator const& y)
{
   CHECK_EQUAL(num_transform_targets(x.TransformsAs(), y.TransformsAs()), 1);
   return TriangularTwoSitePBC(x, y, transform_targets(x.TransformsAs(), y.TransformsAs())[0]);
}

TriangularMPO TriangularTwoSitePBC_Boundary(SimpleOperator const& x, SimpleOperator const& y, QuantumNumber const& Trans)
{
   QuantumNumbers::QuantumNumber Ident(x.GetSymmetryList());
   BasisList b(x.GetSymmetryList());
   b.push_back(Ident);
   b.push_back(y.TransformsAs());
   b.push_back(y.TransformsAs());
   b.push_back(Trans);

   SimpleOperator I = SimpleOperator::make_identity(x.Basis1());

   OperatorComponent data(x.Basis1(), b, b);
   data.set_operator(0,0,I);
   data.set_operator(3,3,I);
   data.set_operator(1,0,y);
   data.set_operator(2,0,y);
   data.set_operator(3,1,x);
   data.set_operator(3,2,x);

   return TriangularMPO(data);
}

TriangularMPO TriangularTwoSitePBC_Boundary(SimpleOperator const& x, SimpleOperator const& y)
{
   CHECK_EQUAL(num_transform_targets(x.TransformsAs(), y.TransformsAs()), 1);
   return TriangularTwoSitePBC_Boundary(x, y, transform_targets(x.TransformsAs(), y.TransformsAs())[0]);
}
#endif

TriangularMPO OnePointStringOperator(std::vector<BasisList> const& Sites, 
					  std::vector<SimpleOperator> const& String,
					  int n, SimpleOperator const& x, double Momentum)
{
   CHECK(n >= 0 && n < int(Sites.size()))(n)(Sites.size());
   CHECK_EQUAL(Sites.size(), String.size());
   CHECK_EQUAL(x.Basis1(), Sites[n]);
   CHECK_EQUAL(x.Basis2(), Sites[n]);
   int const Size = Sites.size();

   // The way we handle momentum here is probably incorrect, we should do the entire momentum
   // on site 0
   std::complex<double> r = std::exp(std::complex<double>(0.0, 1.0) * (Momentum / Size));

   // construct a list of added quantum numbers that are needed for the bond bases
   std::vector<BasisList> BondBasis(Size, BasisList(Sites[0].GetSymmetryList()));

   // Add the identity component to the bond basis
   for (int i = 0; i < Size; ++i)
   {
      BondBasis[i].push_back(QuantumNumber(Sites[0].GetSymmetryList()));
   }

   // the actual operator
   for (int i = 0; i < Size; ++i)
   {
      BondBasis[i].push_back(x.TransformsAs());
   }

   // Assemble the operators
   GenericMPO Result(Size);
   for (int i = 0; i < Size; ++i)
   {
      CHECK_EQUAL(String[i].Basis1(), Sites[i]);
      CHECK_EQUAL(String[i].Basis2(), Sites[i]);

      Result[i] = OperatorComponent(Sites[i], Sites[i], BondBasis[i], BondBasis[smod(i+1,Size)]);
      Result[i](0,0) = i == 0 ? (r * String[i]) : r * String[i];
      Result[i](1,1) = SimpleOperator::make_identity(Sites[i]);
   }

   Result[n](0,1) = x;
   return TriangularMPO(Result.data());
}

TriangularMPO OnePointStringOperatorExclude(std::vector<BasisList> const& Sites, 
					    std::vector<SimpleOperator> const& String,
					    int n, SimpleOperator const& x, double Momentum)
{
   CHECK(n >= 0 && n < int(Sites.size()))(n)(Sites.size());
   CHECK_EQUAL(Sites.size(), String.size());
   CHECK_EQUAL(x.Basis1(), Sites[n]);
   CHECK_EQUAL(x.Basis2(), Sites[n]);
   int const Size = Sites.size();

   // The way we handle momentum here is probably incorrect, we should do the entire momentum
   // on site 0
   std::complex<double> r = std::exp(std::complex<double>(0.0, 1.0) * (Momentum / Size));

   // construct a list of added quantum numbers that are needed for the bond bases
   std::vector<BasisList> BondBasis(Size, BasisList(Sites[0].GetSymmetryList()));

   // Add the identity component to the bond basis
   for (int i = 0; i < Size; ++i)
   {
      BondBasis[i].push_back(QuantumNumber(Sites[0].GetSymmetryList()));
   }

   // Intermediate identity parts
   for (int i = 1; i <= n; ++i)
   {
      BondBasis[i].push_back(QuantumNumber(Sites[0].GetSymmetryList()));
   }

   // the actual operator
   for (int i = 0; i < Size; ++i)
   {
      BondBasis[i].push_back(x.TransformsAs());
   }

   // Assemble the operators
   GenericMPO Result(Size);
   for (int i = 0; i < Size; ++i)
   {
      CHECK_EQUAL(String[i].Basis1(), Sites[i]);
      CHECK_EQUAL(String[i].Basis2(), Sites[i]);

      Result[i] = OperatorComponent(Sites[i], Sites[i], BondBasis[i], BondBasis[smod(i+1,Size)]);
      Result[i](0,0) = i == 0 ? (r * String[i]) : r * String[i];

      if (i == 0 && i < n)
      {
	 Result[i](0,1) = SimpleOperator::make_identity(Sites[i]);
	 Result[i](1,2) = SimpleOperator::make_identity(Sites[i]);
      }
      else if (i == 0 && i == n)
      {
	 Result[i](0,1) = x;
	 Result[i](1,1) = SimpleOperator::make_identity(Sites[i]);
      }
      else if (i < n)
      {
	 Result[i](1,1) = SimpleOperator::make_identity(Sites[i]);
	 Result[i](2,2) = SimpleOperator::make_identity(Sites[i]);
      }
      else if (i == n)
      {
	 Result[i](1,1) = x;
	 Result[i](2,1) = SimpleOperator::make_identity(Sites[i]);
      }
      else
      {
	 Result[i](1,1) = SimpleOperator::make_identity(Sites[i]);
      }
   }
   return TriangularMPO(Result.data());
}

TriangularMPO OneCellStringOperator(std::vector<BasisList> const& Sites, 
				    std::vector<SimpleOperator> const& String,
				    std::vector<SimpleOperator> const& CellOp,
				    double Momentum)
{
   CHECK_EQUAL(Sites.size(), String.size());
   CHECK_EQUAL(CellOp.size(), String.size());
   int const Size = Sites.size();

   // The way we handle momentum here is probably incorrect, we should do the entire momentum
   // on site 0
   std::complex<double> r = std::exp(std::complex<double>(0.0, 1.0) * (Momentum / Size));

   // construct a list of added quantum numbers that are needed for the bond bases
   std::vector<BasisList> BondBasis(Size, BasisList(Sites[0].GetSymmetryList()));

   // Add the identity component to the bond basis
   for (int i = 0; i < Size; ++i)
   {
      BondBasis[i].push_back(QuantumNumber(Sites[0].GetSymmetryList()));
   }

   // Intermediate identity parts
   for (int i = 1; i < Size; ++i)
   {
      BondBasis[i].push_back(QuantumNumber(Sites[0].GetSymmetryList()));
   }

   // the actual operator
   for (int i = 0; i < Size; ++i)
   {
      BondBasis[i].push_back(QuantumNumber(Sites[0].GetSymmetryList()));
   }

   // Assemble the operators
   GenericMPO Result(Size);
   for (int i = 0; i < Size; ++i)
   {
      CHECK_EQUAL(String[i].Basis1(), Sites[i]);
      CHECK_EQUAL(String[i].Basis2(), Sites[i]);

      Result[i] = OperatorComponent(Sites[i], Sites[i], BondBasis[i], BondBasis[smod(i+1,Size)]);
      Result[i](0,0) = i == 0 ? (r * String[i]) : r * String[i];

      if (i == 0)
      {
	 Result[i](0,1) = CellOp[i];
	 Result[i](1,2) = SimpleOperator::make_identity(Sites[i]);
      }
      else if (i < Size-1)
      {
	 Result[i](1,1) = CellOp[i];
	 Result[i](2,2) = SimpleOperator::make_identity(Sites[i]);
      }
      else if (i == Size-1)
      {
	 Result[i](1,1) = CellOp[i];
	 Result[i](2,1) = SimpleOperator::make_identity(Sites[i]);
      }
   }
   return TriangularMPO(Result.data());
}

TriangularMPO OnePointOperator(std::vector<BasisList> const& Sites, 
                                    int n, SimpleOperator const& x, double Momentum)
{
   return OnePointStringOperator(Sites, MakeIdentityUnitCell(Sites), n, x, Momentum);
}

TriangularMPO TwoPointOperator(std::vector<BasisList> const& Sites, 
                                    int n1, SimpleOperator const& x1,
                                    int n2, SimpleOperator const& x2)
{
   DEBUG_TRACE(n1)(x1)(n2)(x2);
   // Normal order the sites.
   if (n1 > n2)
      return TwoPointOperator(Sites, n2, x2, n1, x1);
   else if (n1 == n2)
      return OnePointOperator(Sites, n1, x1*x2);

   int const Size = Sites.size();

   // construct a list of added quantum numbers that are needed for the bond bases
   std::vector<BasisList> BondBasis(Size, BasisList(Sites[0].GetSymmetryList()));

   // Add the identity component to the bond basis
   for (int i = 0; i < Size; ++i)
   {
      BondBasis[i].push_back(QuantumNumber(Sites[0].GetSymmetryList()));
   }

   // the actual operator
   for (int i = n1+1; i <= n2; ++i)
   {
      BondBasis[smod(i,Size)].push_back(x2.TransformsAs());
   }
   
   // Finally, the Hamiltonian component
   for (int i = 0; i < Size; ++i)
   {
      BondBasis[i].push_back(QuantumNumber(Sites[0].GetSymmetryList()));
   }

   // Assemble the operators
   GenericMPO Result(Size);
   for (int i = 0; i < Size; ++i)
   {
      Result[i] = OperatorComponent(Sites[i], Sites[i], BondBasis[i], BondBasis[smod(i+1,Size)]);
      Result[i](0,0) = SimpleOperator::make_identity(Sites[i]);
      Result[i](BondBasis[i].size()-1, BondBasis[smod(i+1,Size)].size()-1) 
         = SimpleOperator::make_identity(Sites[i]);
   }

   // now the operators.  Keep track of which component we insert them into
   std::vector<int> Loc(Size, 0);
   Result[smod(n1,Size)](0,1) = x1;
   ++Loc[smod(n1+1,Size)];
   for (int i = n1+1; i < n2; ++i)
   {
      Result[smod(i,Size)](Loc[smod(i,Size)], Loc[smod(i+1,Size)]+1) 
         = SimpleOperator::make_identity(Sites[smod(i,Size)]);
      ++Loc[smod(i+1,Size)];
   }
   //   TRACE(n1)(n2)(LinearAlgebra::Vector<int>(Loc));
   Result[smod(n2,Size)](Loc[smod(n2,Size)],Loc[smod(n2+1,Size)]+1) = x2;
   TriangularMPO TriResult(Result.data());
   //TRACE(TriResult);
   TriResult.debug_check_structure();
   return TriResult;
}

TriangularMPO TwoPointExponentialOperator(std::vector<BasisList> const& Sites, 
                                          int n1, SimpleOperator const& x1,
                                          int n2, SimpleOperator const& x2,
                                          std::complex<double> Factor)
{
   DEBUG_TRACE(n1)(x1)(n2)(x2);
   PRECONDITION(n1 < n2)(n1)(n2);
   PRECONDITION(n1 < int(Sites.size()));
   PRECONDITION(n2 < 2*int(Sites.size()));

   int const Size = Sites.size();

   // construct a list of added quantum numbers that are needed for the bond bases
   std::vector<BasisList> BondBasis(Size, BasisList(Sites[0].GetSymmetryList()));

   // Add the identity component to the bond basis
   for (int i = 0; i < Size; ++i)
   {
      BondBasis[i].push_back(QuantumNumber(Sites[0].GetSymmetryList()));
   }

   // The string component, transforms the same way as x1
   for (int i = 0; i < Size; ++i)
   {
      BondBasis[i].push_back(x1.TransformsAs());
   }
   
   // Finally, the Hamiltonian component
   for (int i = 0; i < Size; ++i)
   {
      BondBasis[i].push_back(QuantumNumber(Sites[0].GetSymmetryList()));
   }

   // Assemble the operators
   GenericMPO Result(Size);
   for (int i = 0; i < Size; ++i)
   {
      Result[i] = OperatorComponent(Sites[i], Sites[i], BondBasis[i], BondBasis[smod(i+1,Size)]);
      Result[i](0,0) = SimpleOperator::make_identity(Sites[i]);
      Result[i](1,1) = (i == n1 ? Factor : 1.0) * SimpleOperator::make_identity(Sites[i]);
      Result[i](2,2) = SimpleOperator::make_identity(Sites[i]);
   }

   Result[n1](0,1) = x1;
   Result[smod(n2,Size)](1,2) = x2;

   TriangularMPO TriResult(Result.data());
   TriResult.debug_check_structure();
   return TriResult;
}

TriangularMPO TwoPointStringOperator(std::vector<BasisList> const& Sites, 
                                     int n1, SimpleOperator const& x1,
                                     SimpleOperator const& String,
                                     int n2, SimpleOperator const& x2)
{
   return TwoPointStringOperator(Sites, n1, x1, std::vector<SimpleOperator>(Sites.size(), String),
                                 n2, x2);
}

TriangularMPO TwoPointStringOperator(std::vector<BasisList> const& Sites, 
                                     int n1, SimpleOperator const& x1,
                                     std::vector<SimpleOperator> const& String,
                                     int n2, SimpleOperator const& x2)
{
   CHECK(n1 < n2)("TwoPointStringOperator: error: sites must be normal ordered.")(n1)(n2);

   int const Size = Sites.size();

   // construct a list of added quantum numbers that are needed for the bond bases
   std::vector<BasisList> BondBasis(Size, BasisList(Sites[0].GetSymmetryList()));

   // Add the identity component to the bond basis
   for (int i = 0; i < Size; ++i)
   {
      BondBasis[i].push_back(QuantumNumber(Sites[0].GetSymmetryList()));
   }

   // the actual operator
   for (int i = n1+1; i <= n2; ++i)
   {
      BondBasis[smod(i,Size)].push_back(x2.TransformsAs());
   }
   
   // Finally, the Hamiltonian component
   for (int i = 0; i < Size; ++i)
   {
      BondBasis[i].push_back(QuantumNumber(Sites[0].GetSymmetryList()));
   }

   // Assemble the operators
   GenericMPO Result(Size);
   for (int i = 0; i < Size; ++i)
   {
      Result[i] = OperatorComponent(Sites[i], Sites[i], BondBasis[i], BondBasis[smod(i+1,Size)]);
      Result[i](0,0) = SimpleOperator::make_identity(Sites[i]);
      Result[i](BondBasis[i].size()-1, BondBasis[smod(i+1,Size)].size()-1) 
         = SimpleOperator::make_identity(Sites[i]);
   }

   // now the operators.  Keep track of which component we insert them into
   std::vector<int> Loc(Size, 0);
   Result[smod(n1,Size)](0,1) = x1;
   ++Loc[smod(n1+1,Size)];
   for (int i = n1+1; i < n2; ++i)
   {
      Result[smod(i,Size)](Loc[smod(i,Size)], Loc[smod(i+1,Size)]+1) = String[smod(i,Size)];
      ++Loc[smod(i+1,Size)];
   }
   Result[smod(n2,Size)](Loc[smod(n2,Size)],Loc[smod(n2+1,Size)]+1) = x2;
   return TriangularMPO(Result.data());
}
