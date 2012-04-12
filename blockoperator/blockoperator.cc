// -*- C++ -*- $Id$

inline
BlockOperator::BlockOperator()
{
}

inline
BlockOperator::BlockOperator(VectorBasis const& b1, VectorBasis const& b2, QuantumNumber const& q)
   : Basis1_(b1), Basis2_(b2), 
     SparsePart_(b1.Basis(), b2.Basis(), q), 
     DensePart_(b1.Basis(), b2.Basis(), q)
{
}

inline
BlockOperator::BlockOperator(VectorBasis const& b1, VectorBasis const& b2, sparse_type const& s)
   : Basis1_(b1), Basis2_(b2), 
     SparsePart_(s), 
     DensePart_(b1.Basis(), b2.Basis(), s.TransformsAs())
{
   DEBUG_PRECONDITION_EQUAL(s.Basis1(), b1.Basis());
   DEBUG_PRECONDITION_EQUAL(s.Basis2(), b2.Basis());
}

inline
BlockOperator::BlockOperator(VectorBasis const& b1, VectorBasis const& b2, dense_type const& d)
   : Basis1_(b1), Basis2_(b2), 
     SparsePart_(b1.Basis(), b2.Basis(), d.TransformsAs()), 
     DensePart_(d)
{
   DEBUG_PRECONDITION_EQUAL(d.Basis1(), b1.Basis());
   DEBUG_PRECONDITION_EQUAL(d.Basis2(), b2.Basis());
}

inline
BlockOperator::BlockOperator(VectorBasis const& b1, VectorBasis const& b2, 
			     sparse_type const& s, dense_type const& d)
   : Basis1_(b1), Basis2_(b2), SparsePart_(s), DensePart_(d)
{
   DEBUG_PRECONDITION_EQUAL(s.Basis1(), b1.Basis());
   DEBUG_PRECONDITION_EQUAL(s.Basis2(), b2.Basis());
   DEBUG_PRECONDITION_EQUAL(d.Basis1(), b1.Basis());
   DEBUG_PRECONDITION_EQUAL(d.Basis2(), b2.Basis());
   DEBUG_PRECONDITION_EQUAL(s.TransformsAs(), d.TransformsAs());
}   

inline
BlockOperator& BlockOperator::operator+=(BlockOperator const& x)
{
   SparsePart_ += x.SparsePart_;
   DensePart_ += x.DensePart_;
   this->canonicalize();
   return *this;
}

inline
BlockOperator& BlockOperator::operator-=(BlockOperator const& x)
{
   SparsePart_ -= x.SparsePart_;
   DensePart_ -= x.DensePart_;
   this->canonicalize();
   return *this;
}

inline
BlockOperator& BlockOperator::operator*=(double a)
{
   SparsePart_ *= a;
   DensePart_ *= a;
   return *this;
}

#if 0
inline
bool operator==(BlockOperator const& x, BlockOperator const& y)
{
   return x.sparse() == y.sparse() && x.dense() == y.dense();
}

inline
bool operator!=(BlockOperator const& x, BlockOperator const& y)
{
   return x.sparse() != y.sparse() || x.dense() != y.dense();
}
#endif

inline
bool equal(BlockOperator const& x, BlockOperator const& y, double tol)
{
   return equal(x.sparse(), y.sparse(), tol) && equal(x.dense(), y.dense(), tol);
}

inline
BlockOperator operator-(BlockOperator const& x)
{
   return BlockOperator(x.Basis1(), x.Basis2(), -x.sparse(), -x.dense());
}

inline
BlockOperator::value_type trace(BlockOperator const& x)
{
   return trace(x.sparse()) + trace(x.dense());
}

inline
double norm_frob_sq(BlockOperator const& x)
{
   return norm_frob_sq(x.sparse()) + norm_frob_sq(x.dense());
}

inline
double norm_frob(BlockOperator const& x)
{
   return std::sqrt(norm_frob_sq(x));
}

inline
BlockOperator operator*(double a, BlockOperator const& x)
{
   return BlockOperator(x.Basis1(), x.Basis2(), a * x.sparse(), a * x.dense());
}

inline
BlockOperator operator*(BlockOperator const& x, double a)
{
   return BlockOperator(x.Basis1(), x.Basis2(), x.sparse() * a, x.dense() * a);
}

inline
BlockOperator::value_type inner_prod(BlockOperator const& x, BlockOperator const& y)
{
   return inner_prod(x.sparse(), y.sparse())
      + inner_prod(x.sparse(), y.dense())
      + inner_prod(x.dense(), y.sparse())
      + inner_prod(x.dense(), y.dense());
}

inline
BlockOperator prod(BlockOperator const& x, BlockOperator const& y, 
		   QuantumNumber const& q)
{
   BlockOperator Result(x.Basis1(), y.Basis2(),
			prod(x.sparse(), y.sparse(), q),
			prod(x.dense(), y.dense(), q)
			+ prod(x.sparse(), y.dense(), q)
			+ prod(x.dense(), y.sparse(), q));
   return Result;
}

inline
HermitianProxy<BlockOperator const> herm(BlockOperator const& x)
{
   return HermitianProxy<BlockOperator const>(x);
}

#if 0
inline
BlockOperator scalar_prod(BlockOperator const& x, HermitianProxy<BlockOperator const> y)
{
   BlockOperator Result(x.Basis1(), y.base().Basis1(),
			scalar_prod(x.sparse(), herm(y.base().sparse())),
			scalar_prod(x.dense(), herm(y.base().dense()))
			+ scalar_prod(x.sparse(), herm(y.base().dense()))
			+ scalar_prod(x.dense(), herm(y.base().sparse())));
   Result.canonicalize();
   return Result;
}

inline
BlockOperator scalar_prod(HermitianProxy<BlockOperator const> x, BlockOperator const& y)
{
    BlockOperator Result(x.base().Basis2(), y.Basis2(),
			scalar_prod(herm(x.base().sparse()), y.sparse()),
			scalar_prod(herm(x.base().dense()), y.dense())
			+ scalar_prod(herm(x.base().sparse()), y.dense())
			 + scalar_prod(herm(x.base().dense()), y.sparse()));
   Result.canonicalize();
   return Result;
}
#endif

inline
BlockOperator adjoint(BlockOperator const& x)
{
   return BlockOperator(x.Basis2(), x.Basis1(), adjoint(x.sparse()), adjoint(x.dense()));
}

inline
BlockOperator inv_adjoint(BlockOperator const& x)
{
   return BlockOperator(x.Basis2(), x.Basis1(), inv_adjoint(x.sparse()), inv_adjoint(x.dense()));
}

inline
BlockOperator tensor_prod(BlockOperator const& x, BlockOperator const& y,
			  VectorProductBasis const& b1,
			  VectorProductBasis const& b2,
			  QuantumNumber const& q)
{
   BlockOperator Result(b1.Basis(), b2.Basis(),
			tensor_prod(x.sparse(), y.sparse(), b1.ProdBasis(), b2.ProdBasis(), q),
			tensor_prod(x.dense(), y.sparse(), b1.ProdBasis(), b2.ProdBasis(), q)
			+ tensor_prod(x.dense(), y.dense(), b1.ProdBasis(), b2.ProdBasis(), q)
			+ tensor_prod(x.sparse(), y.dense(), b1.ProdBasis(), b2.ProdBasis(), q));
   Result.canonicalize();
   return Result;
}

inline
BlockOperator tensor_sum(BlockOperator const& x, BlockOperator const& y,
			 VectorSumBasis const& b1, VectorSumBasis const& b2)
{
   BlockOperator Result(b1.Basis(), b2.Basis(),
			tensor_sum(x.sparse(), y.sparse(), b1.SumBasis(), b2.SumBasis()),
			tensor_sum(x.dense(), y.dense(), b1.SumBasis(), b2.SumBasis()));
   Result.canonicalize();
   return Result;
}

inline
BlockOperator tensor_row_sum(BlockOperator const& x, BlockOperator const& y,
			     VectorSumBasis const& b1)
{
   BlockOperator Result(b1.Basis(), x.Basis2(),
			tensor_row_sum(x.sparse(), y.sparse(), b1.SumBasis()),
			tensor_row_sum(x.dense(), y.dense(), b1.SumBasis()));
   Result.canonicalize();
   return Result;
}

inline
BlockOperator tensor_col_sum(BlockOperator const& x, BlockOperator const& y,
			     VectorSumBasis const& b2)
{
   BlockOperator Result(x.Basis1(), b2.Basis(),
			tensor_col_sum(x.sparse(), y.sparse(), b2.SumBasis()),
			tensor_col_sum(x.dense(), y.dense(), b2.SumBasis()));
   Result.canonicalize();
   return Result;
}

