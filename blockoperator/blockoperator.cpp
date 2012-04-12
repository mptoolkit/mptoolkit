// -*- C++ -*- $Id$

#include "blockoperator.h"

BlockOperator operator+(BlockOperator const& x, BlockOperator const& y)
{
   DEBUG_PRECONDITION(x.Basis1() == y.Basis1());
   DEBUG_PRECONDITION(x.Basis2() == y.Basis2());
   BlockOperator Result(x.Basis1(), x.Basis2(), x.sparse() + y.sparse(), x.dense() + y.dense());
   Result.canonicalize();
   return Result;
   }

BlockOperator operator-(BlockOperator const& x, BlockOperator const& y)
{
   DEBUG_PRECONDITION(x.Basis1() == y.Basis1());
   DEBUG_PRECONDITION(x.Basis2() == y.Basis2());
   BlockOperator Result(x.Basis1(), x.Basis2(), x.sparse() - y.sparse(), x.dense() - y.dense());
   Result.canonicalize();
   return Result;
}

std::ostream& operator<<(std::ostream& out, BlockOperator const& x)
{
   out << "BlockOperator:\n" << x.Basis1() << x.Basis2()
       << "sparse() =\n" << x.sparse()
       << "dense() =\n" << x.dense();
   return out;
}

void BlockOperator::canonicalize()
{
}

BlockOperator triple_prod(BlockOperator const& x, 
			  BlockOperator const& E,
			  HermitianProxy<BlockOperator const> y)
{
   BlockOperator Result(x.Basis1(), y.base().Basis1(),
			triple_prod(x.sparse(), E.sparse(), herm(y.base().sparse())), 
			triple_prod(x.dense(), E.sparse(), herm(y.base().sparse()))
			+ triple_prod(x.sparse(), E.sparse(), herm(y.base().dense()))
			+ triple_prod(x.dense(), E.sparse(), herm(y.base().dense()))
			+ triple_prod(x.sparse(), E.dense(), herm(y.base().sparse()))
			+ triple_prod(x.dense(), E.dense(), herm(y.base().sparse()))
			+ triple_prod(x.sparse(), E.dense(), herm(y.base().dense()))
			+ triple_prod(x.dense(), E.dense(), herm(y.base().dense())));
   Result.canonicalize();
   return Result;
}

BlockOperator triple_prod(BlockOperator const& x, 
			  BlockOperator const& E,
			  HermitianProxy<BlockOperator const> y,
			  QuantumNumber qxy,
			  QuantumNumber qEp)
{
   BlockOperator Result(x.Basis1(), y.base().Basis1(),
			triple_prod(x.sparse(), E.sparse(), herm(y.base().sparse()), qxy, qEp), 
			triple_prod(x.dense(), E.sparse(), herm(y.base().sparse()), qxy, qEp)
			+ triple_prod(x.sparse(), E.sparse(), herm(y.base().dense()), qxy, qEp)
			+ triple_prod(x.dense(), E.sparse(), herm(y.base().dense()), qxy, qEp)
			+ triple_prod(x.sparse(), E.dense(), herm(y.base().sparse()), qxy, qEp)
			+ triple_prod(x.dense(), E.dense(), herm(y.base().sparse()), qxy, qEp)
			+ triple_prod(x.sparse(), E.dense(), herm(y.base().dense()), qxy, qEp)
			+ triple_prod(x.dense(), E.dense(), herm(y.base().dense()), qxy, qEp));
   Result.canonicalize();
   return Result;
}

BlockOperator triple_prod(HermitianProxy<BlockOperator const> x, 
			  BlockOperator const& E,
			  BlockOperator const& y)
{
   BlockOperator Result(x.base().Basis2(), y.Basis2(),
			triple_prod(herm(x.base().sparse()), E.sparse(), y.sparse()), 
			triple_prod(herm(x.base().dense()), E.sparse(), y.sparse())
			+ triple_prod(herm(x.base().sparse()), E.sparse(), y.dense())
			+ triple_prod(herm(x.base().dense()), E.sparse(), y.dense())
			+ triple_prod(herm(x.base().sparse()), E.dense(), y.sparse())
			+ triple_prod(herm(x.base().dense()), E.dense(), y.sparse())
			+ triple_prod(herm(x.base().sparse()), E.dense(), y.dense())
			+ triple_prod(herm(x.base().dense()), E.dense(), y.dense()));
   Result.canonicalize();
   return Result;
}

BlockOperator triple_prod(HermitianProxy<BlockOperator const> x, 
			  BlockOperator const& E,
			  BlockOperator const& y,
			  QuantumNumber qxy,
			  QuantumNumber qEp)
{
   BlockOperator Result(x.base().Basis2(), y.Basis2(),
			triple_prod(herm(x.base().sparse()), E.sparse(), y.sparse(), qxy, qEp),
			triple_prod(herm(x.base().dense()), E.sparse(), y.sparse(), qxy, qEp)
			+ triple_prod(herm(x.base().sparse()), E.sparse(), y.dense(), qxy, qEp)
			+ triple_prod(herm(x.base().dense()), E.sparse(), y.dense(), qxy, qEp)
			+ triple_prod(herm(x.base().sparse()), E.dense(), y.sparse(), qxy, qEp)
			+ triple_prod(herm(x.base().dense()), E.dense(), y.sparse(), qxy, qEp)
			+ triple_prod(herm(x.base().sparse()), E.dense(), y.dense(), qxy, qEp)
			+ triple_prod(herm(x.base().dense()), E.dense(), y.dense(), qxy, qEp));
   Result.canonicalize();
   return Result;
}

