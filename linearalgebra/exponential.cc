// -*- C++ -*- $Id$

namespace LinearAlgebra
{

namespace Private
{

void Exponentiate(double t, int Size, std::complex<double> const* H, int ldH, std::complex<double>* R, int ldR);

} // namespace Private

template <typename M, typename Mi>
struct ImplementExponentiate<M, Concepts::ContiguousMatrix<std::complex<double>, RowMajor, Mi>>
{
   typedef Matrix<std::complex<double> > result_type;
   result_type operator()(double t, M const& m) const
   {
      CHECK_EQUAL(size1(m), size2(m));

      double Norm = norm_frob(m);
      if (Norm < 1E-10)
      {
         // 1st order Taylor series approximation to the exponential is exp(m) = 1+m
         result_type r = m;
         for (unsigned i = 0; i < size1(r); ++i)
         {
            r(i,i) += 1.0;
         }
         return r;
      }

      result_type R(size1(m), size2(m));
      Private::Exponentiate(t, size1(m), data(m), size1(m), data(R), size1(R));
      return R;
   }
};

} // namespace LinearAlgebra
