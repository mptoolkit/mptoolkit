// -*- C++ -*- $Id$

namespace statistics
{

template <class FwdIter>
double mean(FwdIter start, FwdIter end)
{
   size_t n = 0;
   double x = 0;
   while (start != end)
   {
      x += *start;
      ++n;
      ++start;
   }
   return n == 0 ? 0 : x / n;
}

template <class FwdIter>
double variance(FwdIter start, FwdIter end, double mean, int dof)
{
   size_t n = 0;
   double v = 0;
   while (start != end)
   {
      v += (mean - *start) * (mean - *start);
      ++n;
      ++start;
   }
   return v / (n - dof);
}

template <class FwdIter1, class FwdIter2, class FwdIter3>
linear_fit_result linear_fit(FwdIter1 x_start, FwdIter1 x_end, FwdIter2 y_start, FwdIter3 variance_start)
{
   linear_fit_result result;
   result.N = 0;
   double delta = 0;
   double X = 0, Y = 0, X2 = 0, Y2 = 0, XY = 0, E = 0;
   // X = sum_i {x_i / var_i}
   // Y = sum_i {y_i / var_i}
   // X2 = sum_i {x_i^2 / var_i}
   // Y2 = sum_i {y_i^2 / var_i}
   // XY = sum_i {x_i y_i / var_i}
   // E = sum_i {1 / var_i}
   FwdIter1 x = x_start;
   FwdIter2 y = y_start;
   FwdIter3 v = variance_start;
   while (x != x_end)
   {
      double e = 1.0 / (*v);  // e = 1.0/variance
      X += (*x) * e;
      Y += (*y) * e;
      X2 += (*x) * (*x) * e;
      Y2 += (*y) * (*y) * e;
      XY += (*x) * (*y) * e;
      E += e;
      ++result.N;
      ++x;
      ++y;
      ++v;
   }

   delta = E * X2 - X * X;
   result.m = (1 / delta) * (E * XY - X * Y);
   result.c = (1 / delta) * (X2 * Y - X * XY);
   result.variance_m = (1 / delta) * E;
   result.variance_c = (1 / delta) * X2;

   // calculate chi2 and the variance
   x = x_start;
   y = y_start;
   v = variance_start;
   result.chi2 = 0;
   result.variance = 0;
   while (x != x_end)
   {
      double t = (*y - (result.m * (*x) + result.c));
      result.chi2 += t * t / (*v);
      result.variance += t * t;
      ++x;
      ++y;
      ++v;
   }
   result.variance /= (result.N - 2);  // 2 dof in the calculation, m & c

   return result;
}

template <class FwdIter1, class FwdIter2>
double chi2(FwdIter1 x_start, FwdIter1 x_end, FwdIter2 y_start)
{
   double result = 0;
   while (x_start != x_end)
   {
      double t = (*x_start - *y_start);
      result += t * t;
      ++x_start;
      ++y_start;
   }
   return result;
}

} // namespace

