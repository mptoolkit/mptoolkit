

/*
  experimental code to implement TensorSum as a functor.

  TODO: the Nested argument is not necessary, and actually unworkable.
*/

template <typename S, typename T, typename B1, typename B2,
          typename Nested = DirectSum<typename interface<S>::value_type,
                                      typename interface<T>::value_type> >
struct TensorSum;

template <typename S, typename T>
typename TensorSum<S, T, B1, B2>::result_type
tensor_sum(S const& x, T const& y, B1 const B1, B2 const& b2)
{
   return TensorSum<S, T>()(x,y,b1,b2);
}

template <typename S, typename T, typename F>
typename TensorSum<S, T, B1, B2, FF>::result_type
tensor_sum(S const& x, T const& y,  B1 const B1, B2 const& b2, F const& f)
{
   return TensorSum<S, T, F>(f)(x,y,b1,b2);
}

template <typename S, typename T, typename B1,
          typename Nested = DirectSum<typename interface<S>::value_type,
                                      typename interface<T>::value_type> >
struct TensorColSum;

template <typename S, typename T, typename B1, typename F>
typename TensorColSum<S, T, B1>::result_type
tensor_col_sum(S const& x, T const& y, B1 const& b1)
{
   return TensorColSum<S, T>(x,y,b1);
}

template <typename S, typename T, typename B1, typename F>
typename TensorColSum<S, T, B1, F>::result_type
tensor_col_sum(S const& x, T const& y, B1 const& b1)
{
   return TensorColSum<S, T>(f)(x,y,b1);
}

template <typename S, typename T, typename B2,
          typename Nested = DirectSum<typename interface<S>::value_type,
                                      typename interface<T>::value_type> >
struct TensorColSum;

template <typename S, typename T, typename B2, typename F>
typename TensorRowSum<S, T, B2>::result_type
tensor_row_sum(S const& x, T const& y, B2 const& b2)
{
   return TensorRowSum<S, T>(x,y,b1);
}

template <typename S, typename T, typename B2, typename F>
typename TensorRowSum<S, T, B2, F>::result_type
tensor_row_sum(S const& x, T const& y, B2 const& b2)
{
   return TensorRowSum<S, T>(f)(x,y,b2);
}

template <typename T1, typename B1, typename B2, typename S1, 
          typename T2, typename S2, typename Nest>
struct TensorSum<IrredTensor<T1, B1, B2, S1>, IrredTensor<T2, B1, B2, S2> >
{
   typedef typename result_value<Nest>::type ResultS;
   typedef IrredTensor<typename ResultS::value_type, B1, B2, ResultS> result_type;
   typedef IrredTensor<T1, B1, B2, S1> const& first_argument_type;
   typedef IrredTensor<T2, B1, B2, S2> const& second_argument_type;
   typedef SumBasis<B1> const& third_argument_type;
   typedef SumBasis<B2> const& fourth_argument_type;

   TensorSum() {}
   TensorSum(Nest const& f) : f_(f) {}

   result_type operator()(first_argument_type x, second_argument_type y,
                          third_argument_type b1, fourth_argument_type b2);
};
