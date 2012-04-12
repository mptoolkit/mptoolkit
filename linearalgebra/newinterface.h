// -*- C++ -*- $Id$

#if !defined(NEWINTERFACE_H_CHUIRFHUIYT)
#define NEWINTERFACE_H_CHUIRFHUIYT

namespace LinearAlgebra
{

template <typename T, typename Enable = void>
struct interface
{
};

// unary functions

template <typename T>
struct Identity
{
   typedef T& result_type;
   typedef T& argument_type;
   T& operator()(T& x) const { return x; }
};

struct IdentityF
{
   template <typename T>
   struct result { typedef typename Identity<T>::result_type type; };

   template <typename T>
   typename Identity<T>::result_type operator()(T const& x) const
   {
      return Identity<T>()(x);
   }
};

template <typename T>
inline
typename Identity<T const>::result_type 
identity(T const& x)
{
   return Identity<T const>()(x);
}

template <typename T>
inline
typename Identity<T>::result_type 
identity(T& x)
{
   return Identity<T>()(x);
}

// is_identity_function

template <typename T, typename Enable = void>
struct is_identity_function : public boost::mpl::false_
{
};

template <typename T>
struct is_identity_function<Identity<T> > : public boost::mpl::true_
{
};

// Inverse

template <typename T, typename Enable = void>
struct Inverse { };

struct InverseF
{
   template <typename T>
   struct result { typedef typename Inverse<T>::result_type type; };

   template <typename T>
   typename Inverse<T>::result_type operator()(T const& x) const
   {
      return Inverse<T>()(x);
   }
};

template <typename T>
inline
typename T::result_type
inverse(T const& x)
{
   return Inverse<T>(x);
}

template <>
struct Inverse<double>
{
   typedef double result_type;
   typedef double argument_type;
   double operator()(double x) const { return 1.0 / x; }
};

// Negate

template <typename T, typename TInterface = interface<T>::type>
struct Negate
{
};

struct NegateF
{
   template <typename T>
   struct result { typedef typename Negate<T>::result_type type; };

   template <typename T>
   typename Negate<T>::result_type operator()(T const& x) const
   {
      return Negate<T>()(x);
   }
};

template <typename T>
inline
typename Negate<T>::result_type 
operator-(T const& x)
{
   return Negate<T>()(x);
}

// The inverse of Negate is itself

template <typename T>
struct Inverse<Negate<T> > : public Identity<Negate<T> >
{
};

// Trace

template <typename T, typename TInterface = interface<T>::type>
struct Trace
{
};

struct TraceF
{
   template <typename T>
   struct result { typedef typename Trace<T>::result_type type; };

   template <typename T>
   typename Trace<T>::result_type operator()(T const& x) const
   {
      return Trace<T>()(x);
   }
};

template <typename T>
inline
typename Trace<T>::result_type
trace(T const& x)
{
   return Trace<T>()(x);
}

// Norm2Sq

template <typename T, typename TInterface = interface<T>::type>
struct Norm2Sq
{
};

struct Norm2SqF
{
   template <typename T>
   struct result { typedef typename Norm2Sq<T>::result_type type; };

   template <typename T>
   typename Norm2Sq<T>::result_type operator()(T const& x) const
   {
      return Norm2Sq<T>()(x);
   }
};

template <typename T>
inline
typename Norm2Sq<T>::result_type
norm_2_sq(T const& x)
{
   return Norm2Sq<T>()(x);
}

// Norm1 - the 1-norm

template <typename T, typename TInterface = interface<T>::type>
struct Norm1
{
};

struct Norm1F
{
   template <typename T>
   struct result { typedef typename Norm1<T>::result_type type; };

   template <typename T>
   typename Norm1<T>::result_type operator()(T const& x) const
   {
      return Norm1<T>()(x);
   }
};

template <typename T>
inline
typename Norm1<T>::result_type
norm_1(T const& x)
{
   return Norm1<T>()(x);
}

// Norm2 - the 2-norm

// Norm2 has a default implementation, as sqrt(norm_2_sq(x)).
template <typename T, typename TInterface = interface<T>::type>
struct Norm2
{
   typedef typename Norm2Sq<T, TInterface>::result_type result_type;
   typedef T argument_type;
   result_type operator()(T const& x) const { return sqrt(norm_2_sq(x)); }
};

struct Norm2F
{
   template <typename T>
   struct result { typedef typename Norm2<T>::result_type type; };

   template <typename T>
   typename Norm2<T>::result_type operator()(T const& x) const
   {
      return Norm2<T>()(x);
   }
};

template <typename T>
inline
typename Norm2<T>::result_type
norm_2(T const& x)
{
   return Norm2<T>()(x);
}

// NormInf - the infinity-norm

template <typename T, typename TInterface = interface<T>::type>
struct NormInf
{
};

struct NormInfF
{
   template <typename T>
   struct result { typedef typename NormInf<T>::result_type type; };

   template <typename T>
   typename NormInf<T>::result_type operator()(T const& x) const
   {
      return NormInf<T>()(x);
   }
};

template <typename T>
inline
typename NormInf<T>::result_type
norm_inf(T const& x)
{
   return NormInf<T>()(x);
}

// Transpose

template <typename T, typename TInterface = interface<T>::type>
struct Transpose
{
};

struct TransposeF
{
   template <typename T>
   struct result { typedef typename Transpose<T>::result_type type; };

   template <typename T>
   typename Transpose<T>::result_type operator()(T const& x) const
   {
      return Transpose<T>()(x);
   }
};

template <typename T>
inline
typename Transpose<T>::result_type
transpose(T const& x)
{
   return Transpose<T>()(x);
}

// The inverse of Transpose is itself

template <typename T>
struct Inverse<Transpose<T> > : public Identity<Transpose<T> >
{
};

// Conj

template <typename T, typename TInterface = interface<T>::type>
struct Conj
{
};

struct ConjF
{
   template <typename T>
   struct result { typedef typename Conj<T>::result_type type; };

   template <typename T>
   typename Conj<T>::result_type operator()(T const& x) const
   {
      return Conj<T>()(x);
   }
};

template <typename T>
inline
typename Conj<T>::result_type
conj(T const& x)
{
   return Conj<T>()(x);
}

// The inverse of conj is itself

template <typename T>
struct Inverse<Conj<T> > : public Identity<Conj<T> >
{
};

// Real

template <typename T, typename TInterface = interface<T>::type>
struct Real
{
};

struct RealF
{
   template <typename T>
   struct result { typedef typename Real<T>::result_type type; };

   template <typename T>
   typename Real<T>::result_type operator()(T const& x) const
   {
      return Real<T>()(x);
   }
};

template <typename T>
inline
typename Real<T const>::result_type
real(T const& x)
{
   return Real<T const>()(x);
}

template <typename T>
inline
typename Real<T>::result_type
real(T& x)
{
   return Real<T>()(x);
}

// Imag

template <typename T, typename TInterface = interface<T>::type>
struct Imag
{
};

struct ImagF
{
   template <typename T>
   struct result { typedef typename Imag<T>::result_type type; };

   template <typename T>
   typename Imag<T>::result_type operator()(T const& x) const
   {
      return Imag<T>()(x);
   }
};

template <typename T>
inline
typename Imag<T const>::result_type
imag(T const& x)
{
   return Imag<T const>()(x);
}

template <typename T>
inline
typename Imag<T>::result_type
imag(T& x)
{
   return Imag<T>()(x);
}

// Herm

template <typename T, typename TInterface = interface<T>::type>
struct Herm
{
};

struct HermF
{
   template <typename T>
   struct result { typedef typename Herm<T>::result_type type; };

   template <typename T>
   typename Herm<T>::result_type operator()(T const& x) const
   {
      return Herm<T>()(x);
   }
};

template <typename T>
inline
typename Herm<T>::result_type
herm(T const& x)
{
   return Herm<T>()(x);
}

// The inverse of herm is itself

template <typename T>
struct Inverse<Herm<T> > : public Identity<Herm<T> >
{
};

// UnaryComposer

template <typename T1, typename T2, typename Enable = void>
struct UnaryComposer
{
   typedef typename T1::value_type        value_type;
   typedef typename T1::result_type       result_type;
   typedef typename T2::argument_type     argument_type;

   UnaryComposer(T1 const& f1_ = T1(), T2 const& f2_ = T2()) : first(f1_), second(f2_) {}

   result_type operator()(argument_type const& x) const { return first(second(x)); }

   T1 first;
   T2 second;
};

//
// binary functions
//

// unary_compose

template <typename T1, typename T2, typename Enable = void>
struct UnaryCompose 
{
   typedef T1 first_argument_type;
   typedef T2 second_argument_type;
   typedef UnaryComposer<T1, T2> result_type;
   result_type operator()(T1 const& x, T2 const& y) const { return result_type(x,y); }
};

template <typename T1, typename T2>
inline
typename UnaryCompose<T1, T2>::result_type
unary_compose(T1 const& x, T2 const& y)
{
   return UnaryCompose<T1, T2>()(x,y);
}

template <typename T>
struct UnaryCompose<Negate<T>, Negate<T> >
{
   typedef Negate<T> first_argument_type;
   typedef Negate<T> second_argument_type;
   typedef Identity<T> result_type;

   result_type operator()(Negate<T> const&, Negate<T> const&) const { return Identity<T>(); }
};

template <typename T1, typename T2>
struct Inverse<UnaryComposer<T1, T2>,
  typename boost::enable_if<
   boost::mpl::and_<exists<typename Inverse<T1>::result_type>,
		    exists<typename Inverse<T2>::result_type> > >::type>
   : public UnaryCompose<typename Inverse<T2>::result_type,
			 typename Inverse<T1>::result_type>
{
   typedef UnaryCompose<typename Inverse<T2>::result_type,
			typename Inverse<T1>::result_type> base;

   Inverse(UnaryComposer<T1, T2> const& x) : base(inverse(x.second), inverse(x.first)) {}

};

// Addition

template <typename S, typename T, typename SInterface = interface<S>::type, 
	  typename TInterface = interface<T>::type>
struct Addition
{
};

struct AdditionF
{
   template <typename S, typename T>
   struct result
   {
      typedef typename Addition<S, T>::result_type type;
   };

   template <typename S, typename T>
   typename Addition<S, T>::result_type operator()(S const& x, T const& y) const
   {
      return Addition<S, T>()(x,y);
   }
};

template <typename S, typename T>
inline
typename boost::disable_if<exists<typename Addition<S, T>::builtin>,
			   typename Addition<S, T>::result_type>
operator+(S const& x, T const& y)
{
   return Addition<S, T>()(x,y);
}

// Subtraction

template <typename S, typename T, typename SInterface = interface<S>::type, 
	  typename TInterface = interface<T>::type>
struct Subtraction
{
};

struct SubtractionF
{
   template <typename S, typename T>
   struct result
   {
      typedef typename Subtraction<S, T>::result_type type;
   };

   template <typename S, typename T>
   typename Subtraction<S, T>::result_type operator()(S const& x, T const& y) const
   {
      return Subtraction<S, T>()(x,y);
   }
};

template <typename S, typename T>
inline
typename boost::disable_if<exists<typename Subtraction<S, T>::builtin>,
			   typename Subtraction<S, T>::result_type>
operator-(S const& x, T const& y)
{
   return Subtraction<S, T>()(x,y);
}

// Multiplication

template <typename S, typename T, typename SInterface = interface<S>::type, 
	  typename TInterface = interface<T>::type>
struct Multiplication
{
};

struct MultiplicationF
{
   template <typename S, typename T>
   struct result
   {
      typedef typename Multiplication<S, T>::result_type type;
   };

   template <typename S, typename T>
   typename Multiplication<S, T>::result_type operator()(S const& x, T const& y) const
   {
      return Multiplication<S, T>()(x,y);
   }
};

template <typename S, typename T>
inline
typename boost::disable_if<exists<typename Multiplication<S, T>::builtin>,
			   typename Multiplication<S, T>::result_type>
operator*(S const& x, T const& y)
{
   return Multiplication<S, T>()(x,y);
}

//
// BindFirst and BindSecond
//
// Make a unary function out of a binary function, by binding a constant
// to one of the arguments.
//
// TODO: optimize for the case where the functor is empty
//

template <typename BinaryOp>
struct BinderFirst
{
   typedef BinaryOp binary_functor_type;
   typedef typename BinaryOp::first_argument_type   bound_type;
   typedef typename BinaryOp::second_argument_type  argument_type;
   typedef typename BinaryOp::value_type            value_type;
   typedef typename BinaryOp::result_type           result_type;

   BinderFirst(bound_type const& b_, binary_functor_type const& f_ = binary_functor_type()) 
     : b(b_), f(f_) {}

   result_type operator()(argument_type const& x) const { return f(b, x); }

   bound_type b;
   binary_functor_type f;
};

template <typename BinaryOp>
struct BinderSecond
{
   typedef BinaryOp binary_functor_type;
   typedef typename BinaryOp::first_argument_type argument_type;
   typedef typename BinaryOp::second_argument_type bound_type;
   typedef typename BinaryOp::value_type value_type;
   typedef typename BinaryOp::result_type result_type;

   BinderSecond(bound_type const& b_, binary_functor_type const& f_ = binary_functor_type()) 
     : b(b_), f(f_) {}

   result_type operator()(argument_type const& x) const { return f(x, b); }

   bound_type b;
   binary_functor_type f;
};

//
// bind_first, bind_second free functions
//

template <typename T>
struct BindFirst
{
   typedef T first_argument_type;
   typedef typename T::first_argument_type second_argument_type;
   typedef BinderFirst<T> result_type;
   result_type operator()(first_argument_type const& f, second_argument_type const& y) const
   { return result_type(y,f); }
};

template <typename T>
struct BindSecond
{
   typedef T first_argument_type;
   typedef typename T::second_argument_type second_argument_type;
   typedef BinderSecond<T> result_type;
   result_type operator()(first_argument_type const& f, second_argument_type const& y) const
   { return result_type(y,f); }
};

template <typename T, typename U>
inline
typename BindFirst<T>::result_type
bind_first(T const& f, U const& x)
{
   return BindFirst<T>()(f,x);
}

template <typename T, typename U>
inline
typename BindSecond<T>::result_type
bind_second(T const& f, U const& x)
{
   return BindSecond<T>()(f,x);
}

// InnerProd

// ParProd



//
// defaults for simple types
//

//
// Scalar
//
// top-level interface for types that act as a scalar
//

template <typename T>
struct Scalar
{
   typedef T type;
};



template <typename T>
struct interface<T, typename 
		 boost::enable_if<exists<typename PromoteTraits<T, T>::value_type> >::type>
{
   typedef Scalar<T> type;
};

template <typename S, typename T>
struct Addition<S, T, Scalar<S>, Scalar<T> >
{
   typedef typename PromoteTraits<S, T>::value_type result_type;

   typedef S first_argument_type;
   typedef T second_argument_type;

   result_type operator()(S x, T y) const { return x+y; }
};

template <typename S, typename T>
struct Subtraction<S, T, Scalar<S>, Scalar<T> >
{
   typedef typename PromoteTraits<S, T>::value_type result_type;

   typedef S first_argument_type;
   typedef T second_argument_type;

   result_type operator()(S x, T y) const { return x-y; }
};

template <typename S, typename T>
struct Multiplication<S, T, Scalar<S>, Scalar<T> >
{
   typedef typename PromoteTraits<S, T>::value_type result_type;

   typedef S first_argument_type;
   typedef T second_argument_type;

   result_type operator()(S x, T y) const { return x*y; }
};




} // namespace LinearAlgebra

#endif
