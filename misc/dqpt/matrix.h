// matrix.h
//
// A simple matrix class for PHYS4070/PHYS7270
// including IO

// include guard: prevent including this file multiple times
#if !defined(PHYS4070_MATRIX_H)
#define PHYS4070_MATRIX_H

#include <vector>      // we use a vector as the underlying storage
#include <cassert>     // for debug error checking
#include <iostream>    // for stream output
#include <iomanip>
#include <utility>     // for std::declval
#include <complex>     // for complex operator<< specialization
#include <initializer_list>

// number of digits of precision to show when printing a matrix
int const matrix_digits = 6;

template <typename T>
class matrix
{
   // public section: declare variables, functions etc
   // that are accessible to the outside world
   public:
      // default constructor
      matrix();

      // constructor for a matrix with a given size, with elements uninitialized
      matrix(unsigned rows, unsigned cols);

      // constructor for a matrix with a given size, with all elements initialized to x
      matrix(unsigned rows, unsigned cols, T x);

      // constructor that also takes an initializer list
      matrix(unsigned rows, unsigned cols, std::initializer_list<T> x);

      // constructor that takes a 2-dimensional initializer list.  This doesn't need
      // the number of rows and columns passed separately, since this is set by the 2D list,
      // eg a 2x3 matrix is matrix{{a,b,c}, {d,e,f}};
      matrix(std::initializer_list<std::initializer_list<T>> m);

      // copy constructor.  The compiler can automatically generate this
      matrix(matrix const&) = default;

      // copy-assignment operator.  The compiler can automatically generate this
      matrix& operator=(matrix const&) = default;

      // destructor.  The compiler can automatically generate this.
      ~matrix() = default;

      // move constructor.  We need to implement this ourselves to correctly
      // maintain the invariant that m_rows*m_cols == m_storage.size()
      matrix(matrix&& x);

      // move assignment operator.  We need to implement this ourselves to correctly
      // maintain the invariant that m_rows*m_cols == m_storage.size()
      matrix& operator=(matrix&& x);

      // member functions

      // resize the matrix to the given size, with elements uninitialized
      void resize(unsigned r, unsigned c);

      // resize the matrix to the given size, initializing the elements to x
      void resize(unsigned r, unsigned c, T x);

      // returns the number of rows in the matrix.
      // This is a 'const' function, as it doesn't need to
      // modify the internal state
      unsigned rows() const;

      // returns the number of columns in the matrix
      unsigned cols() const;

      // get an element from the matrix at position (r,c)
      T get(unsigned r, unsigned c) const;

      // set the matrix element at (r,c) to be x
      void set(unsigned r, unsigned c, T x);

      // set a submatrix starting from (r,c)
      template <typename U>
      void set_submatrix(unsigned r, unsigned c, matrix<U> const& m);

      // we can improve the notation for accessing the matrix.
      // The vector<T> class uses square brackets to access an
      // element, eg v[i].  The operator[] only allows 1 argument,
      // so it is not possible in C++ to use the notation m[i,j],
      // however we can instead implement the round-brackets
      // operator (known as the 'function call operator', because
      // it has the same notation as an ordinary function call).
      // This means we can access the matrix elements using
      // m(i,j) notation.  If we have a version that returns
      // the matrix element by reference, we can use it on the
      // left-hand side of an expression as well.

      // First overload returns by value and doesn't allow
      // modifications of the matrix, so it is a 'const' function
      T operator()(unsigned r, unsigned c) const;

      // Second overload is a non-const version that allows
      // us to modify the matrix element
      T& operator()(unsigned r, unsigned c);

      // reshape the matrix into a vector, with row-major format
      // this allows read-only access.  The non-const version is 'private',
      // since it allows a violation of the class invariants
      // (ie, resizing the underlying vector without updating the rows/columns)
      // so we don't want to call it directly, except under controlled circumstances
      std::vector<T> const& vector_view() const;

      // direct access to the matrix, useful for calling functions that operate
      // on raw memory, but 'dangerous' because there is no bounds checking
      T const* data() const;
      T* data();

   // private section: declare variables, functions etc
   // that are not accessible to the outside world
   private:
      std::vector<T>& vector_view();

      // variables to store the number of rows and columns of the matrix.
      // A common convention (by no means universal!) is to denote
      // member variables of a class with the prefix "m_".
      unsigned m_rows;
      unsigned m_cols;
      // the underlying storage of the matrix elements is in a vector.
      // Note, we refer to std::vector here rather than add a using
      // declaration, since we don't want to force anyone who uses the
      // matrix header to use of std:: namespace.
      std::vector<T> m_storage;
};

// Stream output.
// This writes the size of the matrix on the first line,
// followed by the matrix elements
template <typename T>
std::ostream& operator<<(std::ostream& out, matrix<T> const& m);

// Read in a matrix from the stream
template <typename T>
std::istream& operator>>(std::istream& in, matrix<T>& m);

// Arithmetic operations

// transpose

template <typename T>
matrix<T>
trans(matrix<T> const& x);

template <typename T>
matrix<T>
trans(matrix<T> const& x)
{
   matrix<T> Result(x.cols(), x.rows());
   for (unsigned i = 0; i < x.rows(); ++i)
   {
      for (unsigned j = 0; j < x.cols(); ++j)
      {
         Result(j,i) = x(i,j);
      }
   }
   return Result;
}

template <typename T>
matrix<std::complex<T>>
conj(matrix<std::complex<T>> x);

template <typename T>
matrix<std::complex<T>>
herm(matrix<std::complex<T>> x);

template <typename T>
matrix<std::complex<T>>
herm(matrix<std::complex<T>> x)
{
   return conj(trans(x));
}

// matrix addition

// Some 'magic' here: we use the C++-11 decltype() operator
// to automatically deduce the result type of T + U

template <typename T, typename U>
matrix<decltype(std::declval<T>()+std::declval<U>())>
operator+(matrix<T> const& x, matrix<U> const& y);

template <typename T, typename U>
matrix<decltype(std::declval<T>()-std::declval<U>())>
operator-(matrix<T> const& x, matrix<U> const& y);

template <typename T, typename U>
matrix<T>&
operator+=(matrix<T>& x, matrix<U> const& y);

template <typename T, typename U>
matrix<T>&
operator-=(matrix<T>& x, matrix<U> const& y);

// unary negation
template <typename T>
matrix<T>
operator-(matrix<T> x);

// matrix multiply by scalar

template <typename T, typename U>
matrix<decltype(std::declval<T>()*std::declval<U>())>
operator*(matrix<T> const& x, U const& y);

template <typename T, typename U>
matrix<decltype(std::declval<T>()*std::declval<U>())>
operator*(T const& x, matrix<U> const& y);

template <typename T, typename U>
matrix<T>&
operator*=(matrix<T>& x, U const& y);

// vector addition

template <typename T, typename U>
std::vector<decltype(std::declval<T>()+std::declval<U>())>
operator+(std::vector<T> const& x, std::vector<U> const& y);

template <typename T, typename U>
std::vector<decltype(std::declval<T>()-std::declval<U>())>
operator-(std::vector<T> const& x, std::vector<U> const& y);

template <typename T, typename U>
std::vector<T>&
operator+=(std::vector<T>& x, std::vector<U> const& y);

template <typename T, typename U>
std::vector<T>&
operator-=(std::vector<T>& x, std::vector<U> const& y);

// vector multiply by scalar

template <typename T, typename U>
std::vector<decltype(std::declval<T>()*std::declval<U>())>
operator*(std::vector<T> const& x, U const& y);

template <typename T, typename U>
std::vector<decltype(std::declval<T>()*std::declval<U>())>
operator*(T const& x, std::vector<U> const& y);

template <typename T, typename U>
std::vector<T>&
operator*=(std::vector<T>& x, U const& y);

// matrix-vector multiplication

template <typename T, typename U>
std::vector<decltype(std::declval<T>()*std::declval<U>())>
operator*(matrix<T> const& M, std::vector<U> const& x);

// matrix-matrix multiplication

template <typename T, typename U>
matrix<decltype(std::declval<T>()*std::declval<U>())>
operator*(matrix<T> const& x, matrix<U> const& y);

// vector output
template <typename T>
std::ostream& operator<<(std::ostream& out, std::vector<T> const& x);

// specialize for double and complex

std::ostream& operator<<(std::ostream& out, std::vector<double> const& x);
std::ostream& operator<<(std::ostream& out, std::vector<std::complex<double>> const& x);

// support for complex matrices

template <typename T>
matrix<std::complex<T>>
conj(matrix<std::complex<T>> x)
{
   for (unsigned i = 0; i < x.rows(); ++i)
   {
      for (unsigned j = 0; j < x.cols(); ++j)
      {
         x(i,j) = conj(x(i,j));
      }
   }
   return x;
}

// include template function definitions
#include "matrix.cc"

#endif
// end of include guard
