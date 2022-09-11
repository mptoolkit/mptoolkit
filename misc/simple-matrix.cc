//
// Definitions of the matrix member functions.
// because these are template functions, the definition
// MUST be included in the source file.  We could also put these
// in a different file and #include it here.
// (If they were NOT template functions, we could put the definitions
// in a separate .cpp file and compile them separately.)

// default constuctor.
// Here we initialize the member variables.
// By default, we get a zero by zero matrix.
template <typename T>
matrix<T>::matrix()
   : m_rows(0),
     m_cols(0),
     m_storage()   // this is redundant, it would be called automatically
{
}

// Constructor for a rows x cols matrix.
template <typename T>
matrix<T>::matrix(unsigned rows, unsigned cols)
   : m_rows(rows),
     m_cols(cols),
     m_storage(rows*cols)
{
}

// Constructor for a rows x cols matrix initialized to x
template <typename T>
inline
matrix<T>::matrix(unsigned rows, unsigned cols, T x)
   : m_rows(rows),
     m_cols(cols),
     m_storage(rows*cols, x)
{
}

template <typename T>
inline
matrix<T>::matrix(unsigned rows, unsigned cols, std::initializer_list<T> x)
   : m_rows(rows),
     m_cols(cols),
     m_storage(x)
{
   assert(rows*cols == x.size());
}

template <typename T>
inline
matrix<T>::matrix(std::initializer_list<std::initializer_list<T>> m)
   : m_rows(m.size()),
     m_cols(m.size() == 0 ? 0 : m.begin()->size())
{
   m_storage.reserve(m_rows*m_cols);  // reserve enough memory for the matrix
   for (auto x : m)
   {
      assert(x.size() == m_cols);
      m_storage.insert(m_storage.end(), x);
   }
}

template <typename T>
inline
matrix<T>::matrix(matrix<T>&& x)
   : m_rows(std::move(x.m_rows)),
     m_cols(std::move(x.m_cols)),
     m_storage(std::move(x.m_storage))
{
   x.m_rows = 0;
   x.m_cols = 0;
}

template <typename T>
inline
matrix<T>&
matrix<T>::operator=(matrix<T>&& x)
{
   m_rows = std::move(x.m_rows);
   m_cols = std::move(x.m_cols);
   m_storage = std::move(x.m_storage);
   x.m_rows = 0;
   x.m_cols = 0;
   return *this;
}

// member function

template <typename T>
void
matrix<T>::resize(unsigned r, unsigned c)
{
   m_rows = r;
   m_cols = c;
   m_storage = std::vector<T>(m_rows*m_cols);
}

template <typename T>
void
matrix<T>::resize(unsigned r, unsigned c, T x)
{
   m_rows = r;
   m_cols = c;
   m_storage = std::vector<T>(m_rows*m_cols, x);
}

// member function to get an element from the matrix
template <typename T>
T matrix<T>::get(unsigned r, unsigned c) const
{
   // check that the indices are in range
   assert(r >= 0 && r < m_rows);
   assert(c >= 0 && c < m_cols);
   return m_storage[r*m_cols + c];
}

// member function to set an element of the matrix to x
template <typename T>
void matrix<T>::set(unsigned r, unsigned c, T x)
{
   // check that the indices are in range
   assert(r >= 0 && r < m_rows);
   assert(c >= 0 && c < m_cols);
   m_storage[r*m_cols + c] = x;
}

// member function to assign a submatrix
template <typename T>
template <typename U>
void
matrix<T>::set_submatrix(unsigned r, unsigned c, matrix<U> const& m)
{
   // make sure the range fits within the matrix
   assert(r+m.rows() <= m_rows);
   assert(c+m.cols() <= m_cols);

   for (unsigned i = 0; i < m.rows(); ++i)
   {
      for (unsigned j = 0; j < m.cols(); ++j)
      {
         this->set(i+r, j+c, m(i,j));
      }
   }
}

template <typename T>
unsigned matrix<T>::rows() const
{
   return m_rows;
}

template <typename T>
unsigned matrix<T>::cols() const
{
   return m_cols;
}

// the function call operator
template <typename T>
inline
T matrix<T>::operator()(unsigned r, unsigned c) const
{
   // check that the indices are in range
   assert(r >= 0 && r < m_rows);
   assert(c >= 0 && c < m_cols);
   return m_storage[r*m_cols + c];
}

template <typename T>
inline
T& matrix<T>::operator()(unsigned r, unsigned c)
{
   // check that the indices are in range
   assert(r >= 0 && r < m_rows);
   assert(c >= 0 && c < m_cols);
   return m_storage[r*m_cols + c];
}

template <typename T>
inline
T const* matrix<T>::data() const
{
   return m_storage.data();
}

template <typename T>
inline
T* matrix<T>::data()
{
   return m_storage.data();
}

template <typename T>
inline
std::vector<T> const&
matrix<T>::vector_view() const
{
   return m_storage;
}

template <typename T>
inline
std::vector<T>&
matrix<T>::vector_view()
{
   return m_storage;
}


// definitions of free functions

// Stream output
// utility function to print a complex number
std::string
format_complex(std::complex<double> x, int digits = matrix_digits)
{
  std::ostringstream s;
  s.precision(std::streamsize(digits));
  s << x.real() << std::showpos << x.imag() << "i";
  return s.str();
}

template <typename T>
std::ostream& operator<<(std::ostream& out, matrix<T> const& m)
{
   out << m.rows() << " " << m.cols() << "\n";
   for (unsigned i = 0; i < m.rows(); ++i)
   {
      for (unsigned j = 0; j < m.cols(); ++j)
      {
         out << std::setw(4) << m(i,j) << " ";
      }
      out << std::endl;  // new line after every row
   }
   return out;
}

// specialization for matrix<double>
inline
std::ostream& operator<<(std::ostream& out, matrix<double> const& m)
{
   // save the stream format flags
   std::ios::fmtflags oldflags = out.flags();
   std::streamsize oldprecision = out.precision();

   // set the flags to what we want for matrix output
   out.precision(std::streamsize(matrix_digits));
   out.unsetf(std::ios::floatfield);   // default
   out.setf(std::ios::left);

   out << "size " << m.rows() << " " << m.cols() << "\n";
   for (unsigned i = 0; i < m.rows(); ++i)
   {
      out << "[ ";
      for (unsigned j = 0; j < m.cols(); ++j)
      {
         out << std::setw(std::streamsize(matrix_digits+6)) << m(i,j) << " ";
      }
      out << "]\n"; // new line after every row
   }

   // reset stream flags
   out.precision(oldprecision);
   out.flags(oldflags);
   return out;
}

// specialization for matrix<complex<double>>
inline
std::ostream& operator<<(std::ostream& out, matrix<std::complex<double>> const& m)
{
   // save the stream format flags
   std::ios::fmtflags oldflags = out.flags();
   std::streamsize oldprecision = out.precision();

   // set the flags to what we want for matrix output
   out.precision(std::streamsize(matrix_digits));
   out.unsetf(std::ios::floatfield);   // default
   out.setf(std::ios::left);

   out << "size " << m.rows() << " " << m.cols() << "\n";
   for (unsigned i = 0; i < m.rows(); ++i)
   {
      out << "[ ";
      for (unsigned j = 0; j < m.cols(); ++j)
      {
 	 out << std::setw(std::streamsize(matrix_digits*2+12)) << format_complex(m(i,j)) << " ";
      }
      out << "]\n"; // new line after every row
   }

   // reset stream flags
   out.precision(oldprecision);
   out.flags(oldflags);
   return out;
}


template <typename T>
std::istream& operator>>(std::istream& in, matrix<T>& m)
{
   unsigned Rows, Cols;
   in >> Rows >> Cols;
   m = matrix<T>(Rows, Cols);
   for (unsigned i = 0; i < m.rows(); ++i)
   {
      for (unsigned j = 0; j < m.cols(); ++j)
      {
         in >> m(i,j);
      }
   }
   return in;
}

// matrix arithmetic operations

template <typename T, typename U>
matrix<decltype(std::declval<T>()+std::declval<U>())>
operator+(matrix<T> const& x, matrix<U> const& y)
{
   typedef decltype(std::declval<T>()+std::declval<U>()) ResultType;

   assert(x.rows() == y.rows());
   assert(x.cols() == y.cols());

   matrix<ResultType> Result(x.rows(), x.cols());
   for (unsigned i = 0; i < Result.rows(); ++i)
   {
      for (unsigned j = 0; j < Result.cols(); ++j)
      {
	 Result(i,j) = x(i,j) + y(i,j);
      }
   }
   return Result;
}

template <typename T, typename U>
matrix<decltype(std::declval<T>()-std::declval<U>())>
operator-(matrix<T> const& x, matrix<U> const& y)
{
   typedef decltype(std::declval<T>()-std::declval<U>()) ResultType;

   assert(x.rows() == y.rows());
   assert(x.cols() == y.cols());

   matrix<ResultType> Result(x.rows(), x.cols());
   for (unsigned i = 0; i < Result.rows(); ++i)
   {
      for (unsigned j = 0; j < Result.cols(); ++j)
      {
	 Result(i,j) = x(i,j) - y(i,j);
      }
   }
   return Result;
}

template <typename T>
matrix<T>
operator-(matrix<T> x)
{
   for (unsigned i = 0; i < x.rows(); ++i)
   {
      for (unsigned j = 0; j < x.cols(); ++j)
      {
	 x(i,j) = -x(i,j);
      }
   }
   return x;
}


template <typename T, typename U>
matrix<T>&
operator+=(matrix<T>& x, matrix<U> const& y)
{
   assert(x.rows() == y.rows());
   assert(x.cols() == y.cols());

   for (unsigned i = 0; i < x.rows(); ++i)
   {
      for (unsigned j = 0; j < x.cols(); ++j)
      {
	 // Note brace {} initialization avoids narrowing conversions
	 x(i,j) += T{y(i,j)};
      }
   }
   return x;
}

template <typename T, typename U>
matrix<T>&
operator-=(matrix<T>& x, matrix<U> const& y)
{
   assert(x.rows() == y.rows());
   assert(x.cols() == y.cols());

   for (unsigned i = 0; i < x.rows(); ++i)
   {
      for (unsigned j = 0; j < x .cols(); ++j)
      {
	 x(i,j) -= T{y(i,j)};
      }
   }
   return x;
}

template <typename T, typename U>
matrix<decltype(std::declval<T>()*std::declval<U>())>
operator*(matrix<T> const& x, U const& y)
{
   typedef decltype(std::declval<T>()*std::declval<U>()) ResultType;

   matrix<ResultType> Result(y.rows(), y.cols());
   for (unsigned i = 0; i < Result.rows(); ++i)
   {
      for (unsigned j = 0; j < Result.cols(); ++j)
      {
	 Result(i,j) = x(i,j) * y;
      }
   }
   return Result;
}

template <typename T, typename U>
matrix<decltype(std::declval<T>()*std::declval<U>())>
operator*(T const& x, matrix<U> const& y)
{
   typedef decltype(std::declval<T>()*std::declval<U>()) ResultType;

   matrix<ResultType> Result(y.rows(), y.cols());
   for (unsigned i = 0; i < Result.rows(); ++i)
   {
      for (unsigned j = 0; j < Result.cols(); ++j)
      {
	 Result(i,j) = x * y(i,j);
      }
   }
   return Result;
}

template <typename T, typename U>
matrix<T>&
operator*=(matrix<T>& x, U const& y)
{
   for (unsigned i = 0; i < x.rows(); ++i)
   {
      for (unsigned j = 0; j < x.cols(); ++j)
      {
	 x(i,j) *= y;
      }
   }
   return x;
}

// vector arithmetic operations

template <typename T, typename U>
std::vector<decltype(std::declval<T>()+std::declval<U>())>
operator+(std::vector<T> const& x, std::vector<U> const& y)
{
   typedef decltype(std::declval<T>()+std::declval<U>()) ResultType;

   assert(x.size() == y.size());

   std::vector<ResultType> Result(x.size());
   for (unsigned i = 0; i < x.size(); ++i)
   {
      Result[i] = x[i] + y[i];
   }
   return Result;
}

template <typename T, typename U>
std::vector<decltype(std::declval<T>()-std::declval<U>())>
operator-(std::vector<T> const& x, std::vector<U> const& y)
{
   typedef decltype(std::declval<T>()-std::declval<U>()) ResultType;

   assert(x.size() == y.size());

   std::vector<ResultType> Result(x.size());
   for (unsigned i = 0; i < x.size(); ++i)
   {
      Result[i] = x[i] - y[i];
   }
   return Result;
}

template <typename T, typename U>
std::vector<decltype(std::declval<T>()*std::declval<U>())>
operator*(std::vector<T> const& x, U const& y)
{
   typedef decltype(std::declval<T>()*std::declval<U>()) ResultType;

   std::vector<ResultType> Result(x.size());
   for (unsigned i = 0; i < x.size(); ++i)
   {
      Result[i] = x[i] * y;
   }
   return Result;
}

template <typename T, typename U>
std::vector<decltype(std::declval<T>()*std::declval<U>())>
operator*(T const& x, std::vector<U> const& y)
{
   typedef decltype(std::declval<T>()*std::declval<U>()) ResultType;

   std::vector<ResultType> Result(x.size());
   for (unsigned i = 0; i < x.size(); ++i)
   {
      Result[i] = x * y[i];
   }
   return Result;
}

template <typename T, typename U>
std::vector<T>&
operator*=(std::vector<T>& x, U const& y)
{
   for (unsigned i = 0; i < x.size(); ++i)
   {
      x[i] *= y;
   }
   return x;
}

// matrix-vector multiplication

template <typename T, typename U>
std::vector<decltype(std::declval<T>()*std::declval<U>())>
operator*(matrix<T> const& M, std::vector<U> const& x)
{
   typedef decltype(std::declval<T>()*std::declval<U>()) ResultType;

   assert(M.cols() == x.size());

   std::vector<ResultType> Result(M.rows());
   for (unsigned i = 0; i < Result.size(); ++i)
   {
      ResultType Temp = ResultType{};

      for (unsigned j = 0; j < M.cols(); ++j)
      {
	 Temp += M(i,j) * x[j];
      }
      Result[i] = Temp;
   }
   return Result;
}

// matrix-matrix multiplication

template <typename T, typename U>
matrix<decltype(std::declval<T>()*std::declval<U>())>
operator*(matrix<T> const& x, matrix<U> const& y)
{
   typedef decltype(std::declval<T>()*std::declval<U>()) ResultType;

   assert(x.cols() == y.rows());

   matrix<ResultType> Result(x.rows(), y.cols());
   for (unsigned i = 0; i < Result.rows(); ++i)
   {
      for (unsigned j = 0; j < Result.cols(); ++j)
      {
	 ResultType Temp = ResultType{};

	 for (unsigned k = 0; k < x.cols(); ++k)
	 {
	    Temp += x(i,k) * y(k,j);
	 }

	 Result(i,j) = Temp;
      }
   }
   return Result;
}

// vector output

inline
std::ostream& operator<<(std::ostream& out, std::vector<double> const& x)
{
   // save the stream format flags
   std::ios::fmtflags oldflags = out.flags();
   int oldprecision = out.precision();

   // set the flags to what we want for matrix output
   out.precision(std::streamsize(matrix_digits));
   out.unsetf(std::ios::floatfield);   // default
   out.setf(std::ios::left);

   out << "size " << x.size() << "\n[ ";
   for (unsigned i = 0; i < x.size(); ++i)
   {
      out << std::setw(std::streamsize(matrix_digits+6)) << x[i] << " ";
   }
   out << "]\n";

   // reset stream flags
   out.precision(oldprecision);
   out.flags(oldflags);
   return out;
}

// vector output

inline
std::ostream& operator<<(std::ostream& out, std::vector<std::complex<double>> const& x)
{
   // save the stream format flags
   std::ios::fmtflags oldflags = out.flags();
   int oldprecision = out.precision();

   // set the flags to what we want for matrix output
   out.precision(std::streamsize(matrix_digits));
   out.unsetf(std::ios::floatfield);   // default
   out.setf(std::ios::left);

   out << "size " << x.size() << "\n[ ";
   for (unsigned i = 0; i < x.size(); ++i)
   {
     out << std::setw(std::streamsize(matrix_digits*2+12)) << format_complex(x[i]) << " ";
   }
   out << "]\n";

   // reset stream flags
   out.precision(oldprecision);
   out.flags(oldflags);
   return out;
}
