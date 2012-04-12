// -*- C++ -*- $Id$

//
// iter_assign
//

// row-major

template <typename I1, typename I2>
void iter_assign(I1 i1, I2 i2, matrix_iter_row_dense, matrix_iter_row_major<sparse>)
{
   iter_zero(i1);
   while (i2)
   {
      add(i1[i2.index1()], *i2);
      ++i1;
      ++i2;
   }
}

template <typename I1, typename I2>
void iter_assign(I1 i1, I2 i2, matrix_iter_row_dense, matrix_iter_row_major<injective>)
{
   iter_zero(i1);
   while (i2)
   {
      assign(i1[i2.index1()], *i2);
      ++i1;
      ++i2;
   }
}

template <typename I1, typename I2>
void iter_assign(I1 i1, I2 i2, matrix_iter_row_dense, matrix_iter_row_major<dense>)
{
   while (i2)
   {
      assign(*i1, *i2);
      ++i1;
      ++i2;
   }
}

// row major, coordinate right-hand-side

template <typename I1, typename I2>
void iter_assign(I1 i1, I2 i2, matrix_iter_row_dense, matrix_iter_coord<sparse>)
{
   iter_zero(i1);
   while (i2)
   {
      add(i1[i2.index1()][i2.index2()], *i2);
      ++i2;
   }
}

template <typename I1, typename I2>
void iter_assign(I1 i1, I2 i2, matrix_iter_row_dense, matrix_iter_coord<injective>)
{
   iter_zero(i1);
   while (i2)
   {
      assign(i1[i2.index1()][i2.index2()], *i2);
      ++i2;
   }
}

// col-major

template <typename I1, typename I2>
void iter_assign(I1 i1, I2 i2, matrix_iter_col_dense, matrix_iter_col_major<sparse>)
{
   iter_zero(i1);
   while (i2)
   {
      add(i1[i2.index2()], *i2);
      ++i1;
      ++i2;
   }
}

template <typename I1, typename I2>
void iter_assign(I1 i1, I2 i2, matrix_iter_col_dense, matrix_iter_col_major<injective>)
{
   iter_zero(i1);
   while (i2)
   {
      assign(i1[i2.index2()], *i2);
      ++i1;
      ++i2;
   }
}


template <typename I1, typename I2>
void iter_assign(I1 i1, I2 i2, matrix_iter_col_dense, matrix_iter_col_major<dense>)
{
   while (i2)
   {
      assign(*i1, *i2);
      ++i1;
      ++i2;
   }
}

// col major, coordinate right-hand-side

template <typename I1, typename I2>
void iter_assign(I1 i1, I2 i2, matrix_iter_col_dense, matrix_iter_coord<sparse>)
{
   iter_zero(i1);
   while (i2)
   {
      add(i1[i2.index2()][i2.index1()], *i2);
      ++i2;
   }
}

template <typename I1, typename I2>
void iter_assign(I1 i1, I2 i2, matrix_iter_col_dense, matrix_iter_coord<injective>)
{
   iter_zero(i1);
   while (i2)
   {
      assign(i1[i2.index2()][i2.index1()], *i2);
      ++i2;
   }
}

// mixed row/col

template <typename I1, typename I2>
void iter_assign(I1 i1, I2 i2, matrix_iter_col_dense, matrix_iter_row_major<sparse>)
{
   iter_assign(swap_sort_order(i1), i2);
}

template <typename I1, typename I2>
void iter_assign(I1 i1, I2 i2, matrix_iter_row_dense, matrix_iter_col_major<sparse>)
{
   iter_assign(swap_sort_order(i1), i2);
}

//
// iter_max
//

template <typename Iter>
inline
typename iterate<Iter>::result_type
iter_matrix_max(Iter I)
{
   CHECK(bool(I));
   typedef typename iterator<Iter>::type Inner;
   Inner MaxIter = iter_max(iterate(I));
   ++I;
   while (I)
   {
      Inner This = iter_max(iterate(I));
      if (*MaxIter < *This) MaxIter = This;
      ++I;
   }
   return MaxIter;
}

//
// iter_min
//

template <typename Iter>
inline
typename iterate<Iter>::result_type
iter_matrix_min(Iter I)
{
   CHECK(bool(I));
   typedef typename iterator<Iter>::type Inner;
   Inner MIter = iter_max(iterate(I));
   ++I;
   while (I)
   {
      Inner This = iter_min(iterate(I));
      if (*This < *MIter) MIter = This;
      ++I;
   }
   return MIter;
}
