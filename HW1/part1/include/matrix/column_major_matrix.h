#ifndef __MATRIX_COLUMN_MAJOR_MATRIX_H__
#define __MATRIX_COLUMN_MAJOR_MATRIX_H__

#include <vector>
#include "matrix/row_major_matrix.h"

template <class T>
class Row_Major_Matrix;

template <class T>
class Column_Major_Matrix {
  template <class U> friend class Row_Major_Matrix;
  int number_of_row;
  int number_of_column;
  std::vector<std::vector<T>> all_column;
 public:
  Column_Major_Matrix(int, int, bool isRand = true);
  Column_Major_Matrix(const Column_Major_Matrix&);
  Column_Major_Matrix(Column_Major_Matrix&&);

  std::vector<T>& operator[](int);
  Column_Major_Matrix& operator=(const Column_Major_Matrix &);
  Column_Major_Matrix& operator=(Column_Major_Matrix &&);
  Column_Major_Matrix operator*(Row_Major_Matrix<T> &);
  operator Row_Major_Matrix<T>();

  Column_Major_Matrix operator%(Row_Major_Matrix<T> &);
};

#endif  //  __MATRIX_COLUMN_MAJOR_MATRIX_H__
