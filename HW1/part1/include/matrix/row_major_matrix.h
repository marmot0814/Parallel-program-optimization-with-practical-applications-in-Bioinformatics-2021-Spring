#ifndef __MATRIX_ROW_MAJOR_MATRIX_H__
#define __MATRIX_ROW_MAJOR_MATRIX_H__

#include <vector>
#include "matrix/column_major_matrix.h"

template <class T>
class Column_Major_Matrix;

template <class T>
class Row_Major_Matrix {
  template <class U> friend class Column_Major_Matrix;
  int number_of_row;
  int number_of_column;
  std::vector<std::vector<T>> all_row;
 public:
  Row_Major_Matrix(int, int, bool isRand = true);
  Row_Major_Matrix(const Row_Major_Matrix&);
  Row_Major_Matrix(Row_Major_Matrix&&);

  std::vector<T>& operator[](int);
  Row_Major_Matrix& operator=(const Row_Major_Matrix &);
  Row_Major_Matrix& operator=(Row_Major_Matrix &&);
  Row_Major_Matrix operator*(Column_Major_Matrix<T> &m);
  operator Column_Major_Matrix<T>();

  Row_Major_Matrix operator%(Column_Major_Matrix<T> &);
};

#endif  //  __MATRIX_ROW_MAJOR_MATRIX_H__
