#include "matrix/column_major_matrix.h"

#include <cassert>
#include <functional>
#include <random>

#ifdef VERBOSE
#include <iostream>
#endif

template <class T>
Column_Major_Matrix<T>::Column_Major_Matrix(int number_of_row,
                                            int number_of_column,
                                            bool isRand)
    : number_of_row(number_of_row), number_of_column(number_of_column),
      all_column(std::vector<std::vector<T>>(number_of_column,
                                             std::vector<T>(number_of_row))) {
  if (isRand) {
    auto rand = std::bind(std::uniform_int_distribution<int>(-128, 127),
                          std::mt19937());
    for (auto &column : all_column)
      std::generate(begin(column), end(column), rand);
  }
#ifdef VERBOSE
  std::cout << "Column_Major_Matrix::default constructor\n";
#endif
}

template <class T>
Column_Major_Matrix<T>::Column_Major_Matrix(const Column_Major_Matrix &m)
    : number_of_row(m.number_of_row), number_of_column(m.number_of_column),
      all_column(m.all_column) {
#ifdef VERBOSE
  std::cout << "Column_Major_Matrix::copy constructor\n";
#endif
}

template <class T>
Column_Major_Matrix<T>::Column_Major_Matrix(Column_Major_Matrix &&m)
    : number_of_row(std::move(m.number_of_row)),
      number_of_column(std::move(m.number_of_column)),
      all_column(std::move(m.all_column)) {
#ifdef VERBOSE
  std::cout << "Column_Major_Matrix::move constructor\n";
#endif
}

template <class T>
std::vector<T> &Column_Major_Matrix<T>::operator[](int index) {
  return all_column[index];
}

template <class T>
Column_Major_Matrix<T> &
Column_Major_Matrix<T>::operator=(const Column_Major_Matrix &m) {
  number_of_row = m.number_of_row;
  number_of_column = m.number_of_column;
  all_column = m.all_column;
#ifdef VERBOSE
  std::cout << "Column_Major_Matrix::copy assignment\n";
#endif
  return *this;
}

template <class T>
Column_Major_Matrix<T> &
Column_Major_Matrix<T>::operator=(Column_Major_Matrix<T> &&m) {
  number_of_row = std::move(m.number_of_row);
  number_of_column = std::move(m.number_of_column);
  all_column = std::move(m.all_column);
#ifdef VERBOSE
  std::cout << "Column_Major_Matrix::move assignment\n";
#endif
  return *this;
}

template <class T>
Column_Major_Matrix<T>
Column_Major_Matrix<T>::operator*(Row_Major_Matrix<T> &m) {
  assert(number_of_column == m.number_of_row);
  Column_Major_Matrix<T> ret(number_of_row, m.number_of_column, false);
  for (int i = 0; i < number_of_row; i++) {
    for (int j = 0; j < m.number_of_column; j++) {
      T sum = 0;
      for (int k = 0; k < number_of_column; k++)
        sum += all_column[k][i] * m[j][k];
      ret[j][i] = sum;
    }
  }
  return ret;
}

template <class T> Column_Major_Matrix<T>::operator Row_Major_Matrix<T>() {
  Row_Major_Matrix<T> ret(number_of_row, number_of_column);
  for (int i = 0; i < number_of_row; i++)
    for (int j = 0; j < number_of_column; j++)
      ret[i][j] = all_column[j][i];
#ifdef VERBOSE
  std::cout << "col impl\n";
#endif
  return ret;
}

template class Column_Major_Matrix<int>;
