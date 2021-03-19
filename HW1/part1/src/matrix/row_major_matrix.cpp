#include "matrix/row_major_matrix.h"

#include <cassert>
#include <functional>
#include <random>

#ifdef VERBOSE
#include <iostream>
#endif

#include <thread>

template <class T>
Row_Major_Matrix<T>::Row_Major_Matrix(int number_of_row, int number_of_column,
                                      bool isRand)
    : number_of_row(number_of_row), number_of_column(number_of_column),
      all_row(std::vector<std::vector<T>>(number_of_row,
                                          std::vector<T>(number_of_column))) {
  if (isRand) {
    auto rand = std::bind(std::uniform_int_distribution<int>(-128, 127),
                          std::mt19937());
    for (auto &row : all_row)
      std::generate(begin(row), end(row), rand);
  }

#ifdef VERBOSE
  std::cout << "Row_Major_Matrix::default constructor\n";
#endif
}

template <class T>
Row_Major_Matrix<T>::Row_Major_Matrix(const Row_Major_Matrix &m)
    : number_of_row(m.number_of_row), number_of_column(m.number_of_column),
      all_row(m.all_row) {
#ifdef VERBOSE
  std::cout << "Row_Major_Matrix::copy constructor\n";
#endif
}

template <class T>
Row_Major_Matrix<T>::Row_Major_Matrix(Row_Major_Matrix &&m)
    : number_of_row(std::move(m.number_of_row)),
      number_of_column(std::move(m.number_of_column)),
      all_row(std::move(m.all_row)) {
#ifdef VERBOSE
  std::cout << "Row_Major_Matrix::move constructor\n";
#endif
}

template <class T> std::vector<T> &Row_Major_Matrix<T>::operator[](int index) {
  return all_row[index];
}

template <class T>
Row_Major_Matrix<T> &Row_Major_Matrix<T>::operator=(const Row_Major_Matrix &m) {
  number_of_row = m.number_of_row;
  number_of_column = m.number_of_column;
  all_row = m.all_row;
#ifdef VERBOSE
  std::cout << "Row_Major_Matrix::copy assignment\n";
#endif
  return *this;
}

template <class T>
Row_Major_Matrix<T> &Row_Major_Matrix<T>::operator=(Row_Major_Matrix<T> &&m) {
  number_of_row = std::move(m.number_of_row);
  number_of_column = std::move(m.number_of_column);
  all_row = std::move(m.all_row);
#ifdef VERBOSE
  std::cout << "Row_Major_Matrix::move assignment\n";
#endif
  return *this;
}

template <class T>
Row_Major_Matrix<T> Row_Major_Matrix<T>::operator*(Column_Major_Matrix<T> &m) {
  assert(number_of_column == m.number_of_row);
  Row_Major_Matrix<T> ret(number_of_row, m.number_of_column, false);
  for (int i = 0; i < number_of_row; i++) {
    for (int j = 0; j < m.number_of_column; j++) {
      T sum = 0;
      for (int k = 0; k < number_of_column; k++)
        sum += all_row[i][k] * m[j][k];
      ret[i][j] = sum;
    }
  }
  return ret;
}

template <class T> Row_Major_Matrix<T>::operator Column_Major_Matrix<T>() {
  Column_Major_Matrix<T> ret(number_of_row, number_of_column);
  for (int i = 0; i < number_of_row; i++)
    for (int j = 0; j < number_of_column; j++)
      ret[j][i] = all_row[i][j];
#ifdef VERBOSE
  std::cout << "row impl\n";
#endif
  return ret;
}

template <class T>
void partial_matrix_multiply(Row_Major_Matrix<T> &r, int rbg, int red,
                             Column_Major_Matrix<T> &c, int cbg, int ced,
                             int inner, Row_Major_Matrix<T> &ret) {
#ifdef VERBOSE
  std::cout << rbg << ' ' << red << ' ' << cbg << ' ' << ced << '\n';
#endif
  for (int i = rbg; i < red; i++) {
    for (int j = cbg; j < ced; j++) {
      T sum = 0;
      for (int k = 0; k < inner; k++)
        sum += r[i][k] * c[j][k];
      ret[i][j] = sum;
    }
  }
}

template <class T>
Row_Major_Matrix<T> Row_Major_Matrix<T>::operator%(Column_Major_Matrix<T> &m) {
  assert(number_of_column == m.number_of_row);
  int number_of_thread = 4;

  std::vector<int> range{0};
  for (int i = 0; i < number_of_thread; i++)
    range.push_back(range.back() + number_of_row / number_of_thread +
                    (i < number_of_row % number_of_thread));

  Row_Major_Matrix<T> ret(number_of_row, m.number_of_column);

  std::vector<std::thread> threads;
  for (int i = 0; i < number_of_thread; i++)
    threads.push_back(std::thread(
        partial_matrix_multiply<T>, std::ref(*this), range[i], range[i + 1],
        std::ref(m), 0, m.number_of_column, number_of_column, std::ref(ret)));
  for (auto &t : threads)
    t.join();
  return ret;
}

template class Row_Major_Matrix<int>;
