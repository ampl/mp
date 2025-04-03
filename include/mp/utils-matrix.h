#ifndef UTILS__MATRIX_H
#define UTILS__MATRIX_H

#include "mp/utils-vec.h"

namespace mp {

/// Triangular matrix.
/// Based on https://github.com/fylux/TriangularMatrix
template<class T, unsigned preallocN=0>
class TMatrix {
private:
  SmallVec<T, (preallocN*preallocN + preallocN)/2> matrix;
  unsigned int N;
protected:
  unsigned int computeSize();
  unsigned int computePosition(unsigned int x, unsigned int y);
public:
  /// Construct
  TMatrix(unsigned int N=0);
  /// Construct with prefill
  TMatrix(unsigned int N, T value);  //Initialize all values
  /// Destruct
  ~TMatrix();

  /// Clear
  void clear();

  /// Resize
  void resize(unsigned int s);

  /// Get
  T get(unsigned int x, unsigned int y);
  /// Set
  void set(unsigned int x, unsigned int y, T value);
  /// Add to existing element
  void add_to(unsigned int x, unsigned int y, T value);

  /// Overloading operator ()
  T operator()(unsigned int x, unsigned int y)
  { return get(x,y); }
  /// Overloading operator ()
  void operator()(
      unsigned int x, unsigned int y, T value)
  { set(x,y,value); }

  //It would be nice to overload '[][]' and '=' operators
};


//Code

template<class T, unsigned pre>
TMatrix<T, pre>::TMatrix(unsigned int N) {
  resize(N);
}

template<class T, unsigned pre>
TMatrix<T, pre>::TMatrix(unsigned int N, T value) {
  resize(N);
  std::fill(matrix.begin(), matrix.end(), value);
}

template<class T, unsigned pre>
TMatrix<T, pre>::~TMatrix() { }

template<class T, unsigned pre>
void TMatrix<T, pre>::clear() {
  matrix.clear();
}

template<class T, unsigned pre>
void TMatrix<T, pre>::resize(unsigned int N) {
  this->N=N;
  matrix.resize(computeSize());
}

template<class T, unsigned pre>
unsigned int TMatrix<T, pre>::computeSize() {
  return (N*N+N)/2;
}

template<class T, unsigned pre>
unsigned int TMatrix<T, pre>::computePosition(
    unsigned int x, unsigned int y) {
  assert(x<N);
  assert(y<N);
  if (y > x) {
    unsigned int aux = x;
    x=y;
    y=aux;
  }
  return (x*x+x)/2+y;
}

template<class T, unsigned pre>
T TMatrix<T, pre>::get(unsigned int x, unsigned int y) {
  return matrix[computePosition(x,y)];
}

template<class T, unsigned pre>
void TMatrix<T, pre>::set(unsigned int x, unsigned int y, T value) {
  matrix[computePosition(x,y)] = value;
}

template<class T, unsigned pre>
void TMatrix<T, pre>::add_to(unsigned int x, unsigned int y, T value) {
  matrix[computePosition(x,y)] += value;
}

}  // namespace mp

#endif // UTILS-MATRIX_H
