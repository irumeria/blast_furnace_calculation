
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

template <typename T>
void multable(vector<vector<T>> &matrix_a, vector<vector<T>> &matrix_b) {
  if (matrix_a.size() != matrix_b[0].size()) {
    //
  }
  if (matrix_a[0].size() != matrix_b.size()) {
    //
  }
}

template <typename T>
vector<vector<T>> mDot(vector<vector<T>> &matrix_a,
                       vector<vector<T>> &matrix_b) {
  // multable(matrix_a, matrix_b);
  vector<vector<T>> ret = vector<vector<T>>(matrix_a.size(), matrix_a[0]);
  for (auto i = 0; i < matrix_a.size(); i++) {
    for (auto j = 0; j < matrix_a[0].size(); j++) {
      ret[i][j] = matrix_a[i][j] * matrix_b[j][i];
    }
  }
  return ret;
}

template <typename T> vector<vector<T>> transpose(vector<vector<T>> &matrix) {
  vector<vector<T>> ret = vector<vector<T>>(matrix.size(), matrix[0]);
  for (auto i = 0; i < matrix.size(); i++) {
    for (auto j = 0; j < matrix[0].size(); j++) {
      ret[i][j] = matrix[j][i];
    }
  }
  return ret;
}

template <typename T>
vector<vector<T>> mCross(vector<vector<T>> &matrix_a,
                         vector<vector<T>> &matrix_b) {

  // std::cout << matrix_a.size() << " " << matrix_b[0].size() << "\n";
  // multable(matrix_a, transpose(matrix_b));

  vector<vector<T>> ret = vector<vector<T>>(matrix_a.size(), matrix_b[0]);
  for (size_t i = 0; i < matrix_a.size(); i++) {
    for (size_t j = 0; j < matrix_b[0].size(); j++) {
      T element = {};
      for (size_t k = 0; k < matrix_a[0].size(); k++) {
        element += matrix_a[i][k] * matrix_b[k][j];
      }
      ret[i][j] = element;
    }
  }
  return ret;
}

// 按第一行展开计算A的行列式
template <typename T> T getA(vector<vector<T>> &arcs, int n) {
  if (n == 1) {
    return arcs[0][0];
  }
  T ans = {};
  vector<vector<T>> temp = vector<vector<T>>(arcs.size(), arcs[0]);
  int i, j, k;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n - 1; j++) {
      for (k = 0; k < n - 1; k++) {
        temp[j][k] = arcs[j + 1][(k >= i) ? k + 1 : k];
      }
    }
    double t = getA(temp, n - 1);
    if (i % 2 == 0) {
      ans += arcs[0][i] * t;
    } else {
      ans -= arcs[0][i] * t;
    }
  }
  return ans;
}

// 计算每一行每一列的每个元素所对应的余子式，组成A^*
template <typename T>
void getAStart(vector<vector<T>> &arcs, int n, vector<vector<T>> &ans) {
  if (n == 1) {
    ans[0][0] = 1;
    return;
  }
  int i, j, k, t;
  vector<vector<T>> temp = vector<vector<T>>(arcs.size(), arcs[0]);
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      for (k = 0; k < n - 1; k++) {
        for (t = 0; t < n - 1; t++) {
          temp[k][t] = arcs[k >= i ? k + 1 : k][t >= j ? t + 1 : t];
        }
      }

      ans[j][i] = getA(temp, n - 1); //此处顺便进行了转置
      if ((i + j) % 2 == 1) {
        ans[j][i] = -ans[j][i];
      }
    }
  }
}

// 得到给定矩阵src的逆矩阵保存到des中。
template <typename T>
bool get_M_Inverse(vector<vector<T>> &src, vector<vector<T>> &des) {
  auto n = src.size();
  T flag = getA(src, n);
  vector<vector<T>> t = vector<vector<T>>(src.size(), src[0]);
  if (0 == flag) {
    cerr << "err: |A| = 0" << endl;
    return false;
  } else {
    getAStart(src, n, t);
    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < n; j++) {
        des[i][j] = t[i][j] / flag;
      }
    }
  }

  return true;
}

// 打印矩阵
template <typename T> void printMatrix(vector<vector<T>> &matrix, int reserve_number) {
  for (size_t i = 0; i < matrix.size(); i++) {
    for (size_t j = 0; j < matrix.size(); j++) {
      cout << setprecision(reserve_number) << matrix[i][j] << "\t";
    }
    cout << "\n";
  }
}
