#include "s21_matrix_oop.h"

// КОНСТРУКТОРЫ И ДЕСТРУКТОР

// Базовый конструктор, инициализирующий матрицу некоторой заранее заданной
// размерностью
S21Matrix::S21Matrix() {
  rows_ = 3;
  cols_ = 3;

  matrix_ = new double *[rows_];
  for (int i = 0; i < rows_; i++) matrix_[i] = new double[cols_];

  FillMatrixByZero();
}

// Параметризированный конструктор с количеством строк и столбцов
S21Matrix::S21Matrix(int rows, int cols) {
  if (rows <= 0 || cols <= 0) {
    throw std::out_of_range("Incorrect input, index is out of range.");
  } else {
    rows_ = rows;
    cols_ = cols;

    matrix_ = new double *[rows_];
    for (int i = 0; i < rows_; i++) matrix_[i] = new double[cols_];

    FillMatrixByZero();
  }
}

// Конструктор копирования
S21Matrix::S21Matrix(const S21Matrix &other) { CopyMatrix(other); }

// Конструктор переноса
S21Matrix::S21Matrix(S21Matrix &&other) noexcept {
  matrix_ = other.matrix_;
  rows_ = other.rows_;
  cols_ = other.cols_;

  other.rows_ = 0;
  other.cols_ = 0;
}

// Деструктор
S21Matrix::~S21Matrix() { Free(); }

// МЕТОДЫ

// Сравнивает матрицы на равенство

bool S21Matrix::EqMatrix(const S21Matrix &other) {
  bool result = true;
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    result = false;
  } else {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        if (matrix_[i][j] != other.matrix_[i][j]) {
          result = false;
          break;
        }
      }
    }
  }
  return result;
}

// Складывает матрицы
void S21Matrix::SumMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::out_of_range(
        "Incorrect input: matriсes should have the same "
        "size if you want to sum them.");
  } else {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        matrix_[i][j] = matrix_[i][j] + other.matrix_[i][j];
      }
    }
  }
}

// Вычитает матрицы
void S21Matrix::SubMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::out_of_range(
        "Incorrect input: matriсes should have the same "
        "size if you want to sum them.");
  } else {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        matrix_[i][j] = matrix_[i][j] - other.matrix_[i][j];
      }
    }
  }
}

// Умножает матрицу на число
void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = matrix_[i][j] * num;
    }
  }
}

// Умножает матрицы
void S21Matrix::MulMatrix(const S21Matrix &other) {
  if (cols_ != other.rows_) {
    throw std::out_of_range(
        "Incorrect input. Number of first matrix columns "
        "should be equal for second matrix rows");
  } else {
    S21Matrix result(this->rows_, other.cols_);
    for (int i = 0; i < result.rows_; i++) {
      for (int j = 0; j < result.cols_; j++) {
        for (int k = 0; k < this->cols_; k++) {
          result.matrix_[i][j] += (this->matrix_[i][k] * other.matrix_[k][j]);
        }
      }
    }
    this->Free();
    *this = std::move(result);
  }
}

// Создает новую транспонированную матрицу из текущей и возвращает её
S21Matrix S21Matrix::Transpose() {
  S21Matrix result(cols_, rows_);

  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) result.matrix_[j][i] = matrix_[i][j];

  return result;
}

// Вычисляет матрицу алгебраических дополнений текущей матрицы и возвращает ее
S21Matrix S21Matrix::CalcComplements() {
  S21Matrix result(rows_, cols_);
  double det = Determinant();

  if (rows_ != cols_) {
    throw std::out_of_range(
        "Incorrect input: you can't get a complements matrix for matrix that "
        "is not square.");
  } else if (det == 0) {
    throw std::out_of_range(
        "Incorrect input: you can't get a complements matrix for matrix which "
        "determinant is zero.");
  } else {
    double minor_det = 0;

    for (int i = 0; i < rows_; i++)
      for (int j = 0; j < cols_; j++) {
        S21Matrix minor = GetMinor(i, j);
        minor_det = minor.Determinant();
        minor.Free();
        result.matrix_[i][j] = minor_det * pow(-1, i + j);
      }
  }

  return result;
}

// Вычисляет и возвращает определитель текущей матрицы
double S21Matrix::Determinant() {
  double result;

  if (rows_ != cols_) {
    throw std::out_of_range(
        "Incorrect input: you can't get a determinant for matrix that is not "
        "square.");
  } else if (rows_ == 1) {
    result = matrix_[0][0];
  } else if (rows_ == 2) {
    result = (matrix_[0][0] * matrix_[1][1]) - (matrix_[1][0] * matrix_[0][1]);
  } else {
    double minor_det = 0;
    result = 0;

    for (int x = 0; x < cols_; x++) {
      // Идём по первой строке, четвёртый параметр - вычёркиваемый столбец.
      S21Matrix minor = GetMinor(0, x);
      minor_det = minor.Determinant();
      minor.Free();

      if (x % 2 == 0)
        result += (matrix_[0][x] * minor_det);
      else
        result -= (matrix_[0][x] * minor_det);
    }
  }

  return result;
}

// Вычисляет и возвращает обратную матрицу
S21Matrix S21Matrix::InverseMatrix() {
  S21Matrix result(rows_, cols_);

  if (rows_ != cols_) {
    std::out_of_range(
        "Incorrect input: you can't inverse matrix that is not square.");
  } else {
    double det = Determinant();
    if (det == 0) {
      std::out_of_range(
          "Incorrect input: you can't inverse matrix if it's determinant is "
          "equal to zero.");
    } else {
      S21Matrix complem = CalcComplements();
      S21Matrix transp = complem.Transpose();
      transp.MulNumber(1.0 / det);

      result = std::move(transp);

      complem.Free();
      transp.Free();
    }
  }

  return result;
}

// ПЕРЕГРУЗКА ОПЕРАТОРОВ

// Сложение двух матриц
S21Matrix S21Matrix::operator+(const S21Matrix &other) {
  S21Matrix result = *this;
  result.SumMatrix(other);

  return result;
}

// Вычитание одной матрицы из другой
S21Matrix S21Matrix::operator-(const S21Matrix &other) {
  S21Matrix result = *this;
  result.SubMatrix(other);

  return result;
}

// Умножение матриц
S21Matrix S21Matrix::operator*(const S21Matrix &other) {
  S21Matrix result = *this;
  result.MulMatrix(other);

  return result;
}

// Умножение матрицы на число
S21Matrix S21Matrix::operator*(const double num) {
  S21Matrix result = *this;
  result.MulNumber(num);

  return result;
}

// Проверка на равенство матриц (EqMatrix)
bool S21Matrix::operator==(const S21Matrix &other) {
  return (this->EqMatrix(other));
}

// Присвоение матрице значений другой матрицы
S21Matrix &S21Matrix::operator=(const S21Matrix &other) {
  if (!(*this == other)) {
    Free();
    CopyMatrix(other);
  }

  return *this;
}

// Присвоение сложения (SumMatrix)
S21Matrix &S21Matrix::operator+=(const S21Matrix &other) {
  this->SumMatrix(other);

  return *this;
}

// Присвоение разности (SubMatrix)
S21Matrix &S21Matrix::operator-=(const S21Matrix &other) {
  this->SubMatrix(other);

  return *this;
}

// Присвоение умножения (MulMatrix)
S21Matrix &S21Matrix::operator*=(const S21Matrix &other) {
  this->MulMatrix(other);

  return *this;
}

// Присвоение умножения (MulNumber)
S21Matrix &S21Matrix::operator*=(const double num) {
  this->MulNumber(num);

  return *this;
}

// Индексация по элементам матрицы (строка, колонка)
double S21Matrix::operator()(int i, int j) {
  if (i >= rows_ || j >= cols_ || i < 0 || j < 0) {
    throw std::out_of_range("Incorrect input, index is out of range.");
  }

  return matrix_[i][j];
}

// Accessors и mutators

int S21Matrix::GetRows() { return rows_; }

int S21Matrix::GetCols() { return cols_; }

void S21Matrix::SetValue(int row, int col, double value) {
  if (row >= rows_ || col >= cols_ || row < 0 || col < 0) {
    throw std::out_of_range("Incorrect input, index is out of range.");
  } else {
    matrix_[row][col] = value;
  }
}

// При увеличении размера - матрица дополняется нулевыми элементами, при
// уменьшении - лишнее просто отбрасывается
void S21Matrix::EditSize(int rows, int cols) {
  if (rows == rows_ && cols == cols_) {
    // do nothing
  } else if (rows <= 0 || cols <= 0) {
    throw std::out_of_range("Incorrect input, index is out of range.");
  } else {
    S21Matrix result(rows, cols);

    if (rows <= rows_ && cols <= cols_) {
      for (int y = 0; y < rows; y++)
        for (int x = 0; x < cols; x++) result.SetValue(y, x, matrix_[y][x]);
    } else if (rows > rows_ && cols > cols_) {
      for (int y = 0; y < rows_; y++)
        for (int x = 0; x < cols_; x++) result.SetValue(y, x, matrix_[y][x]);
    } else if (rows <= rows_ && cols > cols_) {
      for (int y = 0; y < rows; y++)
        for (int x = 0; x < cols_; x++) result.SetValue(y, x, matrix_[y][x]);
    } else if (rows > rows_ && cols <= cols_) {
      for (int y = 0; y < rows_; y++)
        for (int x = 0; x < cols; x++) result.SetValue(y, x, matrix_[y][x]);
    }

    Free();
    *this = std::move(result);
  }
}

S21Matrix S21Matrix::SetCols(int cols) {
  S21Matrix temp(rows_, cols);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < ((cols_ < cols) ? cols_ : cols); j++) {
      temp.matrix_[i][j] = matrix_[i][j];
    }
  }
  *this = std::move(temp);
  return *this;
}

S21Matrix S21Matrix::SetRows(int rows) {
  S21Matrix temp(rows, cols_);
  for (int i = 0; i < ((rows_ < rows) ? rows_ : rows); i++) {
    for (int j = 0; j < cols_; j++) {
      temp.matrix_[i][j] = matrix_[i][j];
    }
  }
  *this = std::move(temp);
  return *this;
}
// Вспомогательные

// Заполняет матрицу нулями
void S21Matrix::FillMatrixByZero() {
  for (int y = 0; y < rows_; y++)
    for (int x = 0; x < cols_; x++) matrix_[y][x] = (double)0;
}

// Заполняет матрицу значениями от 0 до (rows_ * cols_ - 1)
void S21Matrix::FillMatrix() {
  int k = 0;
  for (int y = 0; y < rows_; y++)
    for (int x = 0; x < cols_; x++) matrix_[y][x] = (double)k++;
}

// Заполняет матрицу случайными значениями
void S21Matrix::FillMatrixRandom() {
  for (int y = 0; y < rows_; y++)
    for (int x = 0; x < cols_; x++) matrix_[y][x] = (double)(std::rand() % 100);
}

// Выводит матрицу в консоль
void S21Matrix::PrintMatrix() {
  printf("matrix:\n");
  for (int y = 0; y < rows_; y++) {
    for (int x = 0; x < cols_; x++) {
      printf("%5.0lf ", matrix_[y][x]);
    }
    printf("\n");
  }
}

// Освобождает память, выделенную под матрицы
void S21Matrix::Free() {
  if (rows_ != 0 && cols_ != 0) {
    if (matrix_) {
      for (int i = 0; i < rows_; i++) {
        if (matrix_[i]) delete matrix_[i];
      }
      delete matrix_;
    }
  }
  rows_ = 0;
  cols_ = 0;
}

// Возвращает минор матрицы
S21Matrix S21Matrix::GetMinor(int oy, int ox) {
  S21Matrix result(rows_ - 1, cols_ - 1);

  for (int y = 0; y < rows_; y++)
    for (int x = 0; x < cols_; x++) {
      if (y == oy || x == ox)
        continue;
      else if (y < oy && x < ox)
        result.matrix_[y][x] = matrix_[y][x];
      else if (y > oy && x < ox)
        result.matrix_[y - 1][x] = matrix_[y][x];
      else if (y < oy && x > ox)
        result.matrix_[y][x - 1] = matrix_[y][x];
      else if (y > oy && x > ox)
        result.matrix_[y - 1][x - 1] = matrix_[y][x];
    }

  return result;
}

// Копирует содержимое одной матрицы в другую
void S21Matrix::CopyMatrix(const S21Matrix &other) {
  rows_ = other.rows_;
  cols_ = other.cols_;
  matrix_ = new double *[rows_];
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_];
    std::memcpy(matrix_[i], other.matrix_[i], other.cols_ * sizeof(double));
  }
}