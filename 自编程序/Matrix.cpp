#include "Matrix.h"

using namespace Matrix_hx;

/*
**********************************************************
函数名：矩阵初始化函数
参数：rows    矩阵行数
      cols	  矩阵列数
      value   矩阵每个元素的初始值
函数功能：构造并初始化矩阵
**********************************************************
*/
Matrix::Matrix(size_t rows, size_t cols, double value = 0.0)
{
    data.resize(rows, vector<double>(cols, value));
}

/*
**********************************************************
函数名：矩阵行数
返回值：矩阵的行数
函数功能：获得矩阵的行数
**********************************************************
*/
size_t Matrix::rows() const
{
    return data.size();
}

/*
**********************************************************
函数名：矩阵列数
返回值：矩阵的列数
函数功能：获得矩阵的列数
**********************************************************
*/
size_t Matrix::cols() const
{
    return data.empty() ? 0 : data[0].size();
}

/*
**********************************************************
函数名：打印矩阵
函数功能：将矩阵打印输出到屏幕
**********************************************************
*/
void Matrix::printMatrix() const
{
    for (const auto& row : data)
    {
        for (const auto& elem : row)
        {
            cout << elem << " ";
        }
        cout << std::endl;
    }
}

/*
**********************************************************
函数名：重载()
参数：row    矩阵某一行的索引（第1行的索引是0）
      col	 矩阵某一列的索引（第1列的索引是0）
返回值：矩阵第row+1行第col+1列的元素值
函数功能：重载运算符()，通过行列索引获取矩阵某一元素的值，且可改变
**********************************************************
*/
double& Matrix::operator()(size_t row, size_t col)
{
    if (row >= rows() || col >= cols())
    {
        throw out_of_range("Index out of range");
    }
    return data[row][col];
}

/*
**********************************************************
函数名：重载()
参数：row    矩阵某一行的索引（第1行的索引是0）
      col	 矩阵某一列的索引（第1列的索引是0）
返回值：矩阵第row+1行第col+1列的元素值
函数功能：重载运算符()，通过行列索引获取矩阵某一元素的值，但不能改变
**********************************************************
*/
const double& Matrix::operator()(size_t row, size_t col) const
{
    if (row >= rows() || col >= cols())
    {
        throw out_of_range("Index out of range");
    }
    return data[row][col];
}

/*
**********************************************************
函数名：重载+
参数：other    另一个矩阵
返回值：两个矩阵相加得到的矩阵
函数功能：重载运算符+，实现矩阵的相加
**********************************************************
*/
Matrix Matrix::operator+(const Matrix& other) const
{
    if (rows() != other.rows() || cols() != other.cols())
    {
        throw invalid_argument("Matrices dimensions must match");
    }
    Matrix result(rows(), cols(), 0);
    for (size_t i = 0; i < rows(); ++i)
    {
        for (size_t j = 0; j < cols(); ++j)
        {
            result(i, j) = (*this)(i, j) + other(i, j);
        }
    }
    return result;
}

/*
**********************************************************
函数名：重载-
参数：other    另一个矩阵(减数)
返回值：两个矩阵相减得到的矩阵
函数功能：重载运算符-，实现矩阵的相减
**********************************************************
*/
Matrix Matrix::operator-(const Matrix& other) const
{
    if (rows() != other.rows() || cols() != other.cols())
    {
        throw invalid_argument("Matrices dimensions must match");
    }
    Matrix result(rows(), cols(), 0);
    for (size_t i = 0; i < rows(); ++i)
    {
        for (size_t j = 0; j < cols(); ++j)
        {
            result(i, j) = (*this)(i, j) - other(i, j);
        }
    }
    return result;
}

/*
**********************************************************
函数名：重载*
参数：other    乘号右边的矩阵
返回值：两个矩阵相乘得到的矩阵
函数功能：重载运算符*，实现矩阵的相乘
**********************************************************
*/
Matrix Matrix::operator*(const Matrix& other) const
{
    if (cols() != other.rows())
    {
        throw invalid_argument("Matrix A's columns must be equal to Matrix B's rows");
    }
    Matrix result(rows(), other.cols(), 0);
    for (size_t i = 0; i < result.rows(); ++i)
    {
        for (size_t j = 0; j < result.cols(); ++j)
        {
            result(i, j) = 0.0;
            for (size_t k = 0; k < cols(); ++k)
            {
                result(i, j) += (*this)(i, k) * other(k, j);
            }
        }
    }
    return result;
}

/*
**********************************************************
函数名：重载*
参数：num    常数
返回值：一个矩阵乘以一个常数
函数功能：重载运算符*，实现矩阵的相乘
**********************************************************
*/
Matrix Matrix::operator*(const double num) const
{
    Matrix result(rows(), cols(), 0);
    for (size_t i = 0; i < result.rows(); ++i)
    {
        for (size_t j = 0; j < result.cols(); ++j)
        {
            result(i, j) = (*this)(i, j) * num;
        }
    }
    return result;
}

/*
**********************************************************
函数名：转秩
返回值：转秩后的矩阵
函数功能：实现矩阵的转秩
**********************************************************
*/
Matrix Matrix::transpose() const
{
    Matrix result(cols(), rows(), 0);
    for (size_t i = 0; i < rows(); ++i)
    {
        for (size_t j = 0; j < cols(); ++j)
        {
            result(j, i) = (*this)(i, j);
        }
    }
    return result;
}

/*
**********************************************************
函数名：交换行
参数：row1    矩阵某一行的索引
      row2    矩阵另一行的索引
函数功能：交换矩阵的某两行
**********************************************************
*/
void Matrix::swapRows(size_t row1, size_t row2)
{
    swap_ranges(data[row1].begin(), data[row1].end(), data[row2].begin());
}

/*
**********************************************************
函数名：将行乘以一个标量
参数：row     矩阵某一行的索引
      scalar  待乘的标量
函数功能：将矩阵的某一行乘以一个标量
**********************************************************
*/
void Matrix::multiplyRow(size_t row, double scalar)
{
    for (auto& elem : data[row])
    {
        elem *= scalar;
    }
}

/*
**********************************************************
函数名：一行加上另一行的倍数
参数：fromrow    加行
      torow      被改变的行
      scalar     倍数
函数功能：将某一行乘以一个倍数加到另一行上，改变另一行的值
**********************************************************
*/
void Matrix::addRowToRow(size_t fromRow, size_t toRow, double scalar)
{
    for (size_t col = 0; col < data[0].size(); ++col)
    {
        data[toRow][col] += data[fromRow][col] * scalar;
    }
}

/*
**********************************************************
函数名：矩阵求逆
返回值：矩阵的逆矩阵
函数功能：使用高斯消元法，获得矩阵的逆矩阵
**********************************************************
*/
Matrix Matrix::inverse() const
{
    if (rows() != cols())
    {
        throw invalid_argument("Matrix must be square for inversion");
    }
    Matrix A = *this;//创建矩阵的副本，以免改变原矩阵
    size_t n = A.rows();
    if (n <= 0)
    {
        throw invalid_argument("Matrix's rows or cols must more than 0 ");
    }
    Matrix inverse(n, n, 0.0);

    // 创建单位矩阵
    for (size_t i = 0; i < n; ++i)
    {
        inverse(i, i) = 1.0;
    }

    // 高斯-约旦消元法
    for (size_t i = 0; i < n; ++i)
    {
        // 寻找最大的行元素
        size_t maxRow = i;
        for (size_t k = i + 1; k < n; ++k)
        {
            if (abs(A(k, i)) > abs(A(maxRow, i)))
            {
                maxRow = k;
            }
        }
        if (A(maxRow, i) == 0)
        {
            throw runtime_error("Matrix is singular and cannot be inverted");
        }
        if (maxRow != i)
        {
            A.swapRows(i, maxRow);
            inverse.swapRows(i, maxRow);
        }
        //将对角元素规范化为1，并更新逆矩阵
        double div = A(i, i);
        A.multiplyRow(i, 1.0 / div);
        inverse.multiplyRow(i, 1.0 / div);

        for (size_t j = 0; j < n; ++j)
        {
            if (i != j)
            {
                double scalar = -A(j, i);
                A.addRowToRow(i, j, scalar);
                inverse.addRowToRow(i, j, scalar);
            }
        }
    }

    return inverse;
}

/*
**********************************************************
函数名：删除行
参数：rowIndex    矩阵某一行的索引
函数功能：将矩阵的某一行删除，从而减少矩阵的行数
**********************************************************
*/
void Matrix::deletRow(int rowIndex)
{
    if (rowIndex >= 0 && rowIndex < data.size())
    {
        data.erase(data.begin() + rowIndex);
    }
}

/*
**********************************************************
函数名：删除列
参数：colIndex    矩阵某一列的索引
函数功能：将矩阵的某一列删除，从而减少矩阵的列数
**********************************************************
*/
void Matrix::deleteColumn(int colIndex)
{
    if (colIndex >= 0 && colIndex < data[0].size())
    {
    for (auto& row : data)
    {
        row.erase(row.begin() + colIndex);
    }
    }
}

/*
**********************************************************
函数名：重载^
参数：other    另一个三维向量
返回值：两个三维向量叉乘得到的三维向量
函数功能：重载运算符~，实现三维向量的叉乘
**********************************************************
*/
Matrix Matrix::operator^(const Matrix& other) const
{
    if (rows() != 3 || cols() != 1 || other.rows() != 3 || other.cols() != 1)
    {
        throw invalid_argument("only tri-vector can calculate cross multiplication.");
    }
    Matrix result(3, 1, 0);
    result(0, 0) = (*this)(1, 0) * other(2, 0) - (*this)(2, 0) * other(1, 0);
    result(1, 0) = (*this)(2, 0) * other(0, 0) - (*this)(0, 0) * other(2, 0);
    result(2, 0) = (*this)(0, 0) * other(1, 0) - (*this)(1, 0) * other(0, 0);

    return result;
}


/*
**********************************************************
函数名：计算向量的模
返回值：当前向量的模
函数功能：计算当前向量的模并返回
**********************************************************
*/
double Matrix::norm() const
{
    if ((*this).cols() != 1)
    {
        throw invalid_argument("Matrices dimensions must match");
    }
    double norm = 0.0;
    for (int i = 0; i < (*this).rows(); i++)
    {
        norm += (*this)(i, 0) * (*this)(i, 0);
    }
    norm = sqrt(norm);

    return norm;
}

/*
**********************************************************
函数名：重载/
参数：other    另一个矩阵
返回值：当前矩阵与另一个矩阵拼接后的矩阵
函数功能：重载运算符/，实现矩阵的纵向拼接
**********************************************************
*/
Matrix Matrix::operator/(const Matrix& other) const
{
    if ((*this).cols() != other.cols())
    {
        throw invalid_argument("Matrices dimensions must match");
    }
    Matrix result((*this).rows() + other.rows(), other.cols(), 0.0);
    for (int i = 0; i < (*this).rows(); i++)
    {
        for (int j = 0; j < (*this).cols(); j++)
        {
            result(i, j) = (*this)(i, j);
        }
    }

    for (int i = 0; i < other.rows(); i++)
    {
        for (int j = 0; j < other.cols(); j++)
        {
            result((*this).rows() + i, j) = other(i, j);
        }
    }

    return result;
}


/*
**********************************************************
函数名：重载，
参数：other    另一个矩阵
返回值：当前矩阵与另一个矩阵拼接后的矩阵
函数功能：重载运算符“，”，实现矩阵的横向拼接
**********************************************************
*/
Matrix Matrix::operator,(const Matrix& other) const
{
    if ((*this).rows() != other.rows())
    {
        throw invalid_argument("Matrices dimensions must match");
    }

    Matrix result(other.rows(), (*this).cols() + other.cols(), 0.0);
    for (int i = 0; i < other.rows(); i++)
    {
        for (int j = 0; j < (*this).cols(); j++)
        {
            result(i, j) = (*this)(i, j);
        }
    }

    for (int i = 0; i < other.rows(); i++)
    {
        for (int j = 0; j < other.cols(); j++)
        {
            result(i, (*this).cols() + j) = other(i, j);
        }
    }

    return result;
}


/*
**********************************************************
函数名：重载&
参数：other    另一个四元数
返回值：当前四元数与另一个四元数相乘后得到的四元数
函数功能：重载运算符&，实现四元数乘法
**********************************************************
*/
Matrix Matrix::operator&(const Matrix& other) const
{
    if ((*this).rows() != 4 || (*this).cols() != 1 || other.rows() != 4 || other.cols() != 1)
    {
        throw invalid_argument("Matrices dimensions must match");
    }
    Matrix result(4, 1, 0.0);
    result(0, 0) = (*this)(0, 0) * other(0, 0) - (*this)(1, 0) * other(1, 0) - (*this)(2, 0) * other(2, 0) - (*this)(3, 0) * other(3, 0);
    result(1, 0) = (*this)(0, 0) * other(1, 0) + (*this)(1, 0) * other(0, 0) + (*this)(2, 0) * other(3, 0) - (*this)(3, 0) * other(2, 0);
    result(2, 0) = (*this)(0, 0) * other(2, 0) + (*this)(2, 0) * other(0, 0) + (*this)(3, 0) * other(1, 0) - (*this)(1, 0) * other(3, 0);
    result(3, 0) = (*this)(0, 0) * other(3, 0) + (*this)(3, 0) * other(0, 0) + (*this)(1, 0) * other(2, 0) - (*this)(2, 0) * other(1, 0);

    return result;
}

/*
**********************************************************
函数名：生成单位阵
返回值：当前矩阵纬度的单位阵
函数功能：按照当前矩阵的纬度生成一个单位阵
**********************************************************
*/
Matrix Matrix::eye() const
{
    if (rows() != cols())
    {
        throw invalid_argument("Matrices dimensions must match");
    }
    Matrix result(rows(), cols(), 0.0);
    for (int i = 0; i < rows(); i++)
    {
        for (int j = 0; j < cols(); j++)
        {
            if (i != j) result(i, j) = 0.0;
            else if (i == j) result(i, j) = 1.0;
        }
    }
    return result;
}

/*
**********************************************************
函数名：反对称矩阵
返回值：当前三维向量的反对称矩阵
函数功能：将当前的三位向量转换为反对称矩阵
**********************************************************
*/
Matrix Matrix::skew() const
{
    if (rows() != 3 || cols() != 1)
    {
        throw invalid_argument("Matrices dimensions must match");
    }
    Matrix result(3, 3, 0.0);
    double vx = (*this)(0, 0);
    double vy = (*this)(1, 0);
    double vz = (*this)(2, 0);
    result(0, 1) = -vz;
    result(0, 2) = vy;
    result(1, 0) = vz;
    result(1, 2) = -vx;
    result(2, 0) = -vy;
    result(2, 1) = vx;

    return result;
}

/*
**********************************************************
函数名：对角阵
返回值：以当前向量为对角元素的对角阵
函数功能：以当前向量为对角元素，生成对角阵
**********************************************************
*/
Matrix Matrix::diag() const
{
    if (cols() != 1)
    {
        throw invalid_argument("only vector can calculate diag.");
    }
    Matrix result(rows(), rows(), 0.0);
    for (int i = 0; i < rows(); i++)
    {
        result(i, i) = (*this)(i, 0);
    }

    return result;
}

/*
**********************************************************
函数名：四元数相乘后单位化
参  数：q2        另一个四元数
返回值：当前四元数与另一个四元数相乘并单位化后得到的四元数
函数功能：将当前四元数与另一个四元数相乘并单位化
**********************************************************
*/
Matrix Matrix::quatProd(const Matrix q2) const
{
    if (rows() != 4 || cols() != 1 || q2.rows() != 4 || q2.cols() != 1)
    {
        throw invalid_argument("Matrices dimensions must match");
    }
    Matrix q = (*this) & q2;

    if (q(0,0) < 0)
    {
        q = q * (-1);
    }
    double norm = q.norm();
    for (int i = 0; i < 4; i++)
        q(i, 0) /= norm;

    return q;
}