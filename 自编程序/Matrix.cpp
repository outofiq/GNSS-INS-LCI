#include "Matrix.h"

using namespace Matrix_hx;

/*
**********************************************************
�������������ʼ������
������rows    ��������
      cols	  ��������
      value   ����ÿ��Ԫ�صĳ�ʼֵ
�������ܣ����첢��ʼ������
**********************************************************
*/
Matrix::Matrix(size_t rows, size_t cols, double value = 0.0)
{
    data.resize(rows, vector<double>(cols, value));
}

/*
**********************************************************
����������������
����ֵ�����������
�������ܣ���þ��������
**********************************************************
*/
size_t Matrix::rows() const
{
    return data.size();
}

/*
**********************************************************
����������������
����ֵ�����������
�������ܣ���þ��������
**********************************************************
*/
size_t Matrix::cols() const
{
    return data.empty() ? 0 : data[0].size();
}

/*
**********************************************************
����������ӡ����
�������ܣ��������ӡ�������Ļ
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
������������()
������row    ����ĳһ�е���������1�е�������0��
      col	 ����ĳһ�е���������1�е�������0��
����ֵ�������row+1�е�col+1�е�Ԫ��ֵ
�������ܣ����������()��ͨ������������ȡ����ĳһԪ�ص�ֵ���ҿɸı�
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
������������()
������row    ����ĳһ�е���������1�е�������0��
      col	 ����ĳһ�е���������1�е�������0��
����ֵ�������row+1�е�col+1�е�Ԫ��ֵ
�������ܣ����������()��ͨ������������ȡ����ĳһԪ�ص�ֵ�������ܸı�
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
������������+
������other    ��һ������
����ֵ������������ӵõ��ľ���
�������ܣ����������+��ʵ�־�������
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
������������-
������other    ��һ������(����)
����ֵ��������������õ��ľ���
�������ܣ����������-��ʵ�־�������
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
������������*
������other    �˺��ұߵľ���
����ֵ������������˵õ��ľ���
�������ܣ����������*��ʵ�־�������
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
������������*
������num    ����
����ֵ��һ���������һ������
�������ܣ����������*��ʵ�־�������
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
��������ת��
����ֵ��ת�Ⱥ�ľ���
�������ܣ�ʵ�־����ת��
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
��������������
������row1    ����ĳһ�е�����
      row2    ������һ�е�����
�������ܣ����������ĳ����
**********************************************************
*/
void Matrix::swapRows(size_t row1, size_t row2)
{
    swap_ranges(data[row1].begin(), data[row1].end(), data[row2].begin());
}

/*
**********************************************************
�����������г���һ������
������row     ����ĳһ�е�����
      scalar  ���˵ı���
�������ܣ��������ĳһ�г���һ������
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
��������һ�м�����һ�еı���
������fromrow    ����
      torow      ���ı����
      scalar     ����
�������ܣ���ĳһ�г���һ�������ӵ���һ���ϣ��ı���һ�е�ֵ
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
����������������
����ֵ������������
�������ܣ�ʹ�ø�˹��Ԫ������þ���������
**********************************************************
*/
Matrix Matrix::inverse() const
{
    if (rows() != cols())
    {
        throw invalid_argument("Matrix must be square for inversion");
    }
    Matrix A = *this;//��������ĸ���������ı�ԭ����
    size_t n = A.rows();
    if (n <= 0)
    {
        throw invalid_argument("Matrix's rows or cols must more than 0 ");
    }
    Matrix inverse(n, n, 0.0);

    // ������λ����
    for (size_t i = 0; i < n; ++i)
    {
        inverse(i, i) = 1.0;
    }

    // ��˹-Լ����Ԫ��
    for (size_t i = 0; i < n; ++i)
    {
        // Ѱ��������Ԫ��
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
        //���Խ�Ԫ�ع淶��Ϊ1�������������
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
��������ɾ����
������rowIndex    ����ĳһ�е�����
�������ܣ��������ĳһ��ɾ�����Ӷ����پ��������
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
��������ɾ����
������colIndex    ����ĳһ�е�����
�������ܣ��������ĳһ��ɾ�����Ӷ����پ��������
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
������������^
������other    ��һ����ά����
����ֵ��������ά������˵õ�����ά����
�������ܣ����������~��ʵ����ά�����Ĳ��
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
������������������ģ
����ֵ����ǰ������ģ
�������ܣ����㵱ǰ������ģ������
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
������������/
������other    ��һ������
����ֵ����ǰ��������һ������ƴ�Ӻ�ľ���
�������ܣ����������/��ʵ�־��������ƴ��
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
�����������أ�
������other    ��һ������
����ֵ����ǰ��������һ������ƴ�Ӻ�ľ���
�������ܣ������������������ʵ�־���ĺ���ƴ��
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
������������&
������other    ��һ����Ԫ��
����ֵ����ǰ��Ԫ������һ����Ԫ����˺�õ�����Ԫ��
�������ܣ����������&��ʵ����Ԫ���˷�
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
�����������ɵ�λ��
����ֵ����ǰ����γ�ȵĵ�λ��
�������ܣ����յ�ǰ�����γ������һ����λ��
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
�����������Գƾ���
����ֵ����ǰ��ά�����ķ��Գƾ���
�������ܣ�����ǰ����λ����ת��Ϊ���Գƾ���
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
���������Խ���
����ֵ���Ե�ǰ����Ϊ�Խ�Ԫ�صĶԽ���
�������ܣ��Ե�ǰ����Ϊ�Խ�Ԫ�أ����ɶԽ���
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
����������Ԫ����˺�λ��
��  ����q2        ��һ����Ԫ��
����ֵ����ǰ��Ԫ������һ����Ԫ����˲���λ����õ�����Ԫ��
�������ܣ�����ǰ��Ԫ������һ����Ԫ����˲���λ��
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