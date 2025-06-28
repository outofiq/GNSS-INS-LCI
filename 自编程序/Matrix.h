#include <iostream>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
using namespace std;

namespace Matrix_hx
{
	class Matrix
	{
	private:
		vector<vector<double>> data;
	public:
		//构造函数
		Matrix(size_t rows, size_t cols, double value);
		//获取矩阵的行数
		size_t rows() const;
		//获取矩阵的列数
		size_t cols() const;
		//打印矩阵
		void printMatrix() const;
		//重载()运算符，用于访问矩阵元素
		double& operator()(size_t row, size_t col);
		//重载()运算符，用于访问矩阵元素(const版本)
		const double& operator()(size_t row, size_t col) const;
		//重载+运算符，实现矩阵加法
		Matrix operator+(const Matrix& other) const;
		//重载-运算符，实现矩阵减法
		Matrix operator-(const Matrix& other) const;
		//重载*运算符，实现矩阵乘法
		Matrix operator*(const Matrix& other) const;
		Matrix operator*(const double num) const;
		//实现矩阵的转置
		Matrix transpose() const;
		//矩阵求逆(采用高斯消元法，包括交换两行、将一行乘以一个标量、将一行加上另一行的倍数)
		//交换两行
		void swapRows(size_t row1, size_t row2);
		//将一行乘以一个标量
		void multiplyRow(size_t row, double scalar);
		//将一行加上另一行的倍数
		void addRowToRow(size_t fromRow, size_t toRow, double scalar);
		//矩阵求逆
		Matrix inverse() const;
		//删除行
		void deletRow(int rowIndex);
		//删除列
		void deleteColumn(int colIndex);
		//重载^运算符，实现三维向量叉乘
		Matrix operator^(const Matrix& other) const;
		//向量的模
		double norm() const;
		//重载,运算符，实现矩阵横向拼接
		Matrix operator,(const Matrix& other) const;
		//重载/运算符，实现矩阵纵向拼接
		Matrix operator/(const Matrix& other) const;
		//重载$运算符，实现四元数乘法
		Matrix operator&(const Matrix& other) const;
		//获得单位阵
		Matrix eye() const;
		//将三维向量转化为反对称矩阵
		Matrix skew() const;
		//将向量转化为对角阵
		Matrix diag() const;
		//进行四元数乘法后单位化
		Matrix quatProd(const Matrix q2) const;
	};
}
