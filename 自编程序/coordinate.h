#include "Matrix.h"
using namespace Matrix_hx;

// 地球参数
#define R_WGS84  6378137.0                // WGS84椭球的长半轴 m 
#define b_WGS84  6356752.3142             // WGS84椭球的短半轴 m
#define F_WGS84  1.0/298.257223563        // WGS84椭球的扁率
#define Omega_WGS 7.2921151467e-5         // WGS84椭球的自转角速度 rad/s
#define GM_WGS   398600.5e+9              // WGS84椭球地心引力常数GM m3/s2
#define e        0.08181919104            // 地球椭球模型第一偏心率

#define R_CGS2K  6378137.0                // CGCS2000椭球的长半轴 m
#define F_CGS2K  1.0/298.257222101        // CGCS2000椭球的扁率
#define Omega_BDS 7.2921150e-5            // CGCS2000椭球的自转角速度 rad/s
#define GM_BDS   398600.4418e+9           // CGCS2000椭球地心引力常数GM m3/s2

#define pi       3.14159265359

struct XYZ
{
	double x;
	double y;
	double z;
	XYZ()
	{
		x = 0.0;
		y = 0.0;
		z = 0.0;
	}
};

struct BLH
{
	double B;
	double L;
	double H;
	BLH()
	{
		B = 0.0;
		L = 0.0;
		H = 0.0;
	}
};

struct Euler
{
	//roll pitch heading
	double roll;
	double pitch;
	double heading;
};

//XYZ转化为BLH
void XYZToBLH(const XYZ Descartes, BLH* Geodetic, const double R, const double F);
//BLH转换为XYZ
void BLHToXYZ(const BLH Geodetic, XYZ* Descartes, const double R, const double F);
//计算重力加速度g
double getGravity(const BLH blh);
//姿态矩阵Cbn转欧拉角
Euler CbnToEuler(const Matrix Cbn);
//欧拉角转姿态矩阵Cbn
Matrix EulerToCbn(const Euler attitude);
//欧拉角转四元数q
Matrix EulerToQuat(const Euler attitude);
//四元数转方向余弦矩阵
Matrix QuatToCbn(const Matrix q);
//旋转向量转四元数
Matrix RotvecToQuat(const Matrix Rot);
//四元数转旋转向量
Matrix QuatToRotvec(const Matrix q);
//计算卯酉圈曲率半径
double getRm(const double phi);
//计算子午圈曲率半径
double getRn(const double phi);
//BL转qne
void blToqne(const double b, const double l, Matrix& qne);
//qne转BL
void qneTobl(const Matrix qne, double& b, double& l);