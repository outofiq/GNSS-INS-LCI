#include "Matrix.h"
using namespace Matrix_hx;

// �������
#define R_WGS84  6378137.0                // WGS84����ĳ����� m 
#define b_WGS84  6356752.3142             // WGS84����Ķ̰��� m
#define F_WGS84  1.0/298.257223563        // WGS84����ı���
#define Omega_WGS 7.2921151467e-5         // WGS84�������ת���ٶ� rad/s
#define GM_WGS   398600.5e+9              // WGS84���������������GM m3/s2
#define e        0.08181919104            // ��������ģ�͵�һƫ����

#define R_CGS2K  6378137.0                // CGCS2000����ĳ����� m
#define F_CGS2K  1.0/298.257222101        // CGCS2000����ı���
#define Omega_BDS 7.2921150e-5            // CGCS2000�������ת���ٶ� rad/s
#define GM_BDS   398600.4418e+9           // CGCS2000���������������GM m3/s2

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

//XYZת��ΪBLH
void XYZToBLH(const XYZ Descartes, BLH* Geodetic, const double R, const double F);
//BLHת��ΪXYZ
void BLHToXYZ(const BLH Geodetic, XYZ* Descartes, const double R, const double F);
//�����������ٶ�g
double getGravity(const BLH blh);
//��̬����Cbnתŷ����
Euler CbnToEuler(const Matrix Cbn);
//ŷ����ת��̬����Cbn
Matrix EulerToCbn(const Euler attitude);
//ŷ����ת��Ԫ��q
Matrix EulerToQuat(const Euler attitude);
//��Ԫ��ת�������Ҿ���
Matrix QuatToCbn(const Matrix q);
//��ת����ת��Ԫ��
Matrix RotvecToQuat(const Matrix Rot);
//��Ԫ��ת��ת����
Matrix QuatToRotvec(const Matrix q);
//����î��Ȧ���ʰ뾶
double getRm(const double phi);
//��������Ȧ���ʰ뾶
double getRn(const double phi);
//BLתqne
void blToqne(const double b, const double l, Matrix& qne);
//qneתBL
void qneTobl(const Matrix qne, double& b, double& l);