#include "coordinate.h"

/*
**********************************************************
���������ѿ���ϵת�������ϵ
������Descartes    ����ĵѿ�������XYZ
	  Geodetic	   ת���õ��Ĵ������BLH
	  R            ���򳤰���
	  F            �������
�������ܣ����ѿ�������XYZת��Ϊ�������BLH
**********************************************************
*/
void XYZToBLH(const XYZ Descartes, BLH* Geodetic, const double R, const double F)
{
	double x = Descartes.x;
	double y = Descartes.y;
	double z = Descartes.z;

	const double e2 = 2.0 * F - F * F;

	double S = sqrt(x * x + y * y);
	double L = atan2(y, x);
	double B = 0;
	double N = 0;
	double tempB = atan2(z, S);

	int counter = 0;
	while (fabs(B - tempB) > 1e-12 && counter < 25)
	{
		B = tempB;
		N = R / sqrt(1 - e2 * sin(B) * sin(B));
		tempB = atan2(z + N * e2 * sin(B), S);
		counter++;
	}

	Geodetic->B = B;
	Geodetic->L = L;
	Geodetic->H = S / cos(B) - N;
}

/*
**********************************************************
���������������ϵת�ѿ�������ϵ
������Geodetic     ����Ĵ������BLH
	  Descartes	   ת���õ��ĵѿ�������XYZ
	  R            ���򳤰���
	  F            �������
�������ܣ��������������BLHת��Ϊ�ѿ�������XYZ
**********************************************************
*/
void BLHToXYZ(const BLH Geodetic, XYZ* Descartes, const double R, const double F)
{
	double B = Geodetic.B;
	double L = Geodetic.L;
	double H = Geodetic.H;
	double e2 = 2 * F - F * F;
	double N = R / sqrt(1 - e2 * sin(B) * sin(B));
	Descartes->x = (N + H) * cos(B) * cos(L);
	Descartes->y = (N + H) * cos(B) * sin(L);
	Descartes->z = (N * (1 - e2) + H) * sin(B);
}


/*
**********************************************************
����������ȡ��������
��  ����blh         ���صĴ������BLH
����ֵ����������ֵ
�������ܣ����뵱�صĴ�����꣬���㲢���ص��ص�����ֵ
**********************************************************
*/
double getGravity(const BLH blh)
{
	double g;
	double h = blh.H;
	double sinB = sin(blh.B);
	double sinB2 = sinB * sinB;
	double sinB4 = sinB2 * sinB2;
	double g0 = 9.7803267715 * (1 + 0.0052790414 * sinB2 + 0.0000232718 * sinB4);
	g = g0 - (3.087691089e-6 - 4.397731e-9 * sinB2) * h + 0.721e-12 * h * h;
	return g;
}

Euler CbnToEuler(const Matrix Cbn)
{
	if (Cbn.rows() != 3 || Cbn.cols() != 3)
	{
		throw invalid_argument("Matrices dimensions must match");
	}

	Euler attitude;
	double c11 = Cbn(0, 0);
	double c12 = Cbn(0, 1);
	double c13 = Cbn(0, 2);
	double c21 = Cbn(1, 0);
	double c22 = Cbn(1, 1);
	double c23 = Cbn(1, 2);
	double c31 = Cbn(2, 0);
	double c32 = Cbn(2, 1);
	double c33 = Cbn(2, 2);

	attitude.pitch = atan(-c31 / sqrt(c32 * c32 + c33 * c33));
	if (abs(c31) < 0.999)
	{
		attitude.roll = atan2(c32, c33);
		attitude.heading = atan2(c21, c11);
	}
	else if (c31 >= 0.999)
	{
		attitude.roll = NAN;
		attitude.heading = pi + atan2(c23 + c12, c13 - c22);
	}
	else if(c31<=-0.999)
	{
		attitude.roll = NAN;
		attitude.heading = atan2(c23 - c12, c13 + c22);
	}

	return attitude;
}

/*
**********************************************************
��������ŷ����ת��̬��ת����
��  ����attitude     �����ŷ����
����ֵ����̬��ת����
�������ܣ����ݵ�ǰ��ŷ����ת���õ���ǰ����̬��ת����Cbn
**********************************************************
*/
Matrix EulerToCbn(const Euler attitude)
{
	Matrix Cbn(3, 3, 0.0);
	double theta = attitude.pitch;
	double phi = attitude.roll;
	double psi = attitude.heading;

	double cos_theta = cos(theta);
	double sin_theta = sin(theta);
	double cos_phi = cos(phi);
	double sin_phi = sin(phi);
	double cos_psi = cos(psi);
	double sin_psi = sin(psi);

	Cbn(0, 0) = cos_theta * cos_psi;
	Cbn(0, 1) = -cos_phi * sin_psi + sin_phi * sin_theta * cos_psi;
	Cbn(0, 2) = sin_phi * sin_psi + cos_phi * sin_theta * cos_psi;
	Cbn(1, 0) = cos_theta * sin_psi;
	Cbn(1, 1) = cos_phi * cos_psi + sin_phi * sin_theta * sin_psi;
	Cbn(1, 2) = -sin_phi * cos_psi + cos_phi * sin_theta * sin_psi;
	Cbn(2, 0) = -sin_theta;
	Cbn(2, 1) = sin_phi * cos_theta;
	Cbn(2, 2) = cos_phi * cos_theta;

	return Cbn;
}

/*
**********************************************************
��������ŷ����ת��Ԫ��
��  ����attitude     �����ŷ����
����ֵ����Ԫ��
�������ܣ����ݵ�ǰ��ŷ����ת���õ���ǰ����Ԫ��
**********************************************************
*/
Matrix EulerToQuat(const Euler attitude)
{
	Matrix q(4, 1, 0.0);

	double roll = attitude.roll / 2;
	double pitch = attitude.pitch / 2;
	double heading = attitude.heading / 2;

	double cosRoll = cos(roll);
	double sinRoll = sin(roll);
	double cosPitch = cos(pitch);
	double sinPitch = sin(pitch);
	double cosHeading = cos(heading);
	double sinHeading = sin(heading);

	q(0, 0) = cosRoll * cosPitch * cosHeading + sinRoll * sinPitch * sinHeading;
	q(1, 0) = sinRoll * cosPitch * cosHeading - cosRoll * sinPitch * sinHeading;
	q(2, 0) = cosRoll * sinPitch * cosHeading + sinRoll * cosPitch * sinHeading;
	q(3, 0) = cosRoll * cosPitch * sinHeading - sinRoll * sinPitch * cosHeading;

	return q;
}

/*
**********************************************************
����������Ԫ��ת��̬��ת����
��  ����q      �������Ԫ��
����ֵ����̬��ת����
�������ܣ����ݵ�ǰ����Ԫ��ת���õ���ǰ����̬��ת����Cbn
**********************************************************
*/
Matrix QuatToCbn(const Matrix q)
{
	if (q.rows() != 4 || q.cols() != 1)
	{
		throw invalid_argument("Matrices dimensions must match");
	}

	Matrix Cbn(3, 3, 0.0);

	double q1 = q(0, 0), q2 = q(1, 0), q3 = q(2, 0), q4 = q(3, 0);
	double q1_2 = q1 * q1;
	double q2_2 = q2 * q2;
	double q3_2 = q3 * q3;
	double q4_2 = q4 * q4;

	Cbn(0, 0) = q1_2 + q2_2 - q3_2 - q4_2;
	Cbn(0, 1) = 2 * (q2 * q3 - q1 * q4);
	Cbn(0, 2) = 2 * (q2 * q4 + q1 * q3);
	Cbn(1, 0) = 2 * (q2 * q3 + q1 * q4);
	Cbn(1, 1) = q1_2 - q2_2 + q3_2 - q4_2;
	Cbn(1, 2) = 2 * (q3 * q4 - q1 * q2);
	Cbn(2, 0) = 2 * (q2 * q4 - q1 * q3);
	Cbn(2, 1) = 2 * (q3 * q4 + q1 * q2);
	Cbn(2, 2) = q1_2 - q2_2 - q3_2 + q4_2;

	return Cbn;
}

/*
**********************************************************
����������Ч��תʸ��ת��Ԫ��
��  ����Rot     ����ĵ�Ч��תʸ��
����ֵ����Ԫ��
�������ܣ����ݵ�ǰ�ĵ�Ч��תʸ��ת���õ���Ԫ��
**********************************************************
*/
Matrix RotvecToQuat(const Matrix Rot)
{
	if (Rot.rows() != 3 || Rot.cols() != 1)
	{
		throw invalid_argument("Matrices dimensions must match");
	}

	Matrix q(4, 1, 0.0);

	double angel = Rot.norm();
	double index;
	if (angel < 0.000005)
		index = 0.5;
	else
		index = sin(0.5 * angel) / angel;

	q(0, 0) = cos(0.5 * angel);
	q(1, 0) = index * Rot(0, 0);
	q(2, 0) = index * Rot(1, 0);
	q(3, 0) = index * Rot(2, 0);

	/*
	if (q(0, 0) < 0)
		q = q * (-1);
	*/

	return q;
}

/*
**********************************************************
����������Ԫ��ת��Ч��תʸ��
��  ����q       �������Ԫ��
����ֵ����Ч��תʸ��
�������ܣ����ݵ�ǰ����Ԫ��ת���õ���Ч��תʸ��
**********************************************************
*/
Matrix QuatToRotvec(const Matrix q)
{
	if (q.rows() != 4 || q.cols() != 1)
	{
		throw invalid_argument("Matrices dimensions must match");
	}

	Matrix rotvec(3, 1, 0.0);

	double w = q(0, 0);
	double x = q(1, 0);
	double y = q(2, 0);
	double z = q(3, 0);

	double qv_norm = sqrt(x * x + y * y + z * z);

	if (w == 0)
	{
		rotvec(0, 0) = x * pi;
		rotvec(1, 0) = y * pi;
		rotvec(2, 0) = z * pi;
	}
	else if(qv_norm < 1e-10)
	{
		rotvec = rotvec;
	}
	else
	{
		double temp = 2 * atan(qv_norm / w);
		double f = sin(temp / 2) / temp;
		rotvec(0, 0) = x / f;
		rotvec(1, 0) = y / f;
		rotvec(2, 0) = z / f;
	}

	return rotvec;
}

/*
**********************************************************
����������ȡ��������Ȧ���ʰ뾶Rm
��  ����phi     ���ص�γ��
����ֵ������Ȧ���ʰ뾶Rm
�������ܣ����ݵ��ص�γ�ȼ��㲢���ص��ص�����Ȧ���ʰ뾶Rm
**********************************************************
*/
double getRm(const double phi)
{
	double sinphi = sin(phi);
	double Rm = R_WGS84 * (1 - e * e) / sqrt(pow((1 - e * e * sinphi * sinphi), 3));
	return Rm;
}

/*
**********************************************************
����������ȡ����î��Ȧ���ʰ뾶Rn
��  ����phi     ���ص�γ��
����ֵ��î��Ȧ���ʰ뾶Rm
�������ܣ����ݵ��ص�γ�ȼ��㲢���ص��ص�î��Ȧ���ʰ뾶Rn
**********************************************************
*/
double getRn(const double phi)
{
	double sinphi = sin(phi);
	double Rn = R_WGS84 / sqrt(1 - e * e * sinphi * sinphi);
	return Rn;
}


/*
**********************************************************
����������γ��תqne
��  ����b          ���ص�γ��
	    l	       ���صľ���
		qne        ת���õ���qne����
�������ܣ�ͨ����γ��ת���õ�qne����
**********************************************************
*/
void blToqne(const double b, const double l, Matrix& qne)
{
	double s1 = sin(l / 2);
	double c1 = cos(l / 2);
	double s2= sin(-pi / 4 - b / 2);
	double c2 = cos(-pi / 4 - b / 2);

	qne(0, 0) = c1 * c2;
	qne(1, 0) = -s1 * s2;
	qne(2, 0) = c1 * s2;
	qne(3, 0) = c2 * s1;
}


/*
**********************************************************
��������qneת��γ��
��  ����qne        qne����
        b          ת���õ��ĵ���γ��
		l	       ת���õ��ĵ��ؾ���
�������ܣ�ͨ��qne����ת���õ����ؾ�γ��
**********************************************************
*/
void qneTobl(const Matrix qne, double& b, double& l)
{
	b= -2 * atan(qne(2,0) / qne(0,0)) - pi / 2;
	l = 2 * atan2(qne(3, 0), qne(0, 0));
}
