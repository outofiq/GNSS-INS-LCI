#include "KF_GINS.h"


/*
**********************************************************
函数名：滤波权阵更新
参  数：navstate     当前的导航状态
        thisimu      当前惯导数据
		dt           相邻两个惯导历元的时间差
		kf           当前滤波状态
		corrtime     相关时间
函数功能：完成卡尔曼滤波的观测噪声（权阵）更新
**********************************************************
*/
void InsPropagate(navstate_t navstate, IMUData_t thisimu, double dt, kf_t& kf, double corrtime)
{
	//适配代码准备
	Matrix pos(3, 1, 0.0);
	pos(0, 0) = navstate.pva.blh.B;
	pos(1, 0) = navstate.pva.blh.L;
	pos(2, 0) = navstate.pva.blh.H;
	double H = pos(2, 0);
	double sinB = sin(pos(0, 0));
	double cosB = cos(pos(0, 0));
	double tanB = tan(pos(0, 0));
	Matrix vel(3, 1, 0.0);
	vel(0, 0) = navstate.pva.Vn;
	vel(1, 0) = navstate.pva.Ve;
	vel(2, 0) = navstate.pva.Vd;
	double vn = vel(0, 0);
	double ve = vel(1, 0);
	double vd = vel(2, 0);
	Matrix Cbn = navstate.Cbn;
	Matrix Cnb = Cbn.transpose();
	double Rm = navstate.Rm;
	double Rn = navstate.Rn;
	double gravity = navstate.gravity;

	Matrix omega_ib(3, 1, 0.0);
	omega_ib(0, 0) = thisimu.gyro_x / dt;
	omega_ib(1, 0) = thisimu.gyro_y / dt;
	omega_ib(2, 0) = thisimu.gyro_z / dt;
	Matrix acce_ib(3, 1, 0.0);
	acce_ib(0, 0) = thisimu.accel_x / dt;
	acce_ib(1, 0) = thisimu.accel_y / dt;
	acce_ib(2, 0) = thisimu.accel_z / dt;

	Matrix wie_n(3, 1, 0.0);
	wie_n(0, 0) = Omega_WGS * cosB;
	wie_n(1, 0) = 0;
	wie_n(2, 0) = -Omega_WGS * sinB;
	Matrix wen_n(3, 1, 0.0);
	wen_n(0, 0) = ve / (Rn + H);
	wen_n(1, 0) = -vn / (Rm + H);
	wen_n(2, 0) = -ve * tanB / (Rn + H);
	Matrix win_n = wie_n + wen_n;

	Matrix F(kf.Rank, kf.Rank, 0.0);
	Matrix PHI(kf.Rank, kf.Rank, 0.0);
	PHI = PHI.eye();
	Matrix G(kf.Rank, kf.Noise_Rank, 0.0);

	// 以下代码基本为继承代码
	Matrix fvr(3, 3, 0.0);
	Matrix fvv(3, 3, 0.0);
	Matrix fphir(3, 3, 0.0);
	Matrix fphiv(3, 3, 0.0);

	Matrix frr(3, 3, 0.0);
	frr(0, 0) = -vd / (Rm + H);
	frr(0, 1) = 0;
	frr(0, 2) = vn / (Rm + H);
	frr(1, 0) = ve * tanB / (Rn + H);
	frr(1, 1) = -(vd + vn * tanB) / (Rn + H);
	frr(1, 2) = ve / (Rn + H);

	fvr(0, 0) = -2 * ve * Omega_WGS * cosB / (Rm + H) - ve * ve / ((Rn + H) * (Rm + H) * (cosB * cosB));
	fvr(0, 2) = vn * vd / ((Rm + H) * (Rm + H)) - ve * ve * tanB / ((Rn + H) * (Rn + H));
	fvr(1, 0) = 2 * Omega_WGS * (vn * cosB - vd * sinB) / (Rm + H) + vn * ve / ((Rn + H) * (Rm + H) * (cosB * cosB));
	fvr(1, 2) = ve * vd / ((Rn + H) * (Rn + H)) + vn * ve * tanB / ((Rn + H) * (Rn + H));
	fvr(2, 0) = 2 * ve * Omega_WGS * sinB / (Rm + H);
	fvr(2, 2) = -ve * ve / ((Rn + H) * (Rn + H)) - vn * vn / ((Rm + H) * (Rm + H)) + 2 * gravity / (sqrt(Rm * Rn) + H);

	fvv(0, 0) = vd / (Rm + H);
	fvv(0, 1) = -2 * (Omega_WGS * sinB + ve * tanB / (Rn + H));
	fvv(0, 2) = vn / (Rm + H);
	fvv(1, 0) = 2 * Omega_WGS * sinB + ve * tanB / (Rn + H);
	fvv(1, 1) = (vd + vn * tanB) / (Rn + H);
	fvv(1, 2) = 2 * Omega_WGS * cosB + ve / (Rn + H);
	fvv(2, 0) = -2 * vn / (Rm + H);
	fvv(2, 1) = -2 * (Omega_WGS * cosB + ve / (Rn + H));

	fphir(0, 0) = -Omega_WGS * sinB / (Rm + H);
	fphir(0, 2) = ve / ((Rn + H) * (Rn + H));
	fphir(1, 2) = -vn / ((Rm + H) * (Rm + H));
	fphir(2, 0) = -Omega_WGS * cosB / (Rm + H) - ve / ((Rn + H) * (Rm + H) * (cosB * cosB));
	fphir(2, 2) = -ve * tanB / ((Rn + H) * (Rn + H));

	fphiv(0, 1) = 1 / (Rn + H);
	fphiv(1, 0) = -1 / (Rm + H);
	fphiv(2, 1) = -tanB / (Rn + H);

	Matrix Cbn_acceib_skew = (Cbn * acce_ib).skew();
	Matrix win_n_skew = (win_n.skew()) * (-1);
	Matrix Cbn_ = Cbn * (-1);
	Matrix Cbn_omega_ib = Cbn_ * omega_ib.diag();
	Matrix Cbn_acce_ib = Cbn * (acce_ib.diag());

	for (int i = 0; i < 3; i++)
	{
		F(i, i + 3) = 1;
		F(i + 9, i + 9) = -1 / corrtime;
		F(i + 12, i + 12) = -1 / corrtime;
		F(i + 15, i + 15) = -1 / corrtime;
		F(i + 18, i + 18) = -1 / corrtime;
		for (int j = 0; j < 3; j++)
		{
			F(i, j) = frr(i, j);
			F(i + 3, j) = fvr(i, j);
			F(i + 3, j + 3) = fvv(i, j);
			F(i + 6, j) = fphir(i, j);
			F(i + 6, j + 3) = fphiv(i, j);

			F(i + 3, j + 6) = Cbn_acceib_skew(i, j);
			F(i + 6, j + 6) = win_n_skew(i, j);
			F(i + 6, j + 9) = Cbn_(i, j);
			F(i + 3, j + 12) = Cbn(i, j);
			F(i + 6, j + 15) = Cbn_omega_ib(i, j);
			F(i + 3, j + 18) = Cbn_acce_ib(i, j);
		}
	}

	//propagate covariance
	PHI = PHI + (F * dt);
	for (int i = 0; i < 3; i++)
	{
		G(i + 9, i + 6) = 1;
		G(i + 12, i + 9) = 1;
		G(i + 15, i + 12) = 1;
		G(i + 18, i + 15) = 1;
		for (int j = 0; j < 3; j++)
		{
			G(i + 3, j) = Cbn(i, j);
			G(i + 6, j + 3) = Cbn(i, j);
		}
	}
	Matrix Q = (G * kf.Qc * G.transpose()) * dt;
	Q = (PHI * Q * PHI.transpose() + Q) * 0.5;

	kf.P = PHI * kf.P * PHI.transpose() + Q;
}