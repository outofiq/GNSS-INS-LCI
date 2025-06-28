#include "KF_GINS.h"


/*
**********************************************************
函数名：误差补偿
参  数：kf           当前的滤波状态
        navstate     当前的导航状态
函数功能：将当前滤波得到的参数残差补偿到导航状态中，
          更新位置、速度、姿态、零偏与比例因子，最后
		  将滤波状态的参数残差置零。
**********************************************************
*/
void ErrorFeedback(kf_t& kf, navstate_t& navstate)
{
	//position and velocity
	Matrix DR(3, 3, 0.0);
	DR(0, 0) = navstate.Rm + navstate.pva.blh.H;
	DR(1, 1) = (navstate.Rn + navstate.pva.blh.H) * cos(navstate.pva.blh.B);
	DR(2, 2) = -1;
	Matrix DR_inv = DR.inverse();

	Matrix Pos(3, 1, 0.0);
	Pos(0, 0) = navstate.pva.blh.B;
	Pos(1, 0) = navstate.pva.blh.L;
	Pos(2, 0) = navstate.pva.blh.H;
	Matrix kf_x_pos(3, 1, 0.0);
	for (int i = 0; i < 3; i++)
		kf_x_pos(i, 0) = kf.x(i, 0);
	Pos = Pos - (DR_inv * kf_x_pos);
	navstate.pva.blh.B = Pos(0, 0);
	navstate.pva.blh.L = Pos(1, 0);
	navstate.pva.blh.H = Pos(2, 0);

	Matrix Vel(3, 1, 0.0);
	Vel(0, 0) = navstate.pva.Vn;
	Vel(1, 0) = navstate.pva.Ve;
	Vel(2, 0) = navstate.pva.Vd;
	Matrix kf_x_vel(3, 1, 0.0);
	for (int i = 0; i < 3; i++)
		kf_x_vel(i, 0) = kf.x(i + 3, 0);
	Vel = Vel - kf_x_vel;
	navstate.pva.Vn = Vel(0, 0);
	navstate.pva.Ve = Vel(1, 0);
	navstate.pva.Vd = Vel(2, 0);

	//attitude
	Matrix kf_x_phi(3, 1, 0.0);
	for (int i = 0; i < 3; i++)
		kf_x_phi(i, 0) = kf.x(i + 6, 0);
	Matrix qpn = RotvecToQuat(kf_x_phi);
	navstate.qbn = qpn.quatProd(navstate.qbn);
	navstate.Cbn = QuatToCbn(navstate.qbn);
	navstate.pva.rph = CbnToEuler(navstate.Cbn);

	//imu error
	Matrix kf_x_gb(3, 1, 0.0);
	Matrix kf_x_ab(3, 1, 0.0);
	Matrix kf_x_gs(3, 1, 0.0);
	Matrix kf_x_as(3, 1, 0.0);
	for (int i = 0; i < 3; i++)
	{
		kf_x_gb(i, 0) = kf.x(i + 9, 0);
		kf_x_ab(i, 0) = kf.x(i + 12, 0);
		kf_x_gs(i, 0) = kf.x(i + 15, 0);
		kf_x_as(i, 0) = kf.x(i + 18, 0);
	}
	navstate.grybias = navstate.grybias + kf_x_gb;
	navstate.accbias = navstate.accbias + kf_x_ab;
	navstate.gryscale = navstate.gryscale + kf_x_gs;
	navstate.accscale = navstate.accscale + kf_x_as;

	//update some parameters
	navstate.Rm = getRm(navstate.pva.blh.B);
	navstate.Rn = getRn(navstate.pva.blh.B);
	navstate.gravity = getGravity(navstate.pva.blh);

	//reset state vector
	Matrix zero(kf.Rank, 1, 0.0);
	kf.x = zero;
}