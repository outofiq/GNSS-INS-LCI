#include "KF_GINS.h"


/*
**********************************************************
函数名：惯导机械编排
参数：lastimu         上一时刻的惯导数据
	  thisimu         本时刻的惯导数据
	  last2state      上上时刻的导航状态
	  laststate       上一时刻的导航状态
	  navstate        本时刻的导航状态
函数功能：根据惯导机械编排算法，完成纯惯导的导航状态更新
**********************************************************
*/
void INSMesh(const IMUData_t lastimu, const IMUData_t thisimu, const navstate_t last2state, const navstate_t laststate, navstate_t& navstate)
{
	//copy last navstate
	double lat_last = laststate.pva.blh.B;
	double lon_last = laststate.pva.blh.L;
	double h_last = laststate.pva.blh.H;
	BLH blh_last;
	blh_last.B = lat_last;
	blh_last.L = lon_last;
	blh_last.H = h_last;
	double vn_last = laststate.pva.Vn;
	double ve_last = laststate.pva.Ve;
	double vd_last = laststate.pva.Vd;
	Matrix lastqbn = laststate.qbn;
	Matrix lastcbn = laststate.Cbn;
	double gravity = laststate.gravity;

	//copy imu data
	Matrix last_dtheta(3, 1, 0.0);
	last_dtheta(0, 0) = lastimu.gyro_x;
	last_dtheta(1, 0) = lastimu.gyro_y;
	last_dtheta(2, 0) = lastimu.gyro_z;
	Matrix last_dvel(3, 1, 0.0);
	last_dvel(0, 0) = lastimu.accel_x;
	last_dvel(1, 0) = lastimu.accel_y;
	last_dvel(2, 0) = lastimu.accel_z;
	Matrix this_dtheta(3, 1, 0.0);
	this_dtheta(0, 0) = thisimu.gyro_x;
	this_dtheta(1, 0) = thisimu.gyro_y;
	this_dtheta(2, 0) = thisimu.gyro_z;
	Matrix this_dvel(3, 1, 0.0);
	this_dvel(0, 0) = thisimu.accel_x;
	this_dvel(1, 0) = thisimu.accel_y;
	this_dvel(2, 0) = thisimu.accel_z;

	//等效旋转矢量更新b系
	Matrix phi = this_dtheta + ((last_dtheta ^ this_dtheta) * (1.0 / 12.0));
	double phi_mo = phi.norm();
	Matrix qb(4, 1, 0.0);
	qb(0, 0) = cos(0.5 * phi_mo);
	Matrix qb_3 = phi * 0.5 * (sin(0.5 * phi_mo) / (0.5 * phi_mo));
	for (int i = 0; i < 3; i++)
		qb(i + 1, 0) = qb_3(i, 0);

	//等效旋转矢量更新n系
	Matrix wie(3, 1, 0.0);
	wie(0, 0) = Omega_WGS * cos(lat_last);
	wie(1, 0) = 0.0;
	wie(2, 0) = -Omega_WGS * sin(lat_last);
	double Rm = getRm(lat_last);
	double Rn = getRn(lat_last);
	Matrix wen(3, 1, 0.0);
	wen(0, 0) = ve_last / (Rn + h_last);
	wen(1, 0) = -vn_last / (Rm + h_last);
	wen(2, 0) = -ve_last * tan(lat_last) / (Rn + h_last);
	double delta_t = thisimu.time - lastimu.time;
	Matrix zeta = (wie + wen) * delta_t;
	double zeta_mo = zeta.norm();
	Matrix qn(4, 1, 0.0);
	qn(0, 0) = cos(0.5 * zeta_mo);
	Matrix qn_3 = zeta * (-0.5) * (sin(0.5 * zeta_mo) / (0.5 * zeta_mo));
	for (int i = 0; i < 3; i++)
		qn(i + 1, 0) = qn_3(i, 0);

	//计算当前姿态四元数
	Matrix qbn = qn & lastqbn & qb;
	double qbn_mo = qbn.norm();
	for (int i = 0; i < 4; i++)
		qbn(i, 0) /= qbn_mo;
	//四元数转方向余弦矩阵
	Matrix Cbn = QuatToCbn(qbn);
	navstate.qbn = qbn;
	navstate.Cbn = Cbn;
	navstate.pva.rph = CbnToEuler(Cbn);

	//速度更新
	Matrix delta_vb = this_dvel + ((this_dtheta ^ this_dvel) * (1.0 / 2.0)) + (((last_dtheta ^ this_dvel) + (last_dvel ^ this_dtheta)) * (1.0 / 12.0));
	//线性外推
	Matrix wie_last(3, 1, 0.0);
	wie_last(0, 0) = Omega_WGS * cos(last2state.pva.blh.B);
	wie_last(1, 0) = 0.0;
	wie_last(2, 0) = -Omega_WGS * sin(last2state.pva.blh.B);
	Matrix wie_mean = (wie * (3.0 / 2.0)) - (wie_last * (1.0 / 2.0));
	double Rm_last = getRm(last2state.pva.blh.B);
	double Rn_last = getRn(last2state.pva.blh.B);
	Matrix wen_last(3, 1, 0.0);
	wen_last(0, 0) = last2state.pva.Ve / (Rn_last + last2state.pva.blh.H);
	wen_last(1, 0) = -last2state.pva.Vn / (Rm_last + last2state.pva.blh.H);
	wen_last(2, 0) = -last2state.pva.Ve * tan(last2state.pva.blh.B) / (Rn_last + last2state.pva.blh.H);
	Matrix wen_mean = (wen * (3.0 / 2.0)) - (wen_last * (1.0 / 2.0));
	Matrix zeta_mean = (wie_mean + wen_mean) * delta_t;

	Matrix I3(3, 3, 0.0);
	I3 = I3.eye();
	Matrix delta_vf = (I3 - (zeta_mean.skew() * 0.5)) * lastcbn * delta_vb;
	double g = getGravity(blh_last);
	BLH blh_last2;
	blh_last2.B = last2state.pva.blh.B;
	blh_last2.L = last2state.pva.blh.L;
	blh_last2.H = last2state.pva.blh.H;
	double g_last = getGravity(blh_last2);
	double g_mean = (g * (3.0 / 2.0)) - (g_last * (1.0 / 2.0));
	Matrix gpn(3, 1, 0.0);
	gpn(2, 0) = g_mean;
	Matrix v_last(3, 1, 0.0);
	v_last(0, 0) = vn_last;
	v_last(1, 0) = ve_last;
	v_last(2, 0) = vd_last;
	Matrix v_last2(3, 1, 0.0);
	v_last2(0, 0) = last2state.pva.Vn;
	v_last2(1, 0) = last2state.pva.Ve;
	v_last2(2, 0) = last2state.pva.Vd;
	Matrix v_mean = (v_last * (3.0 / 2.0)) - (v_last2 * (1.0 / 2.0));
	Matrix delta_vg = (gpn - ((wie_mean * 2 + wen_mean) ^ v_mean)) * delta_t;
	Matrix v = v_last + delta_vf + delta_vg;
	navstate.pva.Vn = v(0, 0);
	navstate.pva.Ve = v(1, 0);
	navstate.pva.Vd = v(2, 0);

	//位置更新
	double h = h_last - (1.0 / 2.0) * (vd_last + navstate.pva.Vd) * delta_t;
	double h_mean = (h_last + h) / 2;
	double lat = lat_last + (vn_last + navstate.pva.Vn) * delta_t / (2 * (Rm + h_mean));
	double lat_mean = (lat_last + lat) / 2;
	double Rn_mean = getRn(lat_mean);
	double lon = lon_last + (ve_last + navstate.pva.Ve) * delta_t / (2 * (Rn_mean + h_mean) * cos(lat_mean));
	navstate.pva.blh.B = lat;
	navstate.pva.blh.L = lon;
	navstate.pva.blh.H = h;
	navstate.time = thisimu.time;
}