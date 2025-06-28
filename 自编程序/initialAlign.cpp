#include "initialAlign.h"

/*
**********************************************************
函数名：初始对准
参  数：imudata      惯导传感器数据
		cfg          配置参数结构体（包含初始姿态）
函数功能：利用初始静态惯导数据完成静态粗对准
**********************************************************
*/
void initialAlign(vector<IMUData_t> imuData, CFG_t& cfg)
{

	int i = 0;
	double gx_raw, gy_raw, gz_raw;
	gx_raw = gy_raw = gz_raw = 0.0;
	double wx_raw, wy_raw, wz_raw;
	wx_raw = wy_raw = wz_raw = 0.0;
	while (imuData[i].time < cfg.starttime)
	{
		gx_raw += imuData[i].accel_x;
		gy_raw += imuData[i].accel_y;
		gz_raw += imuData[i].accel_z;
		wx_raw += imuData[i].gyro_x;
		wy_raw += imuData[i].gyro_y;
		wz_raw += imuData[i].gyro_z;
		i += 1;
	}
	gx_raw /= i;
	gy_raw /= i;
	gz_raw /= i;
	wx_raw /= i;
	wy_raw /= i;
	wz_raw /= i;

	//b系
	Matrix gb(3, 1, 0.0);
	gb(0, 0) = -gx_raw; gb(1, 0) = -gy_raw; gb(2, 0) = -gz_raw;
	Matrix wb(3, 1, 0.0);
	wb(0, 0) = wx_raw; wb(1, 0) = wy_raw; wb(2, 0) = wz_raw;

	//n系（已单位化）
	Matrix Vg(3, 1, 0.0);
	Vg(0, 0) = 0.0; Vg(1, 0) = 0.0; Vg(2, 0) = 1.0;
	Matrix Vw(3, 1, 0.0);
	Vw(0, 0) = 0.0; Vw(1, 0) = 1.0; Vw(2, 0) = 0.0;
	Matrix Vgw(3, 1, 0.0);
	Vgw(0, 0) = 1.0; Vgw(1, 0) = 0.0; Vgw(2, 0) = 0.0;

	//单位化
	Matrix Wg = gb * (1 / gb.norm());
	Matrix Ww = (gb ^ wb) * (1 / (gb ^ wb).norm());
	Matrix Wgw = (gb ^ wb ^ gb) * (1 / (gb ^ wb ^ gb).norm());

	Matrix V = Vg.transpose() / Vw.transpose() / Vgw.transpose();
	Matrix W = Wg.transpose() / Ww.transpose() / Wgw.transpose();
	Matrix Cbn = V.inverse() * W;

	Euler initAtt = CbnToEuler(Cbn);

	cfg.initatt(0, 0) = initAtt.roll;
	cfg.initatt(1, 0) = initAtt.pitch;
	cfg.initatt(2, 0) = initAtt.heading;
}