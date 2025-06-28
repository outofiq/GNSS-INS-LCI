#include "KF_GINS.h"


/*
**********************************************************
函数名：初始化
参数：kf             当前滤波状态
	  navstate       当前导航状态
	  cfg            配置参数结构体
函数功能：定义滤波的初始过程噪声Qc、权阵P和初始导航状态
**********************************************************
*/
void Initialize(kf_t& kf, navstate_t& navstate, CFG_t cfg)
{
	//Qc过程噪声
	for (int i = 0; i < 3; i++)
	{
		kf.Qc(i, i) = cfg.accvrw(i, 0) * cfg.accvrw(i, 0);
		kf.Qc(i + 3, i + 3) = cfg.gryarw(i, 0) * cfg.gryarw(i, 0);
		kf.Qc(i + 6, i + 6) = 2 * cfg.gyrbiasstd(i, 0) * cfg.gyrbiasstd(i, 0) / cfg.corrtime;
		kf.Qc(i + 9, i + 9) = 2 * cfg.accbiasstd(i, 0) * cfg.accbiasstd(i, 0) / cfg.corrtime;
		kf.Qc(i + 12, i + 12) = 2 * cfg.gyrscalestd(i, 0) * cfg.gyrscalestd(i, 0) / cfg.corrtime;
		kf.Qc(i + 15, i + 15) = 2 * cfg.accscalestd(i, 0) * cfg.accscalestd(i, 0) / cfg.corrtime;
	}

	//P0
	for (int i = 0; i < 3; i++)
	{
		kf.P(i, i) = cfg.initposstd(i, 0) * cfg.initposstd(i, 0);
		kf.P(i+3, i+3) = cfg.initvelstd(i, 0) * cfg.initvelstd(i, 0);
		kf.P(i+6, i+6) = cfg.initattstd(i, 0) * cfg.initattstd(i, 0);
		kf.P(i + 9, i + 9) = cfg.initgyrbiasstd(i, 0) * cfg.initgyrbiasstd(i, 0);
		kf.P(i + 12, i + 12) = cfg.initaccbiasstd(i, 0) * cfg.initaccbiasstd(i, 0);
		kf.P(i + 15, i + 15) = cfg.initgyrscalestd(i, 0) * cfg.initgyrscalestd(i, 0);
		kf.P(i + 18, i + 18) = cfg.initaccscalestd(i, 0) * cfg.initaccscalestd(i, 0);
	}

	//navigation state initialization
	navstate.time = cfg.starttime;
	BLH blh;
	blh.B = cfg.initpos(0, 0);
	blh.L = cfg.initpos(1, 0);
	blh.H = cfg.initpos(2, 0);
	navstate.pva.blh = blh;
	navstate.pva.Vn = cfg.initvel(0, 0);
	navstate.pva.Ve = cfg.initvel(1, 0);
	navstate.pva.Vd = cfg.initvel(2, 0);
	Euler rph;
	rph.roll = cfg.initatt(0, 0);
	rph.pitch = cfg.initatt(1, 0);
	rph.heading = cfg.initatt(2, 0);
	navstate.pva.rph = rph;
	navstate.Cbn = EulerToCbn(rph);
	navstate.qbn = EulerToQuat(rph);
	navstate.grybias = cfg.initgrybias;
	navstate.accbias = cfg.initaccbias;
	navstate.gryscale = cfg.initgryscale;
	navstate.accscale = cfg.initaccscale;
	navstate.Rm = getRm(blh.B);
	navstate.Rn = getRn(blh.B);
	navstate.gravity = getGravity(blh);
}