#include <string>
#include "coordinate.h"

struct CFG_t
{
	string          staticIMU_File;    //静态惯导数据文件
	string          dynamicIMU_File;   //动态惯导数据文件
	string          GNSS_File;         //GNSS定位结果数据文件
	string          Result_File;       //结果输出文件

	const double    acc_scale;         //加速度计转换系数（单位：m/s^2/LSB）
	const double    gry_scale;         //陀螺仪转换系数（单位：°/s/LSB）
	const double    Tau0;              //采样间隔（单位：s）
	double          static_T;          //动态采集初始静止的时间
	double          movetime;          //小车开始移动的时间
	double          starttime;         
	double          endtime;

	Matrix          initpos;           //初始位置(rad,rad,m)
	Matrix          initvel;           //初始速度（m/s）
	Matrix          initatt;           //初始姿态（rad）
	Matrix          initposstd;        //初始位置std
	Matrix          initvelstd;        //初始速度std
	Matrix          initattstd;        //初始姿态std

	Matrix          initgrybias;       //初始陀螺零偏（rad/s）
	Matrix          initaccbias;       //初始加速度计零偏（m/s^2）
	Matrix          initgryscale;      //陀螺比例因子
	Matrix          initaccscale;      //加速度计比例因子
	Matrix          initgyrbiasstd;    //初始陀螺零偏std
	Matrix          initaccbiasstd;    //初始加速度计零偏std
	Matrix          initgyrscalestd;   //初始陀螺比例因子std
	Matrix          initaccscalestd;   //初始加速度计比例因子std
	Matrix          gyrbiasstd;        //初始陀螺零偏std
	Matrix          accbiasstd;        //初始加速度计零偏std
	Matrix          gyrscalestd;       //初始陀螺比例因子std
	Matrix          accscalestd;       //初始加速度计比例因子std

	Matrix          gryarw;            //角度随机游走
	Matrix          accvrw;            //速率随机游走

	Matrix          zupt_noise;        //零速观测噪声
	Matrix          zaru_noise;        //零角速率观测噪声

	Matrix          antlever;          //安装角（m）

	double          corrtime;          //相关时间
	CFG_t() :
		staticIMU_File("E:\\综合实习\\data\\static-17.ASC"),
		dynamicIMU_File("E:\\综合实习\\data\\dynamic-17.ASC"),
		GNSS_File("E:\\综合实习\\GNSSData\\GNSSData\\Export\\dynamic-17.pos"),
		Result_File("E:\\综合实习\\Result\\Result.txt"),

		acc_scale(3.74094009399414e-6),
		gry_scale(3.0517578125e-5),
		Tau0(0.008),
		static_T(300.0),
		movetime(0.0),
		starttime(0.0),
		endtime(0.0),

		initpos(3, 1, 0.0),
		initvel(3, 1, 0.0),
		initatt(3, 1, 0.0),
		initposstd(3, 1, 0.0),
		initvelstd(3, 1, 0.0),
		initattstd(3, 1, 0.0),

		initaccbias(3, 1, 0.0),
		initgrybias(3, 1, 0.0),
		initaccscale(3, 1, 0.0),
		initgryscale(3, 1, 0.0),
		initaccbiasstd(3, 1, 0.0),
		initgyrbiasstd(3, 1, 0.0),
		initaccscalestd(3, 1, 0.0),
		initgyrscalestd(3, 1, 0.0),
		accbiasstd(3, 1, 0.0),
		gyrbiasstd(3, 1, 0.0),
		accscalestd(3, 1, 0.0),
		gyrscalestd(3, 1, 0.0),

		gryarw(3, 1, 0.0),
		accvrw(3, 1, 0.0),

		zupt_noise(3, 1, 0.0),
		zaru_noise(3, 1, 0.0),

		antlever(3, 1, 0.0),

		corrtime(2 * 3600.0)
	{}
};