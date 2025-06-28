#include "initialAlign.h"

struct PVA_t
{
	double time;
	BLH blh;
	Euler rph;
	// NED
	double Vn;
	double Ve;
	double Vd;

	PVA_t()
	{
		time = 0.0;
		blh.B = blh.L = blh.H = 0.0;
		rph.roll = rph.pitch = rph.heading = 0.0;
		Vn = Ve = Vd = 0.0;
	}
};

struct kf_t
{
	int Rank;
	int Noise_Rank;
	Matrix P;
	Matrix Qc;
	Matrix x;
	kf_t() :
		Rank(21),
		Noise_Rank(18),
		P(Rank, Rank, 0.0),
		Qc(Noise_Rank, Noise_Rank, 0.0),
		x(Rank, 1, 0.0)
	{}
};

struct navstate_t
{
	double time;
	PVA_t  pva;
	Matrix Cbn;
	Matrix qbn;
	Matrix grybias;
	Matrix accbias;
	Matrix gryscale;
	Matrix accscale;
	double Rm;
	double Rn;
	double gravity;
	navstate_t() :
		Cbn(3, 3, 0.0),
		qbn(4, 1, 0.0),
		grybias(3, 1, 0.0),
		accbias(3, 1, 0.0),
		gryscale(3, 1, 0.0),
		accscale(3, 1, 0.0)
	{}
};

struct GNSSData_t
{
	double time;
	BLH blh;
	double Vel[3];
	double PosStd[3];
	double VelStd[3];
};

//滤波器和状态初始化
void Initialize(kf_t& kf, navstate_t& navstate, CFG_t cfg);
//读取GNSS数据
vector<GNSSData_t> readGNSSData(const string gnssfile);
//GNSS位置更新
void GNSSUpdate_POS(const navstate_t navstate, const GNSSData_t gnssdata, kf_t& kf, const Matrix antlever);
//GNSS速度更新
void GNSSUpdata_Vel(const navstate_t navstate, const GNSSData_t gnssdata, kf_t& kf, const Matrix antlever, const IMUData_t thisimu, const double dt);
//误差补偿
void ErrorFeedback(kf_t& kf, navstate_t& navstate);
//惯导机械编排
void INSMesh(const IMUData_t lastimu, const IMUData_t thisimu, const navstate_t last2state, const navstate_t laststate, navstate_t& navstate);
//卡尔曼滤波状态更新
void InsPropagate(const navstate_t navstate, const IMUData_t thisimu, const double dt, kf_t& kf, const double corrtime);
//内插imu数据
void interpolate(const IMUData_t lastimu, const IMUData_t thisimu, const double intertime, IMUData_t& firstimu, IMUData_t& secondimu);