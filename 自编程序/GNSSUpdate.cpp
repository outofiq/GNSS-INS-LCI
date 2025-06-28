#include "KF_GINS.h"


/*
**********************************************************
函数名：读取GNSS数据
参  数：gnssfile        GNSS结果数据文件名
返回值：GNSS结果数据序列
函数功能：从GNSS结果数据文件中读取GNSS位置、速度及二者的std信息
**********************************************************
*/
vector<GNSSData_t> readGNSSData(const string gnssfile)
{
	vector<GNSSData_t> gnssdata;
	ifstream file(gnssfile);
	if (!file.is_open())
	{
		cerr << "open error!";
	}

	string line;
	//跳过前两行（变量说明）
	getline(file, line);
	getline(file, line);

	double time, B, L, H, Vel[3], PosStd[3], VelStd[3];
	while (getline(file, line))
	{
		stringstream ss(line);
		try
		{
			GNSSData_t entry;
			ss >> time
				>> B
				>> L
				>> H
				>> Vel[0]
				>> Vel[1]
				>> Vel[2]
				>> PosStd[0]
				>> PosStd[1]
				>> PosStd[2]
				>> VelStd[0]
				>> VelStd[1]
				>> VelStd[2];

			entry.time = time;
			entry.blh.B = B * DegToRad;
			entry.blh.L = L * DegToRad;
			entry.blh.H = H;
			entry.Vel[0] = Vel[0];
			entry.Vel[1] = Vel[1];
			entry.Vel[2] = -Vel[2];
			entry.PosStd[0] = PosStd[0];
			entry.PosStd[1] = PosStd[1];
			entry.PosStd[2] = PosStd[2];
			entry.VelStd[0] = VelStd[0];
			entry.VelStd[1] = VelStd[1];
			entry.VelStd[2] = VelStd[2];
			
			gnssdata.push_back(entry);
		}
		catch (...)
		{
			std::cerr << "Error parsing line: " << line << std::endl;
		}
	}
	file.close();

	return gnssdata;
}

/*
**********************************************************
函数名：惯导数据差值
参  数：lastimu          上一时刻的惯导数据
		thisimu          当前时刻的惯导数据
		intertime	     内插的时间
		firstimu         内插得到的第一个惯导数据
		secondimu        内插得到的第二个惯导数据
函数功能：在上一时刻和当前时刻中内插得到两组惯导数据
**********************************************************
*/
void interpolate(IMUData_t lastimu, IMUData_t thisimu, double intertime, IMUData_t& firstimu, IMUData_t& secondimu)
{
	if (intertime >= lastimu.time && intertime <= thisimu.time)
	{
		double lambda = (intertime - lastimu.time) / (thisimu.time - lastimu.time);

		firstimu.time = intertime;
		firstimu.gyro_x = thisimu.gyro_x * lambda;
		firstimu.gyro_y = thisimu.gyro_y * lambda;
		firstimu.gyro_z = thisimu.gyro_z * lambda;
		firstimu.accel_x = thisimu.accel_x * lambda;
		firstimu.accel_y = thisimu.accel_y * lambda;
		firstimu.accel_z = thisimu.accel_z * lambda;

		secondimu.time = thisimu.time;
		secondimu.gyro_x = thisimu.gyro_x * (1 - lambda);
		secondimu.gyro_y = thisimu.gyro_y * (1 - lambda);
		secondimu.gyro_z = thisimu.gyro_z * (1 - lambda);
		secondimu.accel_x = thisimu.accel_x * (1 - lambda);
		secondimu.accel_y = thisimu.accel_y * (1 - lambda);
		secondimu.accel_z = thisimu.accel_z * (1 - lambda);
	}
}

/*
**********************************************************
函数名：GNSS位置更新
参  数：navstate          当前的导航状态
		gnssdata          当前历元的GNSS结果数据
		kf	              卡尔曼滤波状态（x和P）
		antlever          杆臂（惯导传感器到GNSS天线中心）
函数功能：利用GNSS的位置结果对当前导航状态进行滤波更新
**********************************************************
*/
void GNSSUpdate_POS(navstate_t navstate, GNSSData_t gnssdata, kf_t& kf, Matrix antlever)
{
	if (gnssdata.PosStd[0] > 5 || gnssdata.PosStd[1] > 5 || gnssdata.PosStd[2] > 5)
	{
		//throw("WARNING: Abandon gnss position measurement at:", gnssdata.time);
	}
	else
	{
		// measurement innovation
		Matrix DR(3, 3, 0.0);
		DR(0, 0) = navstate.Rm + navstate.pva.blh.H;
		DR(1, 1) = (navstate.Rn + navstate.pva.blh.H) * cos(navstate.pva.blh.B);
		DR(2, 2) = -1;
		Matrix dP(3, 1, 0.0);
		dP(0, 0) = navstate.pva.blh.B - gnssdata.blh.B;
		dP(1, 0) = navstate.pva.blh.L - gnssdata.blh.L;
		dP(2, 0) = navstate.pva.blh.H - gnssdata.blh.H;
		Matrix Z = (DR * dP) + (navstate.Cbn * antlever);

		// measurement matrix and noise matrix
		Matrix R(3, 3, 0.0);
		for(int i=0;i<3;i++)
			R(i,i)= gnssdata.PosStd[i] * gnssdata.PosStd[i];
		Matrix H(3, kf.Rank, 0.0);
		Matrix H7_9 = (navstate.Cbn * antlever).skew();
		for (int i = 0; i < 3; i++)
		{
			H(i, i) = 1;
			for (int j = 0; j < 3; j++)
			{
				H(i, j + 6) = H7_9(i, j);
			}
		}

		// update covariance and state vector
		Matrix K = kf.P * H.transpose() * ((H * kf.P * H.transpose() + R).inverse());
		kf.x = kf.x + (K * (Z - (H * kf.x)));
		Matrix I(kf.Rank, kf.Rank, 0.0);
		I = I.eye();
		kf.P = (I - (K * H)) * kf.P * ((I - (K * H)).transpose()) + (K * R * K.transpose());
	}
}


/*
**********************************************************
函数名：GNSS速度更新
参  数：navstate          当前的导航状态
		gnssdata          当前历元的GNSS结果数据
		kf	              卡尔曼滤波状态（x和P）
		antlever          杆臂（惯导传感器到GNSS天线中心）
		thisimu           当前历元的惯导数据
		dt                当前历元与上一历元惯导的时间差
函数功能：利用GNSS的速度结果对当前导航状态进行滤波更新
**********************************************************
*/
void GNSSUpdata_Vel(navstate_t navstate, GNSSData_t gnssdata, kf_t& kf, Matrix antlever, IMUData_t thisimu, double dt)
{
	if (gnssdata.VelStd[0] > 0.5 || gnssdata.VelStd[1] > 0.5 || gnssdata.VelStd[2] > 0.5)
	{
		//throw("WARNING: Abandon gnss velocity measurement at:", gnssdata.time);
	}
	else
	{
		Matrix omega_ie_n(3, 1, 0.0);
		omega_ie_n(0, 0) = Omega_WGS * cos(navstate.pva.blh.B);
		omega_ie_n(1, 0) = 0;
		omega_ie_n(2, 0) = -Omega_WGS * sin(navstate.pva.blh.B);
		Matrix omega_en_n(3, 1, 0.0);
		omega_en_n(0, 0) = navstate.pva.Ve / (navstate.Rn + navstate.pva.blh.H);
		omega_en_n(1, 0) = -navstate.pva.Vn / (navstate.Rm + navstate.pva.blh.H);
		omega_en_n(2, 0) = -navstate.pva.Ve * tan(navstate.pva.blh.B) / (navstate.Rn + navstate.pva.blh.H);
		Matrix omega_in_n = omega_ie_n + omega_en_n;
		Matrix omega_ib_b(3, 1, 0.0);
		omega_ib_b(0, 0) = thisimu.gyro_x / dt;
		omega_ib_b(1, 0) = thisimu.gyro_y / dt;
		omega_ib_b(2, 0) = thisimu.gyro_z / dt;
		Matrix VGn(3, 1, 0.0);
		VGn(0, 0) = gnssdata.Vel[0];
		VGn(1, 0) = gnssdata.Vel[1];
		VGn(2, 0) = gnssdata.Vel[2];
		Matrix VIn(3, 1, 0.0);
		VIn(0, 0) = navstate.pva.Vn;
		VIn(1, 0) = navstate.pva.Ve;
		VIn(2, 0) = navstate.pva.Vd;
		Matrix Z = VIn - (omega_in_n.skew() * navstate.Cbn * antlever) - (navstate.Cbn * (antlever.skew() * omega_ib_b)) - VGn;

		// measurement innovation, noise, matrix
		Matrix R(3, 3, 0.0);
		for (int i = 0; i < 3; i++)
			R(i, i) = gnssdata.VelStd[i] * gnssdata.VelStd[i];
		Matrix H(3, kf.Rank, 0.0);
		Matrix H7_9 = (omega_in_n.skew() * ((navstate.Cbn * antlever).skew())) * (-1) - (navstate.Cbn * antlever.skew() * omega_ib_b).skew();
		Matrix H10_12 = (navstate.Cbn * antlever.skew()) * (-1);
		Matrix H16_18 = (navstate.Cbn * antlever.skew() * omega_ib_b.diag()) * (-1);
		for (int i = 0; i < 3; i++)
		{
			H(i, i + 3) = 1;
			for (int j = 0; j < 3; j++)
			{
				H(i, j + 6) = H7_9(i, j);
				H(i, j + 9) = H10_12(i, j);
				H(i, j + 15) = H16_18(i, j);
			}
		}

		// update covariance and state vector
		Matrix K = kf.P * H.transpose() * ((H * kf.P * H.transpose() + R).inverse());
		kf.x = kf.x + (K * (Z - (H * kf.x)));
		Matrix I(kf.Rank, kf.Rank, 0.0);
		I = I.eye();
		kf.P = (I - (K * H)) * kf.P * ((I - (K * H)).transpose()) + (K * R * K.transpose());
	}
}