#include "KF_GINS.h"
#include "AllanVariance.h"


int main()
{
	CFG_t cfg;
	vector<GNSSData_t> gnssdata = readGNSSData(cfg.GNSS_File);
	vector<IMUData_t> imudata = readIMUFile(cfg.dynamicIMU_File, cfg); //rad/s，m/s^2
	//get process time
	cfg.movetime = imudata[0].time + cfg.static_T;
	double gnssstarttime = gnssdata[0].time;
	double gnssendtime = gnssdata[gnssdata.size() - 1].time;
	double imustarttime = imudata[0].time;
	double imuendtime = imudata[imudata.size() - 1].time;
	double starttime, endtime;
	if (imustarttime > gnssstarttime)
		starttime = imustarttime;
	else
		starttime = gnssstarttime;
	if (imuendtime > gnssendtime)
		endtime = gnssendtime;
	else
		endtime = imuendtime;
	cfg.starttime = starttime;
	cfg.endtime = endtime;

	//get process index
	int imustartIndex = 0, imuendIndex = 0;
	int gnssstartIndex = 0, gnssendIndex = 0;
	do
	{
		imustartIndex += 1;
	} while (imudata[imustartIndex].time < cfg.starttime);
	do
	{
		imuendIndex += 1;
	} while (imudata[imuendIndex].time < cfg.endtime);
	do
	{
		gnssstartIndex += 1;
	} while (gnssdata[gnssstartIndex].time < cfg.starttime);
	do
	{
		gnssendIndex += 1;
	} while (gnssdata[gnssendIndex].time < cfg.endtime);

	//Allan方差分析
	vector<IMUData_t> imudata_static = readIMUFile(cfg.staticIMU_File, cfg);
	IMUList_t imulist = getIMUList(imudata_static);
	NoiseParameters imuNoise = extract_noise_parameters(imulist.grydata, imulist.accdata, 1 / cfg.Tau0, 100);

	//CFG部分参数配置
	for (int i = 0; i < 3; i++)
	{	
		cfg.initvel(i, 0) = gnssdata[gnssstartIndex].Vel[i];
		cfg.initposstd(i, 0) = gnssdata[gnssstartIndex].PosStd[i];
		cfg.initvelstd(i, 0) = gnssdata[gnssstartIndex].VelStd[i];
		
		cfg.initaccbiasstd(i, 0) = imuNoise.acc_bias_std[i];
		cfg.initgyrbiasstd(i, 0) = imuNoise.gyro_bias_std[i] * DegToRad;
		cfg.accvrw(i, 0) = imuNoise.acc_vrw[i];
		cfg.gryarw(i, 0) = imuNoise.gyro_arw[i] * DegToRad;

		cfg.accbiasstd = cfg.initaccbiasstd;
		cfg.gyrbiasstd = cfg.initgyrbiasstd;
		

		cfg.initaccscalestd(i, 0) = 0.00003;
		cfg.initgyrscalestd(i, 0) = 0.00003;

		cfg.accscalestd(i, 0) = 0.00003;
		cfg.gyrscalestd(i, 0) = 0.00003;
	}

	cfg.initpos(0, 0) = gnssdata[gnssstartIndex].blh.B;
	cfg.initpos(1, 0) = gnssdata[gnssstartIndex].blh.L;
	cfg.initpos(2, 0) = gnssdata[gnssstartIndex].blh.H;
	cfg.antlever(0, 0) = -0.01;
	cfg.antlever(1, 0) = -0.168;
	cfg.antlever(2, 0) = -0.1;

	//粗对准
	initialAlign(imudata, cfg);
	cfg.initatt(2, 0) = 32 * DegToRad;

	cfg.initattstd(0, 0) = cfg.initattstd(1, 0) = 0.1 * DegToRad;
	cfg.initattstd(2, 0) = 0.2 * DegToRad;

	//initialization 
	kf_t kf;
	navstate_t last2state, laststate, navstate;
	Initialize(kf, navstate, cfg);
	laststate = navstate;
	last2state = laststate;

	//data index preprocess
	IMUData_t thisimu = imudata[imustartIndex];
	thisimu.accel_x = (cfg.Tau0 * thisimu.accel_x - cfg.Tau0 * navstate.accbias(0, 0)) / (1 + navstate.accscale(0, 0));
	thisimu.accel_y = (cfg.Tau0 * thisimu.accel_y - cfg.Tau0 * navstate.accbias(1, 0)) / (1 + navstate.accscale(1, 0));
	thisimu.accel_z = (cfg.Tau0 * thisimu.accel_z - cfg.Tau0 * navstate.accbias(2, 0)) / (1 + navstate.accscale(2, 0));
	thisimu.gyro_x = (cfg.Tau0 * thisimu.gyro_x - cfg.Tau0 * navstate.grybias(0, 0)) / (1 + navstate.gryscale(0, 0));
	thisimu.gyro_y = (cfg.Tau0 * thisimu.gyro_y - cfg.Tau0 * navstate.grybias(1, 0)) / (1 + navstate.gryscale(1, 0));
	thisimu.gyro_z = (cfg.Tau0 * thisimu.gyro_z - cfg.Tau0 * navstate.grybias(2, 0)) / (1 + navstate.gryscale(2, 0));
	IMUData_t lastimu = thisimu;
	double imudt = thisimu.time - lastimu.time;
	int gnssindex = gnssstartIndex;
	while (gnssdata[gnssindex].time < thisimu.time)
		gnssindex += 1;
	int imuindex = imustartIndex + 1;

	ofstream outfile(cfg.Result_File);
	if (!outfile.is_open())
	{
		cerr << "open error" << endl;
	}
	// 输出表头 (列名)
	outfile << std::setw(12) << "Time(S)" << " "
		<< std::setw(12) << "Lat(Deg)" << " "
		<< std::setw(12) << "Lon(Deg)" << " "
		<< std::setw(10) << "H(m)" << " "
		<< std::setw(10) << "Vn(m/s)" << " "
		<< std::setw(10) << "Ve(m/s)" << " "
		<< std::setw(10) << "Vd(m/s)" << " "
		<< std::setw(10) << "Roll(Deg)" << " "
		<< std::setw(11) << "Pitch(Deg)" << " "
		<< std::setw(12) << "Heading(Deg)" << " "
		<< std::setw(10) << "LatStd" << " "
		<< std::setw(10) << "LonStd" << " "
		<< std::setw(8) << "HStd" << " "
		<< std::setw(10) << "VnStd" << " "
		<< std::setw(10) << "VeStd" << " "
		<< std::setw(10) << "VdStd" << " "
		<< std::setw(10) << "RollStd" << " "
		<< std::setw(11) << "PitchStd" << " "
		<< std::setw(12) << "HeadingStd" << "\n";

	/************  MAIN PROCEDD PROCEDURE  *************/
	for (imuindex; imuindex < imuendIndex; imuindex++)
	{
		//set value of last state
		lastimu = thisimu;
		last2state = laststate;
		laststate = navstate;
		thisimu = imudata[imuindex];
		imudt = thisimu.time - lastimu.time;
		//compensate IMU error
		//零偏和比例因子初值为0，在卡尔曼滤波过程中进行补偿
		thisimu.accel_x = (imudt*thisimu.accel_x - imudt * navstate.accbias(0, 0)) / (1 + navstate.accscale(0, 0));
		thisimu.accel_y = (imudt*thisimu.accel_y - imudt * navstate.accbias(1, 0)) / (1 + navstate.accscale(1, 0));
		thisimu.accel_z = (imudt*thisimu.accel_z - imudt * navstate.accbias(2, 0)) / (1 + navstate.accscale(2, 0));
		thisimu.gyro_x = (imudt*thisimu.gyro_x - imudt * navstate.grybias(0, 0)) / (1 + navstate.gryscale(0, 0));
		thisimu.gyro_y = (imudt*thisimu.gyro_y - imudt * navstate.grybias(1, 0)) / (1 + navstate.gryscale(1, 0));
		thisimu.gyro_z = (imudt*thisimu.gyro_z - imudt * navstate.grybias(2, 0)) / (1 + navstate.gryscale(2, 0));

		//adjust GNSS index
		while (gnssindex <= gnssendIndex && gnssdata[gnssindex].time < lastimu.time)
			gnssindex += 1;
		if (gnssindex > gnssendIndex)
		{
			cerr << "GNSS file END!" << endl;
			break;
		}
		GNSSData_t thisgnss;

		//determine whether gnss update is required
		
		if (lastimu.time == gnssdata[gnssindex].time)
		{

			//do gnss update for the current state
			thisgnss = gnssdata[gnssindex];
			GNSSUpdate_POS(navstate, thisgnss, kf, cfg.antlever);
			GNSSUpdata_Vel(navstate, thisgnss, kf, cfg.antlever, thisimu, imudt);
			ErrorFeedback(kf, navstate);
			gnssindex += 1;
			last2state = laststate;
			laststate = navstate;

			//do propagation for current imu data
			imudt = thisimu.time - lastimu.time;
			INSMesh(lastimu, thisimu, last2state, laststate, navstate);
			InsPropagate(navstate, thisimu, imudt, kf, cfg.corrtime);
		}
		else if (lastimu.time<gnssdata[gnssindex].time && thisimu.time>gnssdata[gnssindex].time)
		{
			//ineterpolate imu to gnss time
			IMUData_t firstimu, secondimu;
			interpolate(lastimu, thisimu, gnssdata[gnssindex].time, firstimu, secondimu);

			//do propagation for first imu
			imudt = firstimu.time - lastimu.time;
			INSMesh(lastimu, firstimu, last2state, laststate, navstate);
			InsPropagate(navstate, firstimu, imudt, kf, cfg.corrtime);

			//do gnss update
			thisgnss = gnssdata[gnssindex];
			GNSSUpdate_POS(navstate, thisgnss, kf, cfg.antlever);
			GNSSUpdata_Vel(navstate, thisgnss, kf, cfg.antlever, firstimu, imudt);
			ErrorFeedback(kf, navstate);
			gnssindex += 1;
			last2state = laststate;
			laststate = navstate;
			lastimu = firstimu;

			// % do propagation for second imu
			imudt = secondimu.time - lastimu.time;
			INSMesh(lastimu, secondimu, last2state, laststate, navstate);
			InsPropagate(navstate, secondimu, imudt, kf, cfg.corrtime);
		}
		

		else
		{
			//only do propagation
			//INS mechanization
			INSMesh(lastimu, thisimu, last2state, laststate, navstate);
			//error propagation
			InsPropagate(navstate, thisimu, imudt, kf, cfg.corrtime);
		}

		double time = navstate.time;
		double lat_deg = navstate.pva.blh.B * RadToDeg;
		double lon_deg = navstate.pva.blh.L * RadToDeg;
		double H_m = navstate.pva.blh.H;
		double Vn = navstate.pva.Vn;
		double Ve = navstate.pva.Ve;
		double Vd = navstate.pva.Vd;
		double roll_deg = navstate.pva.rph.roll * RadToDeg;
		double pitch_deg = navstate.pva.rph.pitch * RadToDeg;
		double heading = fmod(navstate.pva.rph.heading + 2 * pi, 2 * pi);
		double heading_deg = heading * RadToDeg;
		double lat_std = sqrt(kf.P(0, 0));
		double lon_std = sqrt(kf.P(1, 1));
		double H_std = sqrt(kf.P(2, 2));
		double Vn_std = sqrt(kf.P(3, 3));
		double Ve_std = sqrt(kf.P(4, 4));
		double Vd_std = sqrt(kf.P(5, 5));
		double roll_std_deg = sqrt(kf.P(6, 6)) * RadToDeg;
		double pitch_std_deg = sqrt(kf.P(7, 7)) * RadToDeg;
		double heading_std_deg = sqrt(kf.P(8, 8)) * RadToDeg;

		//结果输出
		// 设置输出精度和格式
		outfile << fixed << setprecision(6);
		outfile << std::setw(12) << time << " "
			<< std::setw(12) << lat_deg << " "
			<< std::setw(12) << lon_deg << " "
			<< std::setw(10) << H_m << " "
			<< std::setw(10) << Vn << " "
			<< std::setw(10) << Ve << " "
			<< std::setw(10) << Vd << " "
			<< std::setw(10) << roll_deg << " "
			<< std::setw(11) << pitch_deg << " "
			<< std::setw(12) << heading_deg << " "
			<< std::setw(10) << lat_std << " "
			<< std::setw(10) << lon_std << " "
			<< std::setw(8) << H_std << " "
			<< std::setw(10) << Vn_std << " "
			<< std::setw(10) << Ve_std << " "
			<< std::setw(10) << Vd_std << " "
			<< std::setw(10) << roll_std_deg << " "
			<< std::setw(11) << pitch_std_deg << " "
			<< std::setw(12) << heading_std_deg << "\n";
	}
	outfile.close();
	return 0;
}
