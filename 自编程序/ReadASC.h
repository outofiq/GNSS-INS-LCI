#include <fstream>
#include <sstream>
#include "Config.h"

const double DegToRad = 0.01745329252;
const double RadToDeg = 57.29577951308;

// ����IMU���ݽṹ��
struct IMUData_t 
{
    double time;          // GPSʱ�䣨�룩
    double accel_x;       // ���ٶȼ�X��(m/s^2)
    double accel_y;       // Y��
    double accel_z;       // Z��
    double gyro_x;        // ������X��(rad/s)
    double gyro_y;        // Y��
    double gyro_z;        // Z��
};

struct IMUList_t
{
    vector<vector<double>> accdata;
    vector<vector<double>> grydata;
};

//��ȡASC��ʽ�Ĺߵ�����
vector<IMUData_t> readIMUFile(const string& filepath, CFG_t cfg);
//����ȡ���ļ��ٶȼƺ����������ݰ�����ϵ�ֱ����ɸ��Ե�����
IMUList_t getIMUList(const vector<IMUData_t> imuData);