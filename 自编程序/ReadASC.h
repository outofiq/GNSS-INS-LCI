#include <fstream>
#include <sstream>
#include "Config.h"

const double DegToRad = 0.01745329252;
const double RadToDeg = 57.29577951308;

// 定义IMU数据结构体
struct IMUData_t 
{
    double time;          // GPS时间（秒）
    double accel_x;       // 加速度计X轴(m/s^2)
    double accel_y;       // Y轴
    double accel_z;       // Z轴
    double gyro_x;        // 陀螺仪X轴(rad/s)
    double gyro_y;        // Y轴
    double gyro_z;        // Z轴
};

struct IMUList_t
{
    vector<vector<double>> accdata;
    vector<vector<double>> grydata;
};

//读取ASC格式的惯导数据
vector<IMUData_t> readIMUFile(const string& filepath, CFG_t cfg);
//将读取到的加速度计和陀螺仪数据按照轴系分别生成各自的序列
IMUList_t getIMUList(const vector<IMUData_t> imuData);