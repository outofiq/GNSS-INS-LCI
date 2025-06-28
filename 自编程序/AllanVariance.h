#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <limits>
using namespace std;

struct AllanResult
{
	vector<double> tau;
	vector<double> sigma;
};

// 噪声参数结构体
struct NoiseParameters 
{
    vector<double> acc_vrw;         // 加速度计速度随机游走 (m/s/√Hz)
    vector<double> gyro_arw;        // 陀螺仪角度随机游走 (rad/s/√Hz)
    vector<double> gyro_bias_std;   // 陀螺仪零偏不稳定性 (rad/s)
    vector<double> acc_bias_std;    // 加速度计零偏不稳定性 (m/s²)
};

AllanResult allan(const vector<double>& omega, double fs, int pts);
double find_min_sigma_in_range(const AllanResult& res, double tau_low, double tau_high);
NoiseParameters extract_noise_parameters(
    const std::vector<std::vector<double>>& acc_data,  // [m/s²] 三轴加速度数据
    const std::vector<std::vector<double>>& gyro_data, // [rad/s] 三轴陀螺仪数据
    double fs,                                         // 采样频率 (Hz)
    int pts                                            // Allan方差计算点数
);