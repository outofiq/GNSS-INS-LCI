#include "AllanVariance.h"

/*
**********************************************************
函数名：Allan方差计算
参  数：omega           要计算Allan方差的数据序列
        fs              该序列的采样频率
        pts             最终生成的Allan方差的点数
返回值：Allan方差序列
函数功能：计算传入序列的Allan方差，生成pts个点
**********************************************************
*/
AllanResult allan(const vector<double>& omega, double fs, int pts)
{
	AllanResult result;
	const int N = omega.size(); // 样本数

	// 生成聚类因子m
    vector<int> m;
	{
        int maxCluster = (int)(N / 2);
        vector<double> logspace;
        for (int i = 0; i < pts; ++i) 
        {
            logspace.push_back(pow(10, i * log10(maxCluster) / (pts - 1)));
        }
        sort(logspace.begin(), logspace.end());
        for (auto& v : logspace) 
        {
            int val = (int)ceil(v);
            if (val <= maxCluster && (m.empty() || val != m.back()))
                m.push_back(val);
        }
	}

    // 计算θ
   vector<double> theta(N);
   double sum = 0;
   for (int i = 0; i < N; ++i)
   {
       sum += omega[i];
       theta[i] = sum / fs;
   }

   // Allan方差计算
   result.tau.resize(m.size());
   result.sigma.resize(m.size());

   for (int i = 0; i < m.size(); ++i) 
   {
       int mi = m[i];
       result.tau[i] = mi / fs;

       double sum = 0;
       for (int k = 0; k < N - 2 * mi; ++k)
       {
           double diff = theta[k + 2 * mi] - 2 * theta[k + mi] + theta[k];
           sum += diff * diff;
       }
       double sigma2 = sum / (2 * pow(mi / fs, 2) * (N - 2 * mi));
       result.sigma[i] = sqrt(sigma2);
   }

   return result;
}


// 在指定tau范围内查找最小sigma值
double find_min_sigma_in_range(const AllanResult& res, double tau_low, double tau_high) {
    double min_sigma = std::numeric_limits<double>::max();
    bool found = false;

    for (size_t i = 0; i < res.tau.size(); ++i) {
        if (res.tau[i] >= tau_low && res.tau[i] <= tau_high) {
            if (res.sigma[i] < min_sigma) {
                min_sigma = res.sigma[i];
                found = true;
            }
        }
    }

    // 如果指定范围内没找到，扩大搜索范围
    if (!found) {
        for (size_t i = 0; i < res.tau.size(); ++i) {
            if (res.tau[i] >= tau_low) {
                if (res.sigma[i] < min_sigma) {
                    min_sigma = res.sigma[i];
                }
            }
        }
    }

    return min_sigma;
}

// 提取噪声参数的主函数
NoiseParameters extract_noise_parameters(
    const std::vector<std::vector<double>>& acc_data,  // [m/s²] 三轴加速度数据
    const std::vector<std::vector<double>>& gyro_data, // [rad/s] 三轴陀螺仪数据
    double fs,                                         // 采样频率 (Hz)
    int pts = 100                                      // Allan方差计算点数
) {
    NoiseParameters params;
    std::vector<double> acc_vrw, gyro_arw, gyro_bias, acc_bias;

    // 处理加速度计数据 (每轴)
    for (const auto& data : acc_data) {
        AllanResult res = allan(data, fs, pts);

        // 提取VRW (τ=1s处的值)
        double vrw = 0.0;
        double min_diff = std::numeric_limits<double>::max();
        for (size_t i = 0; i < res.tau.size(); ++i) {
            double diff = std::abs(res.tau[i] - 1.0);
            if (diff < min_diff) {
                min_diff = diff;
                vrw = res.sigma[i];
            }
        }
        acc_vrw.push_back(vrw);

        // 提取加速度零偏不稳定性 (τ=100-1000s范围内的最小值)
        double bias_std = find_min_sigma_in_range(res, 100.0, 1000.0);
        acc_bias.push_back(bias_std);
    }

    // 处理陀螺仪数据 (每轴)
    for (const auto& data : gyro_data) {
        AllanResult res = allan(data, fs, pts);

        // 提取ARW (τ=1s处的值)
        double arw = 0.0;
        double min_diff = std::numeric_limits<double>::max();
        for (size_t i = 0; i < res.tau.size(); ++i) {
            double diff = std::abs(res.tau[i] - 1.0);
            if (diff < min_diff) {
                min_diff = diff;
                arw = res.sigma[i];
            }
        }
        gyro_arw.push_back(arw);

        // 提取陀螺零偏不稳定性 (τ=10-100s范围内的最小值)
        double bias_std = find_min_sigma_in_range(res, 10.0, 100.0);
        gyro_bias.push_back(bias_std);
    }

    // 计算三轴平均值作为最终参数
    params.acc_vrw = acc_vrw;
    params.gyro_arw = gyro_arw;
    params.gyro_bias_std = gyro_bias;
    params.acc_bias_std = acc_bias;

    return params;
}