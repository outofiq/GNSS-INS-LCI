#include "ReadASC.h"


/*
**********************************************************
����������ȡ�ߵ�����
��  ����filepath          IMU�����ļ���
        cfg               ���ò����ṹ��
����ֵ���ߵ��������У�����Ԫ��
�������ܣ���ȡASC�ļ�������ߵ�����
**********************************************************
*/
vector<IMUData_t> readIMUFile(const string& filepath, CFG_t cfg)
{
	vector<IMUData_t> data;
    double acc_x, acc_y, acc_z;
    double gry_x, gry_y, gry_z;
	ifstream file(filepath);
    if (!file.is_open())
    {
        cerr << "open error!";
    }
	string line;

    while (getline(file, line)) 
    {
        // �滻�ֺ�Ϊ�����Ա�ͳһ�ָ�
        size_t semicolon_pos = line.find(';');
        size_t asterisk_pos = line.find('*');
        if (semicolon_pos != string::npos) 
        {
            line.replace(semicolon_pos, 1, ",");
            line.replace(asterisk_pos, 1, ",");
        }

        stringstream ss(line);
        string token;
        vector<string> tokens;

        // �ָ�ŷָ�������
        while (getline(ss, token, ',')) 
        {
            tokens.push_back(token);
        }

        // ȷ�����㹻�ֶΣ�����15�У�
        if (tokens.size() < 15 || tokens[0]!="%RAWIMUSXA") continue;

        try 
        {
            IMUData_t entry;
            entry.time = stod(tokens[2]); // �߾���GPSʱ��
            acc_x = stod(tokens[10]) * cfg.acc_scale;
            acc_y = -stod(tokens[9]) * cfg.acc_scale;
            acc_z = stod(tokens[8]) * cfg.acc_scale;
            gry_x = stod(tokens[13]) * cfg.gry_scale * DegToRad;
            gry_y = -stod(tokens[12]) * cfg.gry_scale * DegToRad;
            gry_z = stod(tokens[11]) * cfg.gry_scale * DegToRad;

            //��ϵ����
            entry.accel_x = acc_y;
            entry.accel_y = acc_x;
            entry.accel_z = -acc_z;
            entry.gyro_x = gry_y;
            entry.gyro_y = gry_x;
            entry.gyro_z = -gry_z;

            data.push_back(entry);
        }
        catch (...) 
        {
            std::cerr << "Error parsing line: " << line << std::endl;
        }
    }
    file.close();
    return data;
}


/*
**********************************************************
�����������ߵ����ݰ���ϵ������������
��  ����imudata          IMU����
����ֵ���ߵ��������нṹ��
�������ܣ�����ȡ���Ĺߵ����ݰ���ϵ�ָ�ɸ��Ե����У����ڽ���Allan�������
**********************************************************
*/
IMUList_t getIMUList(const vector<IMUData_t> imuData)
{
    IMUList_t list;
    vector<double> acc_x;
    vector<double> acc_y;
    vector<double> acc_z;
    vector<double> gry_x;
    vector<double> gry_y;
    vector<double> gry_z;
    for (const auto& entry : imuData)
    {
        acc_x.push_back(entry.accel_x);
        acc_y.push_back(entry.accel_y);
        acc_z.push_back(entry.accel_z);
        gry_x.push_back(entry.gyro_x);
        gry_y.push_back(entry.gyro_y);
        gry_z.push_back(entry.gyro_z);
    }
    list.accdata.push_back(acc_x);
    list.accdata.push_back(acc_y);
    list.accdata.push_back(acc_z);
    list.grydata.push_back(gry_x);
    list.grydata.push_back(gry_y);
    list.grydata.push_back(gry_z);

    return list;
}
