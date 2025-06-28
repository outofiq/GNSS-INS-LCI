#include <string>
#include "coordinate.h"

struct CFG_t
{
	string          staticIMU_File;    //��̬�ߵ������ļ�
	string          dynamicIMU_File;   //��̬�ߵ������ļ�
	string          GNSS_File;         //GNSS��λ��������ļ�
	string          Result_File;       //�������ļ�

	const double    acc_scale;         //���ٶȼ�ת��ϵ������λ��m/s^2/LSB��
	const double    gry_scale;         //������ת��ϵ������λ����/s/LSB��
	const double    Tau0;              //�����������λ��s��
	double          static_T;          //��̬�ɼ���ʼ��ֹ��ʱ��
	double          movetime;          //С����ʼ�ƶ���ʱ��
	double          starttime;         
	double          endtime;

	Matrix          initpos;           //��ʼλ��(rad,rad,m)
	Matrix          initvel;           //��ʼ�ٶȣ�m/s��
	Matrix          initatt;           //��ʼ��̬��rad��
	Matrix          initposstd;        //��ʼλ��std
	Matrix          initvelstd;        //��ʼ�ٶ�std
	Matrix          initattstd;        //��ʼ��̬std

	Matrix          initgrybias;       //��ʼ������ƫ��rad/s��
	Matrix          initaccbias;       //��ʼ���ٶȼ���ƫ��m/s^2��
	Matrix          initgryscale;      //���ݱ�������
	Matrix          initaccscale;      //���ٶȼƱ�������
	Matrix          initgyrbiasstd;    //��ʼ������ƫstd
	Matrix          initaccbiasstd;    //��ʼ���ٶȼ���ƫstd
	Matrix          initgyrscalestd;   //��ʼ���ݱ�������std
	Matrix          initaccscalestd;   //��ʼ���ٶȼƱ�������std
	Matrix          gyrbiasstd;        //��ʼ������ƫstd
	Matrix          accbiasstd;        //��ʼ���ٶȼ���ƫstd
	Matrix          gyrscalestd;       //��ʼ���ݱ�������std
	Matrix          accscalestd;       //��ʼ���ٶȼƱ�������std

	Matrix          gryarw;            //�Ƕ��������
	Matrix          accvrw;            //�����������

	Matrix          zupt_noise;        //���ٹ۲�����
	Matrix          zaru_noise;        //������ʹ۲�����

	Matrix          antlever;          //��װ�ǣ�m��

	double          corrtime;          //���ʱ��
	CFG_t() :
		staticIMU_File("E:\\�ۺ�ʵϰ\\data\\static-17.ASC"),
		dynamicIMU_File("E:\\�ۺ�ʵϰ\\data\\dynamic-17.ASC"),
		GNSS_File("E:\\�ۺ�ʵϰ\\GNSSData\\GNSSData\\Export\\dynamic-17.pos"),
		Result_File("E:\\�ۺ�ʵϰ\\Result\\Result.txt"),

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