#include"Ctime.h"
using namespace std;

CTimer::CTimer(void)
{
}


CTimer::~CTimer(void)
{
}

int CTimer::time_in()
{
	QueryPerformanceFrequency(&litmp);//���ʱ��Ƶ��
	dff = (double)litmp.QuadPart;

	QueryPerformanceCounter(&litmp);//��ó�ʼֵ
	qt1 = litmp.QuadPart;

	return 1;
}

double CTimer::time_out()
{
	QueryPerformanceCounter(&litmp);//�����ֵֹ
	qt2 = litmp.QuadPart;

	dfm = (double)(qt2 - qt1);
	dft = dfm / dff;//��ö�Ӧ��ʱ��ֵ

	return dft;
}