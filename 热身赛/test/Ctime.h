// #include<Windows.h>
#pragma once
class CTimer
{
public:
	CTimer(void);
	~CTimer(void);

	int time_in();
	double time_out();

private:
	LARGE_INTEGER litmp;
	LONGLONG qt1, qt2;
	double dft, dff, dfm;
};