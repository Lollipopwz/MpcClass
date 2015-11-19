#pragma once
#include "MpcClass.h"

public ref class Managed
{
public:
	Managed();
	Managed(int NpValue, int NcValue, double Q1, double Q2, int RValue);
	~Managed();
	void ManSendValues(double time, double Previous, double u0, double u1, double u2, double u3, double u4, double u5);
	void ManGetPath(double PointX[], double PointY[]);
	double ManCalcu();
private:
	MpcClass * mMpc;
};