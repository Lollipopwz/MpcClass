#pragma once
#include "MpcClass.h"

public ref class Managed
{
public:
	Managed();
	Managed(int NpValue, int NcValue, double Q1, double Q2, int RValue);
	~Managed();
	void ManSendValues(double TInter, double Previous, double XDot, double Phi, double PhiDot);
	void ManGetPath(double PointX[], double PointY[]);
	double ManCalcu();
private:
	MpcClass * mMpc;
};