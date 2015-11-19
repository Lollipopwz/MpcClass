#include "stdafx.h"
#include "Managed.h"
Managed::Managed()
{
	mMpc = new MpcClass();
}
Managed::Managed(int NpValue, int NcValue, double Q1, double Q2, int RValue)
{
	mMpc = new MpcClass(NpValue, NcValue, Q1, Q2, RValue);
}
Managed::~Managed()
{
	delete mMpc;
}
void Managed::ManSendValues(double time, double Previous, double u0, double u1, double u2, double u3, double u4, double u5)
{
	mMpc->SendValues(time, Previous, u0, u1, u2, u3, u4, u5);
}
void Managed::ManGetPath(double PointX[], double PointY[])
{
	mMpc->GetDesignPath(PointX, PointY);
}
double Managed::ManCalcu()
{
	return mMpc->Calculate();
}