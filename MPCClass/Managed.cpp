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
void Managed::ManSendValues(double TInter, double Previous, double XDot, double Phi, double PhiDot)
{
	mMpc->SendValues( TInter,  Previous,  XDot,  Phi,  PhiDot);
}
void Managed::ManGetPath(double PointX[], double PointY[])
{
	mMpc->GetDesignPath(PointX, PointY);
}
double Managed::ManCalcu()
{
	return mMpc->Calculate();
}