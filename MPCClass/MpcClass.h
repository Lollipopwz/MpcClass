#ifndef CONTROLLER_H
#define CONTROLLER_H
#define PI 3.141592653589793
#define EXPORT _declspec(dllexport)
#include "QuadProg++.hh"
// #include <Eigen/Dense>
#include <Core>
class EXPORT MpcClass
{
public:

	MpcClass();
	MpcClass(int NpValue,int NcValue,double Q1,double Q2,int RValue);
	~MpcClass();
	void SendValues(double T_in, double time, double Previous, double u0, double u1, double u2, double u3, double u4, double u5);
	void GetDesignPath(double x, double y, double PointX[], double PointY[]);
	double Calculate();

private:

	double U[1];
	double t;	
	Eigen::MatrixXd Yita_ref;
	//Ŀǰ���޸Ĳ���
	double Q1;
	double Q2;
	double RValue;
	//��Ҫ״̬��
	double y_dot;
	double x_dot;
	double phi;
	double phi_dot;
	double Y;
	double X;
	double Y_dot;
	double X_dot;


	//ʱ������
	double T_inter;
	const double T_all;
	
	//Constraints Setting
	const double umin;
	const double umax;
	const double delta_umin;
	const double delta_umax;
	const int IndexTh;

	//�㷨������ʼ��
	const int Nx;
	const int Nu;
	const int Ny;
	const int Np;
	const int Nc;
	const int Row;

	//��������
	const double Sf;
	const double Sr;
	const double lf;
	const double lr;
	const double Ccf;
	const double Ccr;
	const double Clf;
	const double Clr;
	const double Mass;
	const double Gravity;
	const double Inertia;

	//��������˫���߹켣�Ĳ���
	const double shape;
	const double Dx1;
	const double Dx2;
	const double Dy1;
	const double Dy2;
	const double Xs1;
	const double Xs2;

	//Road shap parameters
	const double road_amp;
	const double road_fre;

	Eigen::MatrixXd Kron(Eigen::MatrixXd A, Eigen::MatrixXd B);
	void Diag_Mat(Eigen::MatrixXd &rst, Eigen::MatrixXd innerMat, int num);
	Eigen::MatrixXd Pow_Mat(Eigen::MatrixXd x, int y);
	void Mer_Mat(Eigen::MatrixXd &dsttmp, Eigen::MatrixXd &src11, Eigen::MatrixXd &src12, Eigen::MatrixXd &src21, Eigen::MatrixXd &src22);
	void Mer_Vertical(Eigen::MatrixXd &dsttmp, Eigen::MatrixXd &src1, Eigen::MatrixXd &src2, Eigen::MatrixXd &src3, Eigen::MatrixXd &src4);
	void Mer_Demo(void);
	void EigenToMatrix(Eigen::MatrixXd &src, QuadProgPP::Matrix<double> &dst);
	void EigenToVector(Eigen::MatrixXd &src, QuadProgPP::Vector<double> &dst);

};
#endif