#ifndef CONTROLLER_H
#define CONTROLLER_H
#include <QuadProg++.hh>
#include <Eigen/Dense>
// #include <Eigen/LU>

class MyController
{

public:
	double U[1];
	double t;

	MyController();
	~MyController();
	_declspec(dllexport) void SendValues(double time, double Previous, double u0, double u1, double u2, double u3, double u4, double u5);
	_declspec(dllexport) Eigen::MatrixXd GetDesignPath(double x, double y, double PointX[], double PointY[]);
	_declspec(dllexport) double Calculate();

private:
	Eigen::MatrixXd Yita_ref = Eigen::MatrixXd::Zero(2, 2);
	//主要状态量
	double y_dot;
	double x_dot;
	double phi;
	double phi_dot;
	double Y;
	double X;
	double Y_dot;
	double X_dot;


	//时间设置
	const double T_inter;
	const double T_all;
	
	//Constraints Setting
	const double umin;
	const double umax;
	const double delta_umin;
	const double delta_umax;

	//算法变量初始化
	const int Nx;
	const int Nu;
	const int Ny;
	const int Np;
	const int Nc;
	const int Row;

	//车辆参数
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

	//用于生成双移线轨迹的参数
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