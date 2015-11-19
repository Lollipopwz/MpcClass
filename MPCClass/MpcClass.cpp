#include "stdafx.h"
#include <iostream>
#include "MpcClass.h"

using namespace std;
using namespace Eigen;


MpcClass::MpcClass()
:
y_dot(0), x_dot(0), phi(0), phi_dot(0), Y(0), X(0), Y_dot(0), X_dot(0),
road_amp(0.5), road_fre(48),
Nx(6), Nu(1), Ny(2), Row(1000),Np(20),Nc(5),Q1(0),Q2(0),RValue(0),//Original Value Np 20,Nc 5.
T_inter(0.02), T_all(60),
umin(-1.744), umax(1.744), delta_umin(1.48), delta_umax(1.48),
Sf(0.2), Sr(0.2), lf(1.232), lr(1.468), Ccf(66900), Ccr(62700), Clf(66900), Clr(62700), Mass(1732), Gravity(9.8), Inertia(4175),
shape(2.4), Dx1(25), Dx2(21.95), Dy1(4.05), Dy2(5.7), Xs1(27.19), Xs2(56.46)
{}

MpcClass::MpcClass(int NpValue, int NcValue, double Q1, double Q2, int RValue)
:
Np(NpValue), Nc(NcValue), Q1(Q1), Q2(Q2), RValue(RValue),
x_dot(0), phi(0), phi_dot(0), 
road_amp(0.5), road_fre(48),
Nx(6), Nu(1), Ny(2), Row(1000),//Original Value Np 20,Nc 5.
T_inter(0.02), T_all(60),
umin(-1.744), umax(1.744), delta_umin(1.48), delta_umax(1.48),
Sf(0.2), Sr(0.2), lf(1.232), lr(1.468), Ccf(66900), Ccr(62700), Clf(66900), Clr(62700), Mass(1732), Gravity(9.8), Inertia(4175),
shape(2.4), Dx1(25), Dx2(21.95), Dy1(4.05), Dy2(5.7), Xs1(27.19), Xs2(56.46)
{}

void MpcClass::SendValues(double Previous, double Xdot, double Phi, double PhiDot)
{
	U[0] = Previous;//上一次输出
	x_dot = Xdot;
	phi = Phi;
	phi_dot = PhiDot;
}
/****************************************************************************************************/
/* 算法描述：通过当前XY坐标获取路径上对应参考点，以其为起始点，向前取Np个点，间距为X_Dot*T_inter    */
/****************************************************************************************************/
void MpcClass::GetDesignPath(double PointX[], double PointY[])
{
	//实际使用中，输入量就是部分参考路径在车体坐标系里边的坐标

	Yita_ref = MatrixXd::Zero(Np, 2);
	MatrixXd Yita_ref_Homo = MatrixXd::Ones(Yita_ref.rows(), Yita_ref.cols() + 1);
	MatrixXd Homo(3, 3);
	MatrixXd Yita_ref_cell(1, 2);

	for (int i = 0; i < Np; i++)
	{
		Yita_ref_cell << PointY[i], PointX[i];
		Yita_ref.block(i, 0, 1, 2) = Yita_ref_cell;
	}
	Yita_ref_Homo.block(0, 0, Yita_ref.rows(), Yita_ref.cols()) = Yita_ref;
	Homo << cos(phi), sin(phi), 0,
		-sin(phi), cos(phi), 0,
		X*sin(phi) - Y*cos(phi), -X*cos(phi) - Y*sin(phi), 1;
	Yita_ref_Homo = Yita_ref_Homo * Homo; //Rotate & Translation
	Yita_ref.resize(2 * Np, 1);
	for (int i = 0; i < Np; i++)
	{
		Yita_ref.block(2 * i, 0, 2, 1) = Yita_ref_Homo.block(i, 0, 1, 2).transpose();
	}
}

double MpcClass::Calculate()
{
	//计算Kesi,将状态量与控制量结合在一起
	MatrixXd Kesi(Nx + Nu, 1);
	Kesi(0, 0) = -x_dot*(tan(U[0] * PI * PI / 180));//Original value y_dot
	Kesi(1, 0) = x_dot;
	Kesi(2, 0) = phi;
	Kesi(3, 0) = phi_dot;
	Kesi(4, 0) = 0;
	Kesi(5, 0) = 0;
	Kesi(6, 0) = U[0];

	double delta_f = U[0];

	//权重矩阵的设置
	MatrixXd Q_cell(2, 2);
	Q_cell << Q1, 0, 0, Q2;
	MatrixXd Q = MatrixXd::Zero(Np * 2, Np * 2);
	Diag_Mat(Q, Q_cell, Np);

	MatrixXd R = MatrixXd::Identity(Nu*Nc, Nu*Nc);
	R = R*(5 * pow(10, RValue));//Original value： 5 * pow(10,5)

	MatrixXd a(6, 6);
	a <<
		(double)(1 - (259200 * T_inter) / (1723 * x_dot)),
		(double)(-T_inter*(phi_dot + (2 * ((460218 * phi_dot) / 5 - 62700 * y_dot)) / (1723 * pow(x_dot, 2)) - (133800 * ((154 * phi_dot) / 125 + y_dot)) / (1723 * pow(x_dot, 2)))),
		0,
		(double)(-T_inter*(x_dot - 96228 / (8615 * x_dot))),
		0,
		0,//First row

		(double)(T_inter*(phi_dot - (133800 * delta_f) / (1723 * x_dot))),
		(double)((133800 * T_inter*delta_f*((154 * phi_dot) / 125 + y_dot)) / (1723 * pow(x_dot, 2)) + 1),
		0,
		(double)(T_inter*(y_dot - (824208 * delta_f) / (8615 * x_dot))),
		0,
		0,//Second row
		0, 0, 1, T_inter, 0, 0,//Third row
		(double)((33063689036759 * T_inter) / (7172595384320 * x_dot)),
		(double)(T_inter*(((2321344006605451863 * phi_dot) / 8589934592000 - (6325188028897689 * y_dot) / 34359738368) / (4175 * pow(x_dot, 2)) + (5663914248162509 * ((154 * phi_dot) / 125 + y_dot)) / (143451907686400 * pow(x_dot, 2)))),
		0,
		(double)(1 - (813165919007900927 * T_inter) / (7172595384320000 * x_dot)),
		0,
		0,//Fourth row
		T_inter, 0, 0, 0, 1, 0,//Fifth row
		0, T_inter, 0, 0, 0, 1;//Sixth row

	MatrixXd b(6, 1);
	b <<
		(double)(133800 * T_inter / 1723),
		(double)(T_inter*((267600 * delta_f) / 1723 - (133800 * ((154 * phi_dot) / 125 + y_dot)) / (1723 * x_dot))),
		0,
		(double)(5663914248162509 * T_inter / 143451907686400),
		0,
		0;
	// 	cout << "Matrix b:\n" << b << endl;
	MatrixXd d_k = MatrixXd::Zero(Nx, 1);
	MatrixXd state_k1 = d_k;
	state_k1(0, 0) = y_dot + T_inter*(-x_dot*phi_dot + 2 * (Ccf*(delta_f - (y_dot + lf*phi_dot) / x_dot) + Ccr*(lr*phi_dot - y_dot) / x_dot) / Mass);
	state_k1(1, 0) = x_dot + T_inter*(y_dot*phi_dot + 2 * (Clf*Sf + Clr*Sr + Ccf*delta_f*(delta_f - (y_dot + phi_dot*lf) / x_dot)) / Mass);
	state_k1(2, 0) = phi + T_inter*phi_dot;
	state_k1(3, 0) = phi_dot + T_inter*((2 * lf*Ccf*(delta_f - (y_dot + lf*phi_dot) / x_dot) - 2 * lr*Ccr*(lr*phi_dot - y_dot) / x_dot) / Inertia);
	state_k1(4, 0) = T_inter * y_dot;
	state_k1(5, 0) = T_inter * x_dot;
	d_k = state_k1 - a*(Kesi.block(0, 0, 6, 1)) - b*Kesi(6, 0);

	// 	cout << "状态误差：\n" << d_k << endl;
	MatrixXd d_piao_k = MatrixXd::Zero(Nx + Nu, 1);
	d_piao_k.block(0, 0, 6, 1) = d_k.block(0, 0, Nx, 1);
	d_piao_k(6, 0) = 0;
	// 	cout << "d_piao_k: " << d_piao_k << endl;
	MatrixXd a21 = MatrixXd::Zero(Nu, Nx);
	MatrixXd a22 = MatrixXd::Identity(Nu, Nu);
	MatrixXd A;
	Mer_Mat(A, a, b, a21, a22);
	// 	cout << "Matrix A:\n" << A << endl;

	MatrixXd B(b.rows() + 1, 1);
	B.block(0, 0, b.rows(), b.cols()) = b.block(0, 0, b.rows(), b.cols());
	B(b.rows(), 0) = 1;
	//  cout << "Matrix B:\n" << B << endl;

	MatrixXd C(2, 7);
	C << 0, 0, 0, 0, 1, 0, 0,
		0, 0, 0, 0, 0, 1, 0;

	MatrixXd PHI(Np*d_piao_k.rows(), 1);
	for (int i = 0; i < Np; i++)
	{
		PHI.block(i*d_piao_k.rows(), 0, d_piao_k.rows(), 1) = d_piao_k;
	}

	// 	cout << "Matrix PHI:\n" << PHI << endl;

	MatrixXd GAMMA = MatrixXd::Zero(Np*C.rows(), Np*C.cols());//C.size:2*7
	//下三角元胞数组
	for (int i = 0; i < Np; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			GAMMA.block(i*C.rows(), j*C.cols(), C.rows(), C.cols()) = C*Pow_Mat(A, i - j);
		}
	}
	//	cout << "Matrix GAMMA:\n" << GAMMA << endl;

	MatrixXd PSI = MatrixXd::Zero(Np*C.rows(), C.cols());//size(PSI)=[Ny*Np,Nx*Nu]
	for (int i = 0; i < Np; i++)
	{
		PSI.block(i*C.rows(), 0, C.rows(), C.cols()) = C*Pow_Mat(A, i + 1);
	}
	// 	cout << "Matrix PSI:\n" << PSI.rows()<<" X "<<PSI.cols() << endl;

	MatrixXd THETA = MatrixXd::Zero(Np*Ny, Nc*Nu);//size(THETA)=[Ny*Np,Nu*Nc]
	for (int i = 0; i < Np; i++)
	{
		for (int j = 0; j < Nc; j++)
		{
			if (j <= i)THETA.block(i*Ny, j*Nu, Ny, Nu) = C*Pow_Mat(A, (int)(i - j))*B;
		}
	}
	// 	cout << "Matrix THETA:\n" << THETA << endl;

	//Get the Matrix H
	MatrixXd H;
	H = THETA.transpose() * Q * THETA + R;

#if 0

	MatrixXd H;
	MatrixXd H11 = THETA.transpose()*Q*THETA + R;
	MatrixXd H12 = MatrixXd::Zero(Nu*Nc, 1);
	MatrixXd H21 = MatrixXd::Zero(1, Nu*Nc);
	MatrixXd H22(1, 1);
	H22 << Row;
	Mer_Mat(H, H11, H12, H21, H22);//Merge the Matrix H

#endif

	/////////////////////////////////Get the Matrix f////////////////////////////////


	MatrixXd errorMat = MatrixXd::Zero(Ny*Np, 1);

	errorMat = Yita_ref - PSI*Kesi - GAMMA*PHI;

	MatrixXd f;
	f = 2 * errorMat.transpose() * Q * THETA;
	f = -f;
#if 0
	MatrixXd f11;
	f11 = 2 * errorMat.transpose()*Q*THETA;
	MatrixXd f(1, f11.cols() + 1);
	f.block(0, 0, 1, f11.cols()) = f11;
	f(0, f11.cols()) = 0;
	f = -f;

#endif

	//////////////////////////////////Generation of constraint////////////////////////////////////////

	MatrixXd A_t = MatrixXd::Zero(Nc, Nc);
	for (int i = 0; i < Nc; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			A_t(i, j) = 1;
		}
	}

	MatrixXd A_I = Kron(A_t, MatrixXd::Identity(Nu, Nu));
	MatrixXd Ut = MatrixXd::Ones(Nc, 1)* U[0];
	/////////////////////////////////Constraint of control/////////////////////////////////////////
	MatrixXd Umin = MatrixXd::Ones(Nc, 1)*umin;
	MatrixXd Umax = MatrixXd::Ones(Nc, 1)*umax;

	////////////////////////////////////Constraint of output//////////////////////////////////////
	MatrixXd ycmax(2, 1);//y x
	ycmax << 10, 10;
	MatrixXd ycmin(2, 1);
	ycmin << -10, -10;
	MatrixXd Ycmax = Kron(MatrixXd::Ones(Np, 1), ycmax);
	MatrixXd Ycmin = Kron(MatrixXd::Ones(Np, 1), ycmin);

	////////////////////////////////////Combine the output and the control//////////////////////////////////////
	MatrixXd A_cons_ori = MatrixXd::Zero(2 * (A_I.rows() + THETA.rows()), A_I.cols());
#if 0
	MatrixXd A_cons_ori = MatrixXd::Zero(2 * (A_I.rows() + THETA.rows()), A_I.cols() + 1);
#endif

	A_cons_ori.block(0, 0, A_I.rows(), A_I.cols()) = A_I;
	A_cons_ori.block(A_I.rows(), 0, A_I.rows(), A_I.cols()) = -A_I;
	A_cons_ori.block(2 * A_I.rows(), 0, THETA.rows(), THETA.cols()) = THETA;
	A_cons_ori.block(2 * A_I.rows() + THETA.rows(), 0, THETA.rows(), THETA.cols()) = -THETA;
	// 	cout << "Matrix A_cons_ori:\n" << A_cons_ori << endl;
	// 	cout << "Size of A_cons_pri:\n" << A_cons_ori.rows() << " X " << A_cons_ori.cols() << endl;

	MatrixXd B_cons_ori;
	Mer_Vertical(B_cons_ori, (MatrixXd)(Umax - Ut), (MatrixXd)(-Umin + Ut), (MatrixXd)(Ycmax - PSI*Kesi - GAMMA*PHI), (MatrixXd)(-Ycmin + PSI*Kesi + GAMMA*PHI));
	// 	cout << "Matrix B_cons:\n" << B_cons_ori << endl;
	// 	cout << "Size of B_cons:\n" << B_cons_ori.rows() << " X " << B_cons_ori.cols() << endl;


	////////////////////////////////////Constraint of states//////////////////////////////////////
	double M = 10;
	//MatrixXd delta_Umin = Kron(MatrixXd::Ones(Nc, 1), delta_umin);
	MatrixXd delta_Umin = MatrixXd::Ones(Nc, 1)*delta_umin;
	MatrixXd delta_Umax = MatrixXd::Ones(Nc, 1)*delta_umax;
	MatrixXd lb(delta_Umin.rows() + 1, 1);
	lb << delta_Umin, 0;
	MatrixXd ub(delta_Umax.rows() + 1, 1);
	ub << delta_Umax, M;

	////////////////////////////////////Combination of A_cons & B_cons//////////////////////////////////////
	MatrixXd EmptyMat;
	MatrixXd ub_Identity = MatrixXd::Identity(Nc, Nc);
	MatrixXd lb_Identity = -MatrixXd::Identity(Nc, Nc);
	lb = -lb;

	//注释以下两句，去掉状态量约束
	Mer_Vertical(A_cons_ori, A_cons_ori, lb_Identity, ub_Identity, EmptyMat);
	Mer_Vertical(B_cons_ori, B_cons_ori, delta_Umin, delta_Umax, EmptyMat);


	// 	cout << "Test the A_cons_ori:\n" << A_cons_ori.rows() << " X " << A_cons_ori.cols() << endl;
	// 	cout << "Test the B_cons_ori:\n" << B_cons_ori.rows() << " X " << B_cons_ori.cols() << endl;
	// 	cout << "矩阵lb:\n" << lb << endl;
	// 	cout << "矩阵ub:\n" << ub << endl;

	///////////////////////////////////////Start the solution/////////////////////////////////////

#if 1	
	//Copy the matrix
	MatrixXd A_cons_temp(A_cons_ori.cols(), A_cons_ori.rows());
	A_cons_temp = -A_cons_ori.transpose();

	QuadProgPP::Matrix<double> H_Mat(H.rows(), H.cols());
	QuadProgPP::Vector<double> f_Mat((f.cols() > f.rows() ? f.cols() : f.rows()));
	QuadProgPP::Matrix<double> A_cons_Mat(A_cons_temp.rows(), A_cons_temp.cols());
	QuadProgPP::Vector<double> B_cons_Mat((B_cons_ori.cols() > B_cons_ori.rows() ? B_cons_ori.cols() : B_cons_ori.rows()));
	EigenToMatrix(H, H_Mat);
	EigenToVector(f, f_Mat);
	EigenToMatrix(A_cons_temp, A_cons_Mat);
	EigenToVector(B_cons_ori, B_cons_Mat);
	QuadProgPP::Matrix<double> CE;
	QuadProgPP::Vector<double> ce0;
	QuadProgPP::Vector<double> x(Nc);
	QuadProgPP::Matrix<double> A_cons_new;
	QuadProgPP::Vector<double> b_cons_new;
	double f_value = QuadProgPP::solve_quadprog(H_Mat, f_Mat, CE, ce0, A_cons_Mat, B_cons_Mat, x);


#endif
	double result;
	result = U[0] + x[0];
	U[0] = result;
	return result;

}


MpcClass::~MpcClass()
{
	//cout << "Calculation complete!!!" << endl;
}

///
///计算克罗内克积
///
MatrixXd MpcClass::Kron(MatrixXd A, MatrixXd B)
{
	MatrixXd C(A.rows()*B.rows(), A.cols()*B.cols());
	for (int i = 0; i < A.rows(); i++)
	{
		for (int j = 0; j < A.cols(); j++)
		{
			C.block(i*B.rows(), j*B.cols(), B.rows(), B.cols()) = A(i, j)*B;
		}
	}
	return C;
}
///
///对角元设置
///
void MpcClass::Diag_Mat(MatrixXd &rst, MatrixXd innerMat, int num)
{
	for (int i = 0; i < num; i++)
	{
		rst.block(i*innerMat.rows(), i*innerMat.cols(), innerMat.rows(), innerMat.cols()) = innerMat;
	}
}
///
///矩阵的幂次方
///
MatrixXd MpcClass::Pow_Mat(MatrixXd x, int y)
{
	MatrixXd rst = MatrixXd::Identity(x.rows(), x.rows());
	if (y != 0)
	{
		for (int i = 1; i <= y; i++)
		{
			rst = rst*x;
		}
	}
	return rst;
}
///
///实现四个矩阵合并
///
void MpcClass::Mer_Mat(MatrixXd &dsttmp, MatrixXd &src11, MatrixXd &src12, MatrixXd &src21, MatrixXd &src22)
{
	MatrixXd dst((src11.rows() + src21.rows()), (src11.cols() + src12.cols()));
	//Mer_Mat合并
	dst.block(0, 0, src11.rows(), src11.cols()) = src11;
	dst.block(0, src11.rows(), src12.rows(), src12.cols()) = src12;
	dst.block(src11.rows(), 0, src21.rows(), src21.cols()) = src21;
	dst.block(src11.rows(), src11.cols(), src22.rows(), src22.cols()) = src22;

	dsttmp = dst;

}
///
///纵向四个矩阵合并
///
void MpcClass::Mer_Vertical(MatrixXd &dsttmp, MatrixXd &src1, MatrixXd &src2, MatrixXd &src3, MatrixXd &src4)
{
	MatrixXd dst((src1.rows() + src2.rows() + src3.rows() + src4.rows()), src1.cols());
	dst.block(0, 0, src1.rows(), src1.cols()) = src1;
	dst.block(src1.rows(), 0, src2.rows(), src2.cols()) = src2;
	dst.block(src1.rows() + src2.rows(), 0, src3.rows(), src3.cols()) = src3;
	dst.block(src1.rows() + src2.rows() + src3.rows(), 0, src4.rows(), src4.cols()) = src4;

	dsttmp = dst;
}

void MpcClass::Mer_Demo(void)
{
	//合成 5 X 4 矩阵
	MatrixXd X11(3, 3);
	X11 << 1, 1, 1, 1, 1, 1, 1, 1, 1;
	MatrixXd X12(3, 2);
	X12 << 2, 2, 2, 2, 2, 2;
	MatrixXd X21(1, 3);
	X21 << 3, 3, 3;
	MatrixXd X22(1, 2);
	X22 << 4, 4;
	MatrixXd Mer_Mat_tmp(1, 1);
	Mer_Mat(Mer_Mat_tmp, X11, X12, X21, X22);

}

#if 1
///
///两个库之间模板的转换
///
void MpcClass::EigenToMatrix(Eigen::MatrixXd &src, QuadProgPP::Matrix<double> &dst)
{
	for (int i = 0; i < src.rows(); i++)
	{
		for (int j = 0; j < src.cols(); j++)
		{
			dst[i][j] = src(i, j);
		}
	}
}
void MpcClass::EigenToVector(Eigen::MatrixXd &src, QuadProgPP::Vector<double> &dst)
{
	if (src.rows() > src.cols())
	{
		for (int i = 0; i < src.rows(); i++)
		{
			dst[i] = src(i, 0);
		}
	}
	else
	{
		for (int i = 0; i < src.cols(); i++)
		{
			dst[i] = src(0, i);
		}
	}
}
#endif