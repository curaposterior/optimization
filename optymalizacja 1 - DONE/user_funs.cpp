#include"user_funs.h"
#define _USE_MATH_DEFINES
#include <math.h>


matrix ff0T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);
	int n = get_len(Y[0]);
	double teta_max = Y[1](0, 0);
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));
	Y[0].~matrix();
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);
	double m = 1, l = 0.5, b = 0.5, g = 9.81;
	double I = m*pow(l, 2);
	dY(0) = Y(1);
	dY(1) = ((t <= ud2(1))*ud2(0) - m*g*l*sin(Y(0)) - b*Y(1)) / I;
	return dY;
}

matrix funkcja_celu_lab_2(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = -cos(0.1 * x(0)) * std::exp(-1.0 * pow((0.1 * x(0) - 2 * 3.14), 2)) + 0.002 * pow(0.1 * x(0), 2);
	
	return y;
}

matrix df1(double t, matrix Y, matrix ud1, matrix ud2) {
	double fOutA = 0, fOutB = 0;
	double fInB = 0.01, tInB = 10.0;
	double tA = 90.0;
	matrix dY(3, 1);
	if (Y(0) > 0)fOutA = 0.98 * 0.63 * m2d(ud2) * sqrt(2.0 * 9.81 * Y(0) / 0.7);
	if (Y(1) > 0)fOutB = 0.98 * 0.63 * 0.00365665 * sqrt(2.0 * 9.81 * Y(1) / 1.0);
	dY(0) = -fOutA;
	dY(1) = (fOutA + fInB - fOutB);
	dY(2) = (fInB / Y(1)) * (tInB - Y(2)) + (fOutA / Y(1)) * (tA - Y(2));
	return dY;
}

matrix ff1R(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	matrix Y0 = matrix(3, new double[3]{ 5, 1, 10 });
	matrix* Y = solve_ode(df1, 0, 1, 1000, Y0, ud1, x);
	int n = get_len(Y[0]);
	double max = Y[1](0, 2);
	for (int i = 0; i < n; i++) {
		if (max < Y[1](i, 2)) max = Y[1](i, 2);
	}
	y = abs(max - 50);
	return y;
}


matrix funkcja_celu_lab_3(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * M_PI * x(0)) - cos(2.5 * M_PI * x(1)) + 2;
	return y;
}