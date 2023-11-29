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

	double* point = new double[3]{ 5, 1, 10 };
	matrix Y0 = matrix(3, point);	delete[] point;
	matrix* Y = solve_ode(df1, 0, 1, 1000, Y0, ud1, x);
	int n = get_len(Y[0]);
	double max = Y[1](0, 2);
	for (int i = 0; i < n; i++) {
		if (max < Y[1](i, 2)) max = Y[1](i, 2);
	}
	y = abs(max - 50);

	delete[] Y;
	return y;
}



matrix funkcja_celu_lab_3(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * M_PI * x(0)) - cos(2.5 * M_PI * x(1)) + 2;
	return y;
}

static const double l_robot = 0.6;
static const double m_robot = 1;
static const double m_c = 9.5;
static const double b_t = 0.5;

// y - [alfa(t), omega(t)], ud1 - [alfa_ref, omega_ref], ud2 (x) - [k1, k2]
matrix df2(double t, matrix y, matrix ud1, matrix ud2)
{
	ved point(0, 0);
	matrix dY = matrix(2, point.data());
	dY(0) = y(1);
	double I_b = m_robot * pow(l_robot, 2) / 3 + m_c * pow(l_robot, 2);
	double M_t = ud2(0) * (ud1(0) - y(0)) + ud2(1) * (ud1(1) - y(1));
	dY(1) = (M_t - b_t * y(1)) / I_b;
	return dY;
}

matrix ff2R(matrix x, matrix ud1, matrix ud2)
{
	ved point(0.0, 0.0);
	matrix y0 = matrix(2, point.data());

	point[0] = M_PI; point[1] = 0.0;
	matrix y_ref(2, point.data());

	matrix* Y = solve_ode(df2, 0, 0.1, 100, y0, y_ref, x);
	matrix y = 0.;
	int n = get_len(Y[0]);

	ofstream FILE("output.txt");
	for (int i = 0; i < n; i++)
	{
		double alfa = y_ref(0) - Y[1](i, 0);
		double omega = y_ref(1) - Y[1](i, 1);

		FILE << alfa << " " << omega << endl;

		double M_t = x(0) * alfa + x(1) * omega;
		y = y + 10 * pow(alfa, 2) + pow(omega, 2) + pow(M_t, 2);
	}
	FILE.close();

	delete[] Y;

	return y * 0.1;
}

matrix funkcja_celu_lab_4(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = ((std::sin(M_PI * sqrt(((pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2)))))) / (sqrt(((pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2))))));
	return y;
}