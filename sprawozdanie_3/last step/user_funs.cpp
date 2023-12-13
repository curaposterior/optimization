#include"user_funs.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include "opt_alg.h"

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
	double* point = new double[2]{};
	matrix dY = matrix(2, point); delete[] point;
	dY(0) = y(1);
	double I_b = m_robot * pow(l_robot, 2) / 3 + m_c * pow(l_robot, 2);
	double M_t = ud2(0) * (ud1(0) - y(0)) + ud2(1) * (ud1(1) - y(1));
	dY(1) = (M_t - b_t * y(1)) / I_b;
	return dY;
}

matrix ff2R(matrix x, matrix ud1, matrix ud2)
{
	double* point = new double[2]{ 0., 0. };
	matrix y0 = matrix(2, point); delete[] point;

	point = new double[2]{ M_PI, 0. };
	matrix y_ref(2, point); delete[] point;

	matrix* Y = solve_ode(df2, 0, 0.1, 100, y0, y_ref, x);
	matrix y = 0.;
	int n = get_len(Y[0]);

	ofstream FILE("output.txt");
	FILE << Y[1];
	FILE.close();

	delete[] Y;

	return y * 0.1;
}

matrix Fun3(matrix x1, matrix x2, matrix ud1) {
	return (sin(M_PI * sqrt(pow(x1() / M_PI, 2) + pow(x2() / M_PI, 2)))) / (M_PI * sqrt(pow(x1() / M_PI, 2) + pow(x2() / M_PI, 2)));
}

bool g1(matrix x1) {
	if (-x1() + 1 <= 0)
		return true;
	else
		return false;
}
bool g2(matrix x2) {
	if (-x2() + 1 <= 0)
		return true;
	else
		return false;
}
bool g3(matrix x1, matrix x2, double alpha) {
	if (sqrt(pow(x1(), 2) + pow(x2(), 2)) - alpha <= 0)
		return true;
	else
		return false;
}

matrix df3(double t, matrix Y, matrix ud1, matrix ud2) {
	double C = 0.47;
	double r = 0.12; // 12 cm
	double m = 0.6; // 600 g
	double ro = 1.2; // 1.2 kg/m3
	double g = 9.81; // 9.81 m/s2

	double S = M_PI * pow(r, 2);

	double Dx = 0.5 * C * ro * S * Y(1) * abs(Y(1));
	double Dy = 0.5 * C * ro * S * Y(3) * abs(Y(3));

	double Fmx = M_PI * ro * Y(3) * m2d(ud2) * pow(r, 3);
	double Fmy = M_PI * ro * Y(1) * m2d(ud2) * pow(r, 3);

	matrix dY(4, 1);
	dY(0) = Y(1);
	dY(1) = (-Dx - Fmx) / m;
	dY(2) = Y(3);
	dY(3) = (-m * g - Dy - Fmy) / m;

	return dY;
}

matrix fR3(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	matrix Y0(4, new double[4] { 0, x(0), 100, 0 });
	matrix* Y = solve_ode(df3, 0, 0.01, 7, Y0, ud1, x(1)); // to ud1?
	int n = get_len(Y[0]);
	int i0 = 0, i50 = 0;

	//cout << x(0) << " " << x(1) << endl;

	ofstream FILE("out.txt");
	double t = 0;
	for (int i = 0; i < n; i++) {
		if (abs(Y[1](i, 2) - 50) < abs(Y[1](i50, 2) - 50))
			i50 = i;
		if (abs(Y[1](i, 2)) < abs(Y[1](i0, 2)))
			i0 = i;

		cout << t << ";";
		FILE << t << ";";
		t += 0.01;
		cout << Y[1](i, 0) << ";";
		FILE << Y[1](i, 0) << ";";
		//cout << Y[1](i, 1) << ";";
		cout << Y[1](i, 2) << ";";
		FILE << Y[1](i, 2) << ";";
		//cout << Y[1](i, 3) << ";";
		cout << endl;
		FILE << endl;
	}
	FILE.close();

	y = -Y[1](i0, 0);

	if (abs(x(0)) - 10 > 0)
		y = y + ud2 * pow(abs(x(0)) - 10, 2);
	if (abs(x(1)) - 25 > 0)
		y = y + ud2 * pow(abs(x(1)) - 25, 2);
	if (abs(Y[1](i50, 0) - 5) - 1 > 0)
		y = y + ud2 * pow(abs(Y[1](i50, 0) - 5) - 1, 2);

	//cout << Y0(0) << " " << Y0(1) << " " << Y0(2) << " " << Y0(3) << endl;
	//cout << x(0) << " " << x(1) << endl; // z optymalnego dla tego uruchomic symulacje

	//cout << "y = " << y << endl;

	return y;
}

matrix fun3(matrix x, matrix ud1, matrix ud2) {
	double arg = M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2));
	matrix y = sin(arg) / arg;
	// y = pow(x(0),2) + pow(x(1),2);

	if (ud2(1) > 1) { //kara zew
		if (-x(0) + 1 > 0)		// g1
			y = y + (ud2)(0) * pow(-x(0) + 1, 2);
		if (-x(1) + 1 > 0)		// g2
			y = y + (ud2)(0) * pow(-x(1) + 1, 2);
		if (norm(x) - (ud1)(0) > 0)  // g3
			y = y + (ud2)(0) * pow(norm(x) - (ud1)(0), 2);
	}
	else { //kara wew
		if (-x(0) + 1 > 0)
			y = 1e10;
		else
			y = y - (ud2)(0) / (-x(0) + 1);

		if (-x(1) + 1 > 0)
			y = 1e10;
		else
			y = y - (ud2)(0) / (-x(1) + 1);

		if (norm(x) - (ud1)(0) > 0)
			y = 1e10;
		else
			y = y - (ud2)(0) / (sqrt(pow(x(0), 2) + pow(x(1), 2)) - (ud1)(0));
	}
	return y;
}



matrix fun4(matrix x, matrix ud1, matrix ud2) {
	return pow(x(0) + 2 * x(1) - 7, 2) + pow(2 * x(0) + x(1) - 5, 2);
}

matrix grad4(matrix x, matrix ud1, matrix ud2) {
	matrix g(2, 1);
	g(0) = 10 * x(0) + 8 * x(1) - 34;
	g(1) = 8 * x(0) + 10 * x(1) - 38;
	return g;
}

matrix hesj4(matrix x, matrix ud1, matrix ud2) {
	matrix H(2, 2);
	H(0, 0) = 10;
	H(0, 1) = 8;
	H(1, 0) = 8;
	H(1, 1) = 10;
	return H;
}

matrix fT4(matrix x, matrix ud1, matrix ud2) {
	matrix y;

	if (isnan(ud2(0, 0)))
		y = pow(x(0) + 2 * x(1) - 7, 2) + pow(2 * x(0) + x(1) - 5, 2);
	else
		y = fT4(ud2[0] + x * ud2[1], 0, ud1);
	return y;
}

matrix fR4(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	int m = 100;
	int n = get_len(x);
	static matrix X(n, m), Y(1, m);
	if (solution::f_calls == 1) {
		ifstream in("XData.txt");
		in >> X;
		in.close();
		in.open("YData.txt");
		in >> Y;
		in.close();
	}
	int P = 0;
	double h;
	y = 0;
	for (int i = 0; i < m; i++) {
		h = m2d(trans(x) * X[i]);
		h = 1.0 / (1.0 + exp(-h));
		y = y - Y(0, i) * log(h) - (1 - Y(0, i)) * log(1 - h);
	}
	y = y / m;
	return y;
}
