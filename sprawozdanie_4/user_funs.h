#pragma once

#include"ode_solver.h"
#include <cmath>


matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);

matrix funkcja_celu_lab_2(matrix , matrix = NAN, matrix = NAN);

matrix df1(double t, matrix Y, matrix = NAN, matrix = NAN);
matrix ff1R(matrix x, matrix = NAN, matrix = NAN);

matrix funkcja_celu_lab_3(matrix, matrix = NAN, matrix = NAN);

matrix df2(double t, matrix y, matrix ud1 = NAN, matrix ud2 = NAN);

matrix ff2R(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);

matrix Fun3(matrix x1, matrix x2, matrix ud1);

bool g1(matrix x1);

bool g2(matrix x2);

bool g3(matrix x1, matrix x2, double alpha);

matrix df3(double t, matrix Y, matrix ud1, matrix ud2);

matrix fR3(matrix x, matrix ud1, matrix ud2);

matrix fun3(matrix x, matrix ud1, matrix ud2);

matrix fun4(matrix x, matrix ud1, matrix ud2);
matrix grad4(matrix x, matrix ud1, matrix ud2);
matrix hesj4(matrix x, matrix ud1, matrix ud2);

matrix fT4(matrix x, matrix ud1, matrix ud2);

matrix fR4(matrix x, matrix ud1, matrix ud2);
matrix gf(matrix x, matrix ud1, matrix ud2);

matrix fT5(matrix x, matrix ud1, matrix ud2);
matrix fR5(matrix x, matrix ud1, matrix ud2);

matrix f5(matrix x, matrix ud1, matrix ud2);
matrix f5_1(double a, matrix x, matrix ud1, matrix ud2);
matrix f5_2(double a, matrix x, matrix ud1, matrix ud2);

matrix fT6(matrix x, matrix ud1, matrix ud2);