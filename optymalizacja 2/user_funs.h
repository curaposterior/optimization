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