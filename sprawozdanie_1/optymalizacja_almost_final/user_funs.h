#pragma once

#include"ode_solver.h"
#include <cmath>


matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);

matrix funkcja_celu_lab_2(matrix , matrix = NAN, matrix = NAN);
matrix f_celu_test(matrix, matrix = NAN, matrix = NAN);