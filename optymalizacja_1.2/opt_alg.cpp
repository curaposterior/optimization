#include"opt_alg.h"
#include <vector>

using std::vector;

#define ERROR { Xopt.flag = 0; break; }
#define var(x) cout << #x << " = " << x << '\n';

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);
			for (int i = 0; i < N; ++i)	Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);	// przeskalowanie do 1 (chyba)

			Xopt.fit_fun(ff, ud1, ud2);		// fit_fun	-> in solution.cpp
											// ff		-> in user_funs.cpp
											// (x - x1)^2 + (x - x2)^2
			if (Xopt.y < epsilon)
			{
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = 0;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(*ff)(matrix, matrix, matrix), double _x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2]{ 0,0 };	// przedzial

		solution x1(_x0 + d);
		solution x0(_x0);

		x1.fit_fun(ff, ud1, ud2);
		x0.fit_fun(ff, ud1, ud2);

		if (x1.y == x0.y)
		{
			p[0] = m2d(x0.x); p[1] = m2d(x1.x);
			return p;
		}
		if (x1.y > x0.y)
		{
			d = -d;
			x1.x = x0.x + d;

			x1.fit_fun(ff, ud1, ud2);
			x0.fit_fun(ff, ud1, ud2);
			if (x1.y >= x0.y)
			{
				p[0] = m2d(x1.x); p[1] = m2d(x0.x - d);
				return p;
			}
		}


		int i = 0;
		solution x_i_plus;
		solution x_i = x1;
		solution x_i_minus = x0;
		do
		{
			if (solution::f_calls > Nmax)
			{
				delete[] p;
				return NULL;
			}


			i++;	// pierwsze okrazenie i == 1
			if (i != 1) // przesuwamy punkty o 1 do przodu
			{
				x_i_minus = x_i;
				x_i = x_i_plus;
			}			
			x_i_plus.x = x0.x + ( pow(alpha, i) * d );

			var(i);
			var(x_i_minus.x);
			var(x_i.x);
			var(x_i_plus.x);

			x_i_plus.fit_fun(ff, ud1, ud2);	// pierwsze okrazenie - x_plus_1 == x_2
			x_i.fit_fun(ff, ud1, ud2);	// pierwsze okrazenie - x_i == x_1

		} while ( x_i.y <= x_i_plus.y );

		if (d > 0)
		{
			p[0] = m2d(x_i_minus.x); p[1] = m2d(x_i_plus.x);
			return p;
		}
		else
		{
			p[0] = m2d(x_i_plus.x); p[1] = m2d(x_i_minus.x);
			return p;
		}
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

vector<double> fibonachi(int limit)
{
	vector<double> res(2, 1);
	res[0] = 1;
	res[1] = 2;

	for (long long i = 1; i < limit; i++)
		res.push_back(res[i] + res[i - 1]);

	//for (auto& v : res) cout << v << endl;

	return res;
}

#define var(x) ;

solution fib(matrix(*ff)(matrix, matrix, matrix), double _a, double _b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;

		int k{};
		vector<double> fib_tab = fibonachi(80);
		for (int i = 0; i < fib_tab.size(); i++) if (fib_tab[i] > ((_b - _a) / epsilon)) { k = i; break; }
		cout << "k = " << k << "\n";

		
		solution a (_a);
		solution b (_b);
		solution c;
		solution d;


		var(a.x);
		var(b.x);
		c.x = b.x - (fib_tab[k - 2] / fib_tab[k - 1]) * (b.x - a.x);


		d.x = (a.x + b.x - c.x);

		var(c.x);
		var(d.x);

		//cout << "\n\n\n";
							// 13
		for (int i = 0; i <= k - 3; ++i)
		{
			c.fit_fun(ff, ud1, ud2);
			d.fit_fun(ff, ud1, ud2);

			var((c.y < d.y));
			if (c.y < d.y)
			{
				//cout << "b - changed\n";
				b = d;
			}
			else
			{
				//cout << "a - changed\n";
				a = c;
			}


			int adder = 0;
			//					  16   0 -> 13
			
			//(fib_tab[k - i - 2] / fib_tab[k - i - 1]
			c.x = b.x - (fib_tab[k - i- 2] / fib_tab[k - i- 1]) * (b.x - a.x);
			d.x = a.x + b.x - c.x;

			var(i);
			var(a.x);
			var(b.x);
			var(c.x);
			var(d.x);
			cout << '\n';
		}
		cout << "DONE" << "\n";

		Xopt = c;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double _a, double _b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		solution l;
		solution m;
		
		solution a(_a);	solution b(_b); 
		solution c(1); solution d(1);
		solution d_i_minus(1);

		solution condition_resolution;
		int i = 0;	bool after_first = false;
		do
		{
			a.fit_fun(ff);
			b.fit_fun(ff);
			c.fit_fun(ff);

			l.x = (a.y * ((pow(b.x, 2) - pow(c.x, 2)))) + (b.y * ((pow(c.x, 2) - pow(a.x, 2)))) + (c.y * ((pow(a.x, 2) - pow(b.x, 2))));
			m.x = (a.y * (b.x - c.x)) + (b.y * (c.x - a.x)) + (c.y * (a.x - b.x));
			
			if (m.x <= 0) ERROR;

			if (after_first) d_i_minus = d;		// sprawdziæ
			d.x = 0.5 * (l.x / m.x);
			
			if (a.x < d.x < c.x)
			{
				d.fit_fun(ff);
				c.fit_fun(ff);

				if (d.y < c.y)
				{
					b = c;
					c = d;
				}
				else
				{
					a = d;
				}
			}
			else
			{
				if (c.x < d.x < b.x)
				{
					d.fit_fun(ff);
					c.fit_fun(ff);
					if (d.y < c.y)
					{
						a = c;
						c = d;
					}
					else
					{
						b = d;
					}
				}
				else ERROR;
			}

			i++;
			if (solution::f_calls > Nmax) ERROR;
			
			condition_resolution = d.x - d_i_minus.x; if (condition_resolution.x < 0) condition_resolution.x = -condition_resolution.x;
		} while ( ((b.x - a.x) < epsilon) || (condition_resolution.x < gamma) );

		Xopt = d;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}