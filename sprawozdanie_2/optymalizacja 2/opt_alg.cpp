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
		double* p = new double[2] { 0, 0 };
		solution x0(_x0);
		solution x1(_x0 + d);

		std::vector<solution> x_val;
		x_val.push_back(x0.x);
		x_val.push_back(x1.x);

		x_val[0].fit_fun(ff, ud1, ud2);
		x_val[1].fit_fun(ff, ud1, ud2);

		if (x_val[1].y(0) == x_val[0].y(0))
		{
			p[0] = x_val[0].x(0); p[1] = x_val[1].x(0);
			return p;
		}
		if (x_val[1].y(0) > x_val[0].y(0))
		{
			d = -d;
			x_val[1].x = _x0 + d;

			x_val[1].fit_fun(ff, ud1, ud2);
			if (x_val[1].y(0) >= x_val[0].y(0))
			{
				p[0] = x_val[1].x(0); p[1] = x_val[0].x(0) - d;
				return p;
			}
		}

		int i = 0;

		while(true)
		{
			if (solution::f_calls > Nmax)
			{
				delete[] p;
				return NULL;
			}

			i++;

			x_val.push_back(_x0 + pow(alpha, i) * d);
			x_val[i + 1].fit_fun(ff, ud1, ud2);
			if (x_val[i].y(0) <= x_val[i + 1].y(0))
			{
				break;
			}
		}
		

		if (d > 0)
		{
			p[0] = x_val[i - 1].x(0); p[1] = x_val[i + 1].x(0);
			return p;
		}
		else
		{
			p[0] = x_val[i + 1].x(0); p[1] = x_val[i - 1].x(0);
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

	return res;
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double _a, double _b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		int k{};
		vector<double> fib_tab = fibonachi(80);
		for (int i = 0; i < fib_tab.size(); i++) if (fib_tab[i] > ((_b - _a) / epsilon)) { k = i; break; }

		solution a(_a);
		solution b(_b);
		solution c;
		solution d;

		c.x = b.x - (fib_tab[k - 2] / fib_tab[k - 1]) * (b.x - a.x);		

		d.x = (a.x + b.x - c.x);
		for (int i = 0; i <= k - 3; ++i)
		{
			c.fit_fun(ff, ud1, ud2);
			d.fit_fun(ff, ud1, ud2);
			
			if (c.y < d.y)
			{
				b = d;
			}
			else
			{	
				a = c;
			}

			c.x = b.x - (fib_tab[k - i - 2] / fib_tab[k - i - 1]) * (b.x - a.x);			
			d.x = a.x + b.x - c.x;

			//cout << b.x - a.x << endl;
		}		

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
		double l;
		double m;

		solution a(_a);	
		solution b(_b);
		solution c((_a + _b) / 2);
		solution d_i_minus(c);
		solution d(0);


		a.fit_fun(ff);
		b.fit_fun(ff);
		c.fit_fun(ff);

		int i = 0;	bool after_first = false;
		while(true)
		{
			l = m2d((a.y * ((pow(b.x, 2) - pow(c.x, 2)))) 
					+ (b.y * ((pow(c.x, 2) - pow(a.x, 2)))) 
					+ (c.y * ((pow(a.x, 2) - pow(b.x, 2)))));

			m = m2d((a.y * (b.x - c.x)) 
					+ (b.y * (c.x - a.x)) 
					+ (c.y * (a.x - b.x)));


			if (m <= 0) ERROR;
			

			d.x = 0.5 * (l / m);
			d.fit_fun(ff);
			if (a.x(0) < d.x(0) && d.x(0) < c.x(0))
			{
				if (d.y(0) < c.y(0))
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
				if (c.x(0) < d.x(0) && d.x(0) < b.x(0))
				{
					if (d.y(0) < c.y(0))
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

			if (((b.x(0) - a.x(0)) < epsilon) || (abs(d.x(0) - d_i_minus.x(0)) < gamma))
			{
				d_i_minus = d;
				d_i_minus.flag = 1;
				break;
			}

			d_i_minus = d;

			//cout << b.x - a.x << endl;
		}

		return d_i_minus;		 
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}
//linia 9 w pseudokodzie, powinno byc PROBUJ(x, s)

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	solution Xopt;

	try
	{
		solution xb = x0;
		solution x = x0;
		solution xb_;

		while (true) 
		{
			xb = x;
			x = HJ_trial(ff, xb, s);
			x.fit_fun(ff);
			xb.fit_fun(ff);

			if (x.y < xb.y)
			{
				while (true) 
				{
					xb_ = xb;
					xb = x;
					x = 2 * xb.x - xb_.x;
					x = HJ_trial(ff, xb, s);
					if (solution::f_calls > Nmax) ERROR;
					
					x.fit_fun(ff);
					xb.fit_fun(ff);
					if (x.y < xb.y) break;
				}
				x = xb;
			}
			else
			{
				s = alpha * s;
			}
			if (solution::f_calls > Nmax) ERROR;
			if (s >= epsilon) break;
		}
		Xopt = xb;
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
		matrix e1(2, 1);	e1(0) = 1; e1(1) = 0;
		matrix e2(2, 1);	e2(0) = 0; e2(1) = 1;

		std::vector<matrix> e = { e1, e2 };
		XB.fit_fun(ff);

		for (int j = 0; j < 2; j++)
		{
			solution pom1(XB.x + s * e[j]);
			solution pom2(XB.x - s * e[j]);

			pom1.fit_fun(ff);
			pom2.fit_fun(ff);

			if (pom1.y < XB.y)
			{
				XB.x = XB.x + s * e[j];
			}
			else
			{
				if (pom2.y < XB.y)
				{
					XB.x = XB.x - s * e[j];
				}
			}
		}
		return XB;		
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

double max(matrix* tab, int length)
{
	double biggest = tab[0](0); if (biggest < 0) biggest = -biggest;

	double incoming;
	for (int i = 1; i < length; i++)
	{
		incoming = tab[i](0); if (incoming < 0) incoming = -incoming;

		if (biggest < incoming) biggest = incoming;
	}
	return biggest;
}

// vj(i + 1) = Q*j(i)

// jak rozumiem (i) jest naszą wcześniejszą iteracją czyli nie powinniśmy inkrementować "i" przed tymi obliczeniami
//

matrix dj_khaled(int j)
{
	matrix upper = return_vj(j);
	double lower = norm(upper);

	return upper / lower;
}

matrix return_vj(int j)	// j może byc tylko 1 lub 2 // (jako indexy 0 i 1)
{
	matrix res;
	if (j == 1)
	{
		// res = Q(i)1 	tutaj pytanie, czy to jest "j"ota kolumna czy "j"ty wiersz ???
	}
	else if (j == 2)
	{
		// res = Q(i)2 - ((trans(Q(i)2 * d[0])) * d[0])				d[0] -> to jest to co obliczyliśmy we wcześniejszej iteracji
	}

	return res;
}

// Matrix to (wiersze, kolumny)

matrix funkcja_co_wysiaga_z_matrixa_2x2___1_wiersz_albo_1_kolumne(matrix _2x2_, int j) //zmienic nazwe
{
	//// kolumna
	//{
	//	matrix res(2, 1); // 2 wiersze i 1 kolumna

	//	res(0, 0) = _2x2_(0, 0);
	//	res(1, 0) = _2x2_(1, 0);

	//	return res;
	//}

	//// wiersz
	//{
	//	matrix res(1, 2); // 1 wiersz i 2 kolumny

	//	res(0, 0) = _2x2_(0, 0);
	//	res(0, 1) = _2x2_(0, 1);

	//	return res;
	//}
}

//				wiersz  kolumna
//				(0, 0)	(0, 1)
//				(1, 0)	(1, 1)

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	solution Xopt;

	try
	{
		const int n = 2;
		int i = 0;

		matrix e1(2, 1);	e1(0) = 1; e1(1) = 0;
		matrix e2(2, 1);	e2(0) = 0; e2(1) = 1;

		matrix* d = new matrix[n];		d[0] = e1;			d[1] = e2;
		matrix* lambda = new matrix[n];	lambda[0] = 0;		lambda[1] = 0;
		matrix* p = new matrix[n];		p[0] = 0;			p[1] = 0;
		matrix* s = new matrix[n];		s[0] = s0;			s[1] = s0;
		solution xb = x0;		


		solution tmp1, tmp2;
		while (true)
		{
			for (int j = 0; j < n; j++)
			{
				tmp1 = xb.x + (s[j] * d[j]);				
				
				tmp1.fit_fun(ff);
				xb.fit_fun(ff);
				if (tmp1.y < xb.y)
				{
					xb = tmp1;
					lambda[i + 1] = lambda[i] + s[i];
					s[i + 1] = alpha * s[i];
				}
				else
				{
					s[i + 1] = (-beta) * s[i];
					p[i + 1] = p[i] + 1;
				}
			}

			i++;
			Xopt = xb;

			bool value = true;
			for (int _i_ = 0; _i_ < n; _i_++)
			{
				if (!(lambda[_i_] != 0 && p[_i_] != 0))
				{
					value = false;
					break;
				}
			}

			// zmiana bazy kierunków
			if (value)
			{
				matrix D(2, 2);
				
				D(0, 0) = d[0](0);
				D(1, 0) = d[0](1);
				
				D(0, 1) = d[1](0);
				D(1, 1) = d[1](1);


				matrix tmp_lambda_matrix(2, 2);

				tmp_lambda_matrix(0, 0) = lambda[0](0);
				tmp_lambda_matrix(1, 0) = lambda[1](0);
				
				tmp_lambda_matrix(0, 1) = 0;
				tmp_lambda_matrix(1, 1) = lambda[1](0);

				matrix Q = D * tmp_lambda_matrix;



				for (int _i_ = 0; _i_ < n; _i_++)
				{
					lambda[_i_](0) = 0;
					p[_i_](0) = 0;
					s[_i_](0) = s[0](0);
				}
			}

			// duza petla
			if (solution::f_calls > Nmax) ERROR;
			if (!(max(s, n) < epsilon)) break;
		}

		delete[] d;
		delete[] lambda;
		delete[] p;

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