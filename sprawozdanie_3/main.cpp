/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

#include"opt_alg.h"
#include <cstdlib>
#include <ctime>
#include <windows.h>

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();
void recur();

int main()
{	
	//recur();

	try
	{
		lab3();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	return 0;
}

void recur()
{
	try
	{
		lab2();
	}
	catch (string EX_INFO)
	{
		recur();
	}
}

void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;		// liczba mniejsza ni¿ epsilon jest uznawana za zero
	int Nmax = 10000;			// maksymalna liczba iteracji, po której funkcja przerywa dzia³anie



	// wektory o 2 zmiennych		lb, ub
	// [0]   (-5) (5)
	// [1]   (-5) (5)

	matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1); // teraz a ma 0
	solution opt;
	a(0) = -1;		// x1		sta³e okreœlaj¹ce funkcje celu
	a(1) = 2;		// x2
	// to jak rozumiem, s¹ wspó³rzêdne startowe od których licz¹c gradienty (chyba jeszcze nie)
	// one teraz s¹ czym w³aœciwie ? - punktem pocz¹tkowym


	// MC (Monte Carlo) to sposób rozwi¹zywania zadania optymalizacyjnego
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
	// ff0T - funkcja celu (to jej ekstremum globalnego szukamy)
	// - zwraca jak¹œ wartoœæ a my na jej podstawie albo zatrzymujemy obliczenia, bo tutaj wszystko jest
	// skalowane tak, ¿e wartoœæ pomiêdzy 0, a epsilon (bardzo blisko 0) jest uznawana jako rozwi¹zanie
	
	// w 1D, szuka siê miejsca gdzie pochodna jest równa 0, bo to jest ekstremum
	// wiêcej wymiarów to szukamy gdzie gradient jest równy 0
	//             (wektor pokazuj¹cy kierunek najszybszego wzrostu wartoœci funkcji)


	cout << opt << endl << endl;
	solution::clear_calls();

	//Wahadlo
	Nmax = 1000;
	epsilon = 1e-2;
	lb = 0;
	ub = 5;
	double teta_opt = 1;
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1);
	matrix MT = matrix(2, new double[2]{ m2d(opt.x),0.5 });

	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
	ofstream Sout("symulacja_lab0.csv");
	Sout << hcat(Y[0], Y[1]);
	Sout.close();
	Y[0].~matrix();
	Y[1].~matrix();
}

void lab1()
{	
	double epsilon = 1e-2;		// liczba mniejsza ni¿ epsilon jest uznawana za zero
	double gamma = 1e-7;
	int Nmax = 1000;

	solution opt;
	
	// Sprawdzanie Fib (nie dzia³a dla pierwszych parametrów)
	{
		double a, b;
		a = -10; b = 10;
						//x = -2,10655e-12;
						//y = -7,44939e-18;
						//f_calls = 62					
						//Exit flag: 0

				/* k - 1, k				k - 2 , k -1
			x = 5, 67444e-06;		x = -1, 37721e-05;
			y = 6, 36536e-16;		y = 4, 14043e-16;
			f_calls = 56			f_calls = 56
			Exit flag : -1			Exit flag : -1*/



		a = 10; b = 100;			//	te same wyniki

						//x = 62,727;
						//y = -0,921198;
						//f_calls = 68
						//Exit flag: 0

		opt = fib(funkcja_celu_lab_2, a, b, 1e-5);	
	}



	// Tabela nr 1 -> ró¿ne alpha
	{
		double  X0, d, alpha, zmiana;

		alpha = 2.0;
		alpha = 1.75;
		alpha = 1.42;

		epsilon = 0.001;

		X0 = -12; d = 1.0; Nmax = 1000;						// [36, 84]
		zmiana = 0.65;

		for (int i = 0; i < 100; i++)
		{
			double* borders = expansion(funkcja_celu_lab_2, X0, d, alpha, Nmax);
			cout << X0 << " " << borders[0] << " " << borders[1] << " " << solution::f_calls;
			X0 += zmiana;
			solution::clear_calls();


			opt = fib(funkcja_celu_lab_2, borders[0], borders[1], epsilon);
			opt.my_pick();
			solution::clear_calls();


			opt = lag(funkcja_celu_lab_2, borders[0], borders[1], epsilon, gamma, Nmax);
			opt.my_pick();
			solution::clear_calls();

			cout << endl;
			delete[] borders;
		}
	}
	
	// Tabela nr 2 -> œrednie z wczeœniejszej tabeli

	// Wykres -> przejœcie krok po kroku interacji i sprawdzanie b.x - a.x
	{
		cout << "Fibonachi" << endl;
		solution::clear_calls();
		opt = fib(funkcja_celu_lab_2, -100, 100, epsilon);
		opt.my_pick();										cout << endl;
		cout << "Lagrange" << endl;
		solution::clear_calls();
		opt = lag(funkcja_celu_lab_2, -100, 100, epsilon, gamma, Nmax);
		opt.my_pick();
		solution::clear_calls();
	}
	
	// Tabela nr 3
	{
		solution fib_result = fib(ff1R, 0.0001, 0.01, 0.0000001, 0.0000001, 1000);
		//solution fib_result = lag(ff1R, 0.0001, 0.01, 0.0000001, 0.0000001, 1000);
		cout << fib_result.x(0) << endl;
		cout << fib_result.y(0) << endl;
		cout << solution::f_calls << endl;
	}

	// Symulacja
	{
		ofstream file("lab1_tabela_wynikow.csv");
		random_device rd;
		mt19937 gen(rd());
		uniform_real_distribution<double> distribution(-100.0, 100.0);
		if (!file.is_open()) {
			cout << "file error" << endl;
			return;
		}		
			
		//solution fib_result = fib(ff1R, 0.0001, 0.01, 0.0000001, 0.0000001, 1000);
		solution fib_result = lag(ff1R, 0.0001, 0.01, 0.0000001, 0.0000001, 1000);
		matrix Y0 = matrix(3, new double[3]{ 5, 1, 10 });
		matrix* Y = solve_ode(df1, 0, 1, 1000, Y0, NAN, fib_result.x(0));
		file << Y[1] << endl;
		file.close();
	}
}

void lab2()
{
	solution Xopt;
	//matrix x(2, 1, 0.0);
	//x(0) = 2.73503;
	//x(1) = 6.89565;

	// x(0) = 0.623824;
	// x(1) = 0.45588;

	//x(0) = 5.10107;
	//x(1) = 8.62952;


	// x(1) = 1;
	// Xopt = HJ(funkcja_celu_lab_3, x, 0.5, 0.5, 1e-7, 1000);
	
	//Xopt = HJ(funkcja_celu_lab_3, x, 0.5, 0.5, 1e-7, 1000);

	double krok{};

	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_real_distribution<double> dist(-1., 1.);
	matrix copy;

	/**/
	double s = 0.85, alphaHJ = 0.5, alphaR = 1.3, beta = 0.5, epsilon = 1e-3, Nmax = 1000;
	matrix s0(2, 1, 0.5);

	//for (int i = 0; i < 100; i++)
	{
		//double* points = new double[2]{ dist(rng), dist(rng) };

		//double* points = new double[2]{ 2.818104, 3.841883};
		double* points = new double[2] {2.817589, 3.841482};
		
		matrix x(2, points);delete[] points;

		copy = x;
		cout << x(0) << " ";
		cout << x(1);

		Xopt.clear_calls();
		//Xopt = HJ(funkcja_celu_lab_3, x, s, alphaHJ, epsilon, Nmax);
		//Xopt = HJ(ff2R, x, s, alphaHJ, epsilon, Nmax);
		//Xopt.my_pick();		
		
		

		cout << "\n";
		cout << "\n";
		cout << "\n";

		x = copy;
		cout << " ";
		Xopt.clear_calls();
		//Xopt = Rosen(funkcja_celu_lab_3, x, s0, alphaR, beta, epsilon, Nmax);
		//Xopt = Rosen(ff2R, x, s0, alphaR, beta, epsilon, Nmax);
		matrix wyn;
		wyn = ff2R(x);
		//Xopt.my_pick();
		cout << "\n";

	}
	

	/*
	krok = 1.57;
	for (int i = 0; i < 100; i++)
	{
		double* points = new double[2]{ dist(rng), dist(rng) };
		matrix x(2, points);		delete[] points;
		copy = x;
		cout << x(0) << " ";
		cout << x(1);

		Xopt.clear_calls();
		Xopt = HJ(funkcja_celu_lab_3, x, krok, 0.69, 1e-7, 1000);
		Xopt.my_pick();

		x = copy;
		cout << " ";
		Xopt.clear_calls();
		Xopt = Rosen(funkcja_celu_lab_3, x, s, krok, 0.69, 0.001, 1000);
		Xopt.my_pick();
		cout << "\n";
	}
	*/

	/*
	krok = 1.69;
	for (int i = 0; i < 100; i++)
	{
		double* points = new double[2]{ dist(rng), dist(rng) };
		matrix x(2, points);		delete[] points;
		copy = x;
		cout << x(0) << " ";
		cout << x(1);

		Xopt.clear_calls();
		Xopt = HJ(funkcja_celu_lab_3, x, krok, 0.69, 1e-7, 1000);
		Xopt.my_pick();

		x = copy;
		cout << " ";
		Xopt.clear_calls();
		Xopt = Rosen(funkcja_celu_lab_3, x, s, krok, 0.69, 0.001, 1000);
		Xopt.my_pick();
		cout << "\n";
	}
	*/
}

/*
void lab3()
{	
	//// ============= arkusz 1 ===============
	//solution opt_3;
	//double r;
	//srand(time(NULL));
	//int c0 = 2;
	//int dc;
	//int eps = 1e-8;
	//int Nmax = 1000;

	//// zmiana w kazdym z 3 przypadkow
	////matrix a = 4;
	////matrix a = 4.4934;
	//matrix a = 5;


	//for (int i = 0; i < 100; i++) {
	//	double x0_1 = (rand() % 10000) / 10000.0f;
	//	double x0_2 = (rand() % 10000) / 10000.0f;
	//	matrix x0(2, 1, 0.0f);
	//	x0.set_row(x0_1, 0);
	//	x0.set_row(x0_2, 1);

	//	// zewnetrzne
	//	dc = 2;
	//	solution::clear_calls();
	//	std::cout << x0_1 << "\t" << x0_2 << "\t";
	//	opt_3 = pen(fun3, x0, c0, dc, eps, Nmax, a);
	//	r = norm(opt_3.x);
	//	std::cout << opt_3.x(0) << "\t" << opt_3.x(1) << "\t" <<
	//	r << "\t" << opt_3.y << "\t" << solution::f_calls << "\t";

	//	// wewnetrzne
	//	dc = 0.5;
	//	solution::clear_calls();
	//	opt_3 = pen(fun3, x0, c0, dc, eps, Nmax, a);
	//	r = norm(opt_3.x);
	//	std::cout << opt_3.x(0) << "\t" << opt_3.x(1) << "\t" <<
	//	r << "\t" << opt_3.y << "\t" << solution::f_calls << "\n";
	//}	
	// ============ problem rzezczywisty =================
std::cout << pen(fR3, matrix(2, 1, 2), 1, 2, 1e-5, 10000, 4);
}
*/

void lab3()
{
	matrix x0 = matrix(2, 1, 1.0);

	double c0 = 1;
	const double epsilon = 1e-3;
	const int Nmax = 10000;

	matrix a[3] = { 4, 4.4934, 5 };
	const int ktory = 2;
	
	double pkt1[101] = { 0 };
	double pkt2[101] = { 0 };

	for (int i = 0; i < 101; i++) {
		do
			x0 = 5 * rand_mat(2, 1) + 1;
		while (norm(x0) > a[ktory]);

		pkt1[i] = x0(0);
		pkt2[i] = x0(1);
	}
	
	
	//for (int i = 0; i < 100; i++) {
	//	x0(0) = pkt1[i];
	//	x0(1) = pkt2[i];

	//	cout << x0(0) << " " << x0(1) << " ";

	//	solution zew_opt = pen(fun3, x0, c0, 2, epsilon, Nmax, a[ktory]);
	//	cout << zew_opt.x(0) << " " << zew_opt.x(1) << " " << norm(zew_opt.x) << " " << zew_opt.y[0] << solution::f_calls << " ";
	//	solution::clear_calls();

	//	solution wew_opt = pen(fun3, x0, c0, 0.5, epsilon, Nmax, a[ktory]);
	//	cout << wew_opt.x(0) << " " << wew_opt.x(1) << " " << norm(wew_opt.x) << " " << wew_opt.y[0] << solution::f_calls << endl;
	//	solution::clear_calls();
	//}
	

	// Problem rzeczywisty
	
	matrix x1 = matrix(2, 1);
	x1(0) = 0.;
	x1(1) = 0.;

	//x0(0) = pkt1[100];
	//x0(1) = pkt2[100];
	//solution pr_rzeczywisty = pen(fR3, x0, c0, 2, epsilon, Nmax, a[ktory]);
	//cout << pr_rzeczywisty.x(0) << " " << pr_rzeczywisty.x(1) << " " << norm(pr_rzeczywisty.x) << " " << pr_rzeczywisty.y[0] << solution::f_calls << " ";
	


	//cout << pen(fR3, x1, c0, 2, epsilon, Nmax) << endl;

	//-3.3473 25.0006
	
	x1(0) = -3.88; // -10  10
	x1(1) = 9.14; // -23  23

	cout << pen(fR3, x1, c0, 2, epsilon, Nmax) << endl;
	//fR3(x1, c0, 2);
	cout << solution::f_calls << endl;;
}


void lab4()
{
	

}

void lab5()
{

}

void lab6()
{

}
