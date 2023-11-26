/*********************************************
Kod stanowi uzupe�nienie materia��w do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udost�pniony na licencji CC BY-SA 3.0
Autor: dr in�. �ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G�rniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

#include"opt_alg.h"

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		lab1();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	system("pause");
	return 0;
}

void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;		// liczba mniejsza ni� epsilon jest uznawana za zero
	int Nmax = 10000;			// maksymalna liczba iteracji, po kt�rej funkcja przerywa dzia�anie



	// wektory o 2 zmiennych		lb, ub
	// [0]   (-5) (5)
	// [1]   (-5) (5)

	matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1); // teraz a ma 0
	solution opt;
	a(0) = -1;		// x1		sta�e okre�laj�ce funkcje celu
	a(1) = 2;		// x2
	// to jak rozumiem, s� wsp�rz�dne startowe od kt�rych licz�c gradienty (chyba jeszcze nie)
	// one teraz s� czym w�a�ciwie ? - punktem pocz�tkowym


	// MC (Monte Carlo) to spos�b rozwi�zywania zadania optymalizacyjnego
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
	// ff0T - funkcja celu (to jej ekstremum globalnego szukamy)
	// - zwraca jak�� warto�� a my na jej podstawie albo zatrzymujemy obliczenia, bo tutaj wszystko jest
	// skalowane tak, �e warto�� pomi�dzy 0, a epsilon (bardzo blisko 0) jest uznawana jako rozwi�zanie
	
	// w 1D, szuka si� miejsca gdzie pochodna jest r�wna 0, bo to jest ekstremum
	// wi�cej wymiar�w to szukamy gdzie gradient jest r�wny 0
	//             (wektor pokazuj�cy kierunek najszybszego wzrostu warto�ci funkcji)


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
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(opt.x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
	ofstream Sout("symulacja_lab0.csv");
	Sout << hcat(Y[0], Y[1]);
	Sout.close();
	Y[0].~matrix();
	Y[1].~matrix();
}

void lab1()
{	
	double epsilon = 1e-2;		// liczba mniejsza ni� epsilon jest uznawana za zero
	int Nmax = 10000;

	solution opt;


	double* borders = expansion(funkcja_celu_lab_2, -50, 1.0, 1.5, Nmax);
	cout << "borders (" << borders[0] << ", " << borders[1] << ")\n";


	//opt = fib(funkcja_celu_lab_2, borders[0], borders[1], epsilon);
	//opt = fib(funkcja_celu_lab_2, 60, 65, epsilon);
	
	//cout << opt << endl << endl;
	solution::clear_calls();
}

void lab2()
{

}

void lab3()
{

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
