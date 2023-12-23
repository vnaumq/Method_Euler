

#include <iostream>
#include <math.h>
#include <iomanip>
#include <time.h>
#include <stdlib.h>
#include <fstream>
#include <conio.h>
using namespace std;

const double Eps = 0.001;
const int itr_max = 100;
double func(double* u, double t, int n)
{
	if (n == 0)
	{
		if (t < Eps)
			return -u[0]*u[1] + 1;
		else return -u[0]*u[1] + (sin(t)/t);
	}
	else if (n == 1)
		return -u[1]*u[1] + (3.5*t)/(1 + t*t);
}

void EulerJava(double* u, int n)
{
	double Tau = 0;
	double T = 1, TauMax = 0.01;
	double tk = 0;
	double* yk = new double[n];
	//Epsilon zadano sverhu

	for (int i = 0; i < n; i++) 
		yk[i] = u[i];

	int vybor;
	cout << "1 - txt, 2 - console" << endl;
	cin >> vybor;

  switch (vybor)
  {
	case 1:
	{
	ofstream out;
	out.open("EulerJavny.txt", ios::trunc);

	for (int i = 0; i < n; i++)
	{
		if (i == 0)
		{
			out << "    t" << setw(15);
			out << "    u" << i + 1 << setw(14);
		}
		else if (i == n - 1)
			out << "    u" << i + 1 << endl;
	}

	int kol = 0;
	do
	{
		double *tmp = new double[n];
		for (int i = 0; i < n; i++)
			tmp[i] = func(yk, tk, i);

		if (Eps/(fabs(tmp[0]) + Eps/TauMax) > Eps/(fabs(tmp[1]) + Eps/TauMax))
			Tau = Eps / (fabs(tmp[1]) + Eps / TauMax);
		else
			Tau = Eps / (fabs(tmp[0]) + Eps / TauMax);  //formula 3.11 i 3.12 (minimalnoje Tau)

		for (int i = 0; i < n; i++)  // funkcija na novom shage
			yk[i] += Tau * tmp[i];

		tk += Tau;

		for (int i = 0; i < n; i++) // otvet
		{
			if (i == 0)
				out << tk << setw(15);
			out << setw(15) << yk[i];
			if (i == n - 1)
				out << endl;
		}
	
		kol++;
	} while (tk < T);

	out << endl << "Iterations quantity is " << kol << endl;
	cout << endl << "Done" << endl;
	out.close();
	break;
	}
	
	case 2:
	{
	for (int i = 0; i < n; i++)
	{
		if (i == 0)
		{
			cout << "    t" << setw(15);
			cout << "    u" << i + 1 << setw(14);
		}
		else if (i == n - 1)
			cout << "    u" << i + 1 << endl;
	}

	int kol = 0;
	do
	{
		double *tmp = new double[n];
		for (int i = 0; i < n; i++)
			tmp[i] = func(yk, tk, i);

		if (Eps/(fabs(tmp[0]) + Eps/TauMax) > Eps/(fabs(tmp[1]) + Eps/TauMax))
			Tau = Eps / (fabs(tmp[0]) + Eps / TauMax);
		else
			Tau = Eps / (fabs(tmp[1]) + Eps / TauMax);

		for (int i = 0; i < n; i++)  // funkcija na novom shage
			yk[i] += Tau * tmp[i];

		tk += Tau;

		for (int i = 0; i < n; i++) // otvet
		{
			if (i == 0)
				cout << tk << setw(15);
			cout << setw(15) << yk[i];
			if (i == n - 1)
				cout << endl;
		}
	
		kol++;
	} while (tk < T);

	cout << endl << "Iterations quantity is " << kol << endl;
	break;
	}
  }
}



double** createMatrix(int x)
{
	double **A = new double*[x];
	for (int i = 0; i < x; i++)
		A[i] = new double[x + 1];
	return A;
}
//освобождение памяти
void deleteMatrix(double **X, int x)
{
	for (int i = 0; i < x; i++)
		delete X[i];
	delete[] X;
}
//копирования для сохранения оригинала
void copyMatrix(double **X, double **copyX, int x)
{
	for (int i = 0; i < x; i++)
	{
		for (int j = 0; j < x + 1; j++)
			copyX[i][j] = X[i][j];
	}
}

double f1(double *uk1, double *uk, double t, double Tau)
{
	return uk1[0] - uk[0] - Tau*(-uk1[0] * uk1[1] + ((t < 1e-9) ? 0.0 : (sin(t) / t)));
}

double f2(double *uk1, double *uk, double t, double Tau)
{
	return uk1[1] - uk[1] - Tau*(-uk1[1] * uk1[1] + (3.125*t) / (1 + t*t));
}
//весь метод Гаусса
bool Gauss(double *An, double **X, int x)
{
	for (int k = 0; k < x; k++)
	{
		double max = fabs(X[k][k]);
		int remeber = k;		//запоминаем строку, чтобы не поменять саму себя
		for (int i = k + 1; i < x; i++)
		{
			if (max < fabs(X[i][k]))		//находим максимальный по модулю элемент в столбце
			{
				max = fabs(X[i][k]);
				remeber = i;
			}
		}

		if (fabs(max - 0) < 1e-6)
		{
			return 0;
		}

		if (k != remeber)				//меняем строки местами
		{
			double *temp = X[k];
			X[k] = X[remeber];
			X[remeber] = temp;
		}

		//viewMatrix(X, x);

		double lead = X[k][k];			//запоминаем ведущий элемент
		for (int r = k; r < x + 1; r++)
		{
			X[k][r] /= lead;
		}
		//начиная со следующей строки приводим исходную матрицу к диагональному виду
		for (int i = k + 1; i < x; i++)
		{
			double temp = X[i][k];
			for (int j = k; j < x + 1; j++)
			{
				X[i][j] -= X[k][j] * temp;
			}
		}
		//viewMatrix(X, x);
	}

	An[x - 1] = X[x - 1][x + 1 - 1];				//обратный ход
	for (int i = x - 2; i >= 0; i--)
	{
		An[i] = X[i][x + 1 - 1];
		for (int j = i + 1; j < x + 1 - 1; j++)
		{
			An[i] -= X[i][j] * An[j];
		}
	}
	return 1;
}

typedef double(*pf)(double*, double*, double, double);
//находим производную по n аргументу
double Differential(pf f, double *uk1, double *uk, double t, double Tau, int n)
{
	double dx = 1e-9;
	double* D = new double[2];
	for (int i = 0; i < 2; i++)
	{
		if (i == n - 1)
		{
			D[i] = uk1[i] + dx;
			i++;
		}
		D[i] = uk1[i];
	}

	double F = f(uk1, uk, t, Tau);
	double dF = f(D, uk, t, Tau);

	delete[] D;
	return (dF - F) / dx;
}
//метод Ньютона
double* neutone(double *yk_plus, double *yk, double tk, double Tau, int n)
{
	double **Jako = createMatrix(n);
	double *An = new double[n];
	double e = 1e-9, b1 = 0, b2 = 0;

	int itr = 1;

	do{
		//заполнение матрицы Якоби в ячейки F[0][0], F[0][1], F[1][0], F[1][1]
		Jako[0][0] = Differential(f1, yk_plus, yk, tk, Tau, 1);
		Jako[0][1] = Differential(f1, yk_plus, yk, tk, Tau, 2);
		Jako[0][2] = -f1(yk_plus, yk, tk, Tau);					//вектор столбец с -
		Jako[1][0] = Differential(f2, yk_plus, yk, tk, Tau, 1);
		Jako[1][1] = Differential(f2, yk_plus, yk, tk, Tau, 2);
		Jako[1][2] = -f2(yk_plus, yk, tk, Tau);					//вектор столбец с -

		double **copyF = createMatrix(n);					// копия для сохранения оригинала
		copyMatrix(Jako, copyF, n);

		if (!Gauss(An, copyF, n))				//Проверка метода Гаусса
		{
			deleteMatrix(copyF, n);		//очистка памяти
			deleteMatrix(Jako, n);
			delete An;
			return 0;
		}

		yk_plus[0] += An[0];				//уточнение решения
		yk_plus[1] += An[1];

		if (fabs(f1(yk_plus, yk, tk, Tau)) > fabs(f2(yk_plus, yk, tk, Tau))) //подсчет первой погрешности
			b1 = fabs(f1(yk_plus, yk, tk, Tau));
		else
			b1 = fabs(f2(yk_plus, yk, tk, Tau));

		for (int i = 0; i < n; i++)										//подсчет второй погрешности
		{
			if (fabs(An[i]) < 1)
				b2 = fabs(An[i]);
			else if (fabs(An[i]) >= 1)
				b2 = fabs(An[i] / yk_plus[i]);
		}

		deleteMatrix(copyF, n);
		itr++;
	} while ((b1 > e || b2 > e) && (itr < itr_max));

	deleteMatrix(Jako, n);				//очистка памяти
	delete An;

	return yk_plus;
}



void EulerNejava(double* u, int n)
{
	double Tau, Tau_minus, Tau_plus;
	double T = 1, TauMax = 0.01, TauMin = 0.01;
	double tk = 0, tk_plus;
	double Eps_k;  //Epsilon_k
	double* yk = new double[n];               
	double* yk_minus = new double[n];              
	double* yk_plus = new double[n];               
	Tau_minus = Tau = TauMin;
	//Epsilon zadano sverhu

	for (int i = 0; i < n; i++) 
		yk[i] = yk_minus[i] = yk_plus[i] = u[i];

	ofstream out;
	out.open("EulerNejavny.txt");

	int sposob;
	cout << "1 - kvazioptimum, 2 - trjohzon" << endl;
	cin >> sposob;

	for (int i = 0; i < n; i++)
	{
		if (i == 0)
		{
			out << "    t" << setw(13);
			out << "    u" << i + 1 << setw(13);
		}
		else if (i == n - 1)
			out << "    u" << i + 1 << endl;
	}

	int kol = 0;
	do
	{
		do
		{
			tk_plus = tk + Tau;
			yk_plus = neutone(yk_plus, yk, tk, Tau, n);

			for (int k = 0; k < n; k++)  // Eps_k po formule  (3.16)
				Eps_k = -(Tau/(Tau + Tau_minus))*(yk_plus[k] - yk[k] - Tau*(yk[k] - yk_minus[k]) / Tau_minus);

			for (int k = 0; k < n; k++)  //если  Eps_k > Eps, то выполняем и переходим в начало
			if (fabs(Eps_k) > Eps)
			{
				Tau /= 2;
				tk_plus = tk;
				for (int j = 0; j < n; j++)
					yk_plus[j] = yk[j];
			}
		}
		while (fabs(Eps_k) > Eps);
		
		if (sposob == 1)          // kvazioptimum
			Tau_plus = sqrt(Eps/abs(Eps_k))*Tau;

		else if (sposob == 2)     // trjohzon
		{
			for (int i = 0; i < n; i++)  
			{
				if (fabs(Eps_k) > Eps)
					Tau_plus = Tau/2;
				if ((Eps/4 < fabs(Eps_k)) && (fabs(Eps_k) <= Eps))
					Tau_plus = Tau;
				if (fabs(Eps_k) <= Eps / 4)
					Tau_plus = 2*Tau;
			}
		}
		else break;

		if (Tau_plus > TauMax) 
			Tau_plus = TauMax;

		for (int i = 0; i < n; i++) 
		{
			if (i == 0)
				out << tk << setw(15);
			out << setw(15) << yk[i];
			if (i == n - 1)
				out << endl;
		}

		for (int i = 0; i < n; i++)  // sdvig
		{
			yk_minus[i] = yk[i];
			yk[i] = yk_plus[i];
		}
		Tau_minus = Tau;
		Tau = Tau_plus;
		tk = tk_plus;

		kol++;
	} while(tk < T);

	cout << endl << "Iterations quantity is " << kol << endl;
	out << endl << "Iterations quantity is " << kol << endl;
	out.close();
}


int main()
{
	int n = 2;
	double* u = new double[n];

	cout << "Start value of U[0] = 0.0 U[1] = -0.412 " << endl;
	u[0] = 0.0;
	u[1] = -0.412;

	int method;
	cout << "1 - javnyj, 2 - nejavnyj" << endl;
	cin >> method;
	switch (method)
	{
		case 1:
		{
			cout << endl << "Jav Method" << endl;
			EulerJava(u, n);
			break;
		}

		case 2:
		{
			cout << endl << "Nejav Method" << endl;
			EulerNejava(u, n);
			break;
		}
	}
}
 