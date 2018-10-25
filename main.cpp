#include <iostream>
#include <fstream>
#include <iomanip> 
#include <cmath>
#include <cfloat>


using namespace std;
const int w = 8;
const double eps = FLT_EPSILON/10;



template <typename T> 
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}



void mulmatr(int n, double **a,double **t, double cosf, double sinf, int i, int j){
	double tempstri[n], tempstrj[n], tempclni[n], tempclnj[n];
	for (int k = 0 ; k < n; ++k) {
		tempstri[k] = cosf * a[i][k] - sinf * a[j][k];
		tempstrj[k] = sinf * a[i][k] + a[j][k] * cosf;
		//cout << "i" << k << ' ' << tempstri[k] << endl;
		//cout << "j" << k << ' ' << tempstrj[k] << endl;
		//использовать tempstri;
	}
	cout << endl << endl;
	for (int k = 0 ; k < n; ++k) {
		tempclni[k] = a[k][i] * cosf - sinf * a[k][j];
		tempclnj[k] = a[k][i] * sinf + cosf * a[k][j];
		//cout << "i" << k << ' ' << tempstri[k] << endl;
		//cout << "j" << k << ' ' << tempstrj[k] << endl;
	}
	tempclni[i] = tempstri[i] * cosf - sinf * tempstri[j];
	tempclni[j] = tempstri[i] * sinf + cosf * tempstri[j];
	tempclnj[i] = tempstrj[i] * cosf - sinf * tempstrj[j];
	tempclnj[j] = tempstrj[i] * sinf + cosf * tempstrj[j];
	//cout << tempclni[i] << ' ' << tempclni[j] <<' ' << tempclnj[i] << ' ' << tempclnj[j] << " " << tempclnj[1] <<endl;
	tempstri[i] = tempclni[i];
	tempstri[j] = tempclnj[i];
	tempstrj[i] = tempclni[j];
	tempstrj[j] = tempclnj[j];

	for(int k = 0 ; k < n; ++k) {
		a[i][k] = tempstri[k];
		a[j][k] = tempstrj[k];
		a[k][i] = tempclni[k];
		a[k][j] = tempclnj[k];
	}
	
}



void freemem(double **&a, int n){
	int i;
	for(i = 0; i < n; ++i) {
		delete[] a[i];
	}
	delete[] a;
}



void printmatr(int n, double **a) {
	int i,j;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			cout << setw(w) << a[i][j] << ' ';
		}
		cout << endl;
	}
	cout << endl;
}



double func(int i,int j, int n) {
	return 	(2 * i + 2 * j) * (i != j) + 3 * (i == j);
}




void selection(int mode, int &i1, int &j1, double ** a, int n) {
	switch (mode) { 
		case 1:
			{

			double max = 0;
			for ( int i = 0; i < n; ++i) {
				for (int j = 0; j < n; ++j) {
					if (abs(a[i][j]) * (i != j) > max) {
						i1 = i;
						j1 = j;
						max = abs(a[i][j]);
					}
				}
			}
			break;
			}
		case 2:
			break;
		case 3:
			if ( j1 == n - 1) {
				if (i1 == n - 1) {
					i1 = 0;
					j1 = 1;
					return;
				} else {
					++i1;
					j1 = i1 + 1;
					return ;
				}
			} else {
				++j1;
				return;
			}
			break;
		default:
			cout << "error mode";
			exit(1);
			break;
	}
}



void answer(int n, double **a, double **t, int mode) {
	double end = 0;
	for( int i = 0 ; i < n; ++i) {
		for( int j = 0 ; j < n; ++j) {
			end += abs(a[i][j]) * (i != j);
		}
	}	
	double sumstr[n];
	for(int i = 0 ; i < n; ++i) {
		sumstr[i] = 0;
		for (int j = 0 ; j < n ; ++j) {
			sumstr[i] += a[i][j] * (i != j);
		}
	}

	while ( end > eps) {
		int i = 0;
		int j = 0;
		selection(mode, i, j, a, n);
		double x = - 2 * a[i][j];
		double y = a[i][i] - a[j][j];
		double cosf, sinf;
		if ( y == 0 ) {
			cosf = 1.0 / sqrt(2);
			sinf = cosf;
		} else {
			cosf = sqrt(1.0/2 + abs(y)/(2 * sqrt(x * x + y * y)));
			sinf = sgn( x * y) * abs(x) / (2 * cosf * sqrt(x*x + y * y));
		}
		//cout << "cos sin" << cosf << ' ' << sinf << endl;
		mulmatr(n, a, t, cosf, sinf,i, j);
		/*if (abs(a[i][j]) < eps) {
			cout << "a :" << endl;
			printmatr(n,a);
			cout << "t :" << endl;
			printmatr(n,t);
			return ;
		}*/
		cin >> x;
		printmatr(n,a);
		//cout << i << ' ' << j << endl;
	}
	
}



int main(){
	int n;
	double **a;

	cout << "1 - fromfile, 2 - generate" << endl;
	int genmatr;
	cin >> genmatr;
	cout << "mode :" << endl;
	int mode;
	cin >> mode;
	if (genmatr == 1) {
		ifstream fin("input.txt");
		if (!fin.is_open()) {
			cout << "cant open input.txt" << endl;
			return 1;
		}
		fin >> n;
		a = new double *[n];
		int i,j;
		for(i = 0; i < n; ++i) {
			a[i] = new double[n];
		}
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				fin >> a[i][j];
			}
		}
		//printmatr(a,n);
	} else {
		cout << "size = ";
		cin >> n;
		a = new double *[n];
		int i,j;
		for(i = 0; i < n; ++i) {
			a[i] = new double[n];
		}
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				a[i][j] = func(i,j,n);
			}
		}
		printmatr(n,a);
	}
	double **t = new double *[n];
	for (int i = 0; i < n; ++i) {
		t[i] = new double [n];
	}
	for (int i = 0 ; i < n; ++i) {
		for (int j = 0 ; j < n; ++ j) {
			t[i][j] = (i == j);
		}
	}
	printmatr(n,t);
	answer(n, a, t, mode);
	freemem(a,n);
	freemem(t,n);
	return 0;
}