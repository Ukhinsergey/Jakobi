#include <iostream>
#include <fstream>
#include <iomanip> 
#include <cmath>
#include <cfloat>


using namespace std;
const int w = 8;
const double eps = FLT_EPSILON/10000000;



void freemem(double **&a, int n){
	int i;
	for(i = 0; i < n; ++i) {
		delete[] a[i];
	}
	delete[] a;
}



void printmatr(double **a, int n) {
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
	return 	2 * i + j;
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

	while ( end > eps) {
		int i, j;
		selection(mode, i, j, a, n);
		cout << i << ' ' << j << endl;
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
		printmatr(a,n);
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
	printmatr(t, n);
	answer(n, a, t, mode);
	freemem(a,n);
	freemem(t,n);
	return 0;
}