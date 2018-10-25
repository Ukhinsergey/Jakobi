#include <iostream>
#include <fstream>
#include <iomanip> 
#include <cmath>
#include <cfloat>


using namespace std;
const int w = 8;
const double eps = FLT_EPSILON;



template <typename T> 
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}



void mulmatr(int n, double **a,double **t, double cosf, double sinf, int i, int j, double &end, int mode, double *sumstr){
	double tempstri[n], tempstrj[n], tempclni[n], tempclnj[n];
	for (int k = 0 ; k < n; ++k) {
		tempstri[k] = cosf * a[i][k] - sinf * a[j][k];
		tempstrj[k] = sinf * a[i][k] + a[j][k] * cosf;
		//cout << "i" << k << ' ' << tempstri[k] << endl;
		//cout << "j" << k << ' ' << tempstrj[k] << endl;
	}
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

	//пересчитываем массив сумм строк
	if (mode == 2) {
		sumstr[i] = 0;
		sumstr[j] = 0;
		for (int k = 0 ; k < n; ++k) {
			sumstr[i] += tempstri[k] * tempstri[k] * (k != i);
			sumstr[j] += tempstrj[k] * tempstrj[k] * (k != j) ;
			if ( k!= i && k != j) {
				sumstr[k] -= a[k][i] * a[k][i] ;
				sumstr[k] -= a[k][j] * a[k][j] ;
				sumstr[k] += tempclni[k] * tempclni[k] ;
				sumstr[k] += tempclnj[k] * tempclnj[k] ;
			}
			
		}

	}

	//считаем end
	for(int k = 0 ; k < n; ++k) {
		end -= abs(a[i][k]) * ( i != k);
		end -= abs(a[j][k]) * ( j != k);
		end -= abs(a[k][i]) * ( i != k ) * (k != j);
		end -= abs(a[k][j]) * ( j != k ) * (k != i);
		end += abs(tempstri[k]) * ( i != k);
		end += abs(tempstrj[k]) * ( j != k);
		end += abs(tempclni[k]) * ( i != k ) * (k != j);;
		end += abs(tempclnj[k]) * ( j != k ) * (k != i);
	}


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




void selection(int mode, int &i1, int &j1, double ** a, int n, double *sumstr) {
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
			{
				double maxsum = 0;
				int maxind = 0;
				for(int k = 0; k < n; ++k) {
					if (sumstr[k] > maxsum) {
						maxsum = sumstr[k];
						maxind = k;
					}
				}	
				i1 = maxind;
				double max = 0;
				for (int k = 0 ; k < n; ++k) {
					if ( abs(a[i1][k]) * (k != i1) > max) {
						j1 = k;
						max = abs(a[i1][k]);
				}
			}
			break;
			}
		case 3:
			if ( j1 == n - 1) {
				if (i1 == n - 2) {
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
	long long iter = 0;
	for( int i = 0 ; i < n; ++i) {
		for( int j = 0 ; j < n; ++j) {
			end += abs(a[i][j]) * (i != j);
		}
	}	
	double sumstr[n];
	for(int i = 0 ; i < n; ++i) {
		sumstr[i] = 0;
		for (int j = 0 ; j < n ; ++j) {
			sumstr[i] += a[i][j] * a[i][j] * (i != j);
		}
	}


	int i = 0;
	int j = 0;
	while ( end > eps) {
		++iter;
		selection(mode, i, j, a, n, sumstr);
		//printmatr(n,a);
		//for (int k = 0 ; k < n; ++k) {
		//	cout<< sumstr[k] << endl;
		//}
		//cout << endl << i << j << endl;
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

		mulmatr(n, a, t, cosf, sinf,i, j, end, mode, sumstr);
		//cin >> x;
		
	}
	cout << "number of iterations = " << iter << endl << "lambda:" << endl;
	for (int k = 0 ; k < n; ++k) {
		cout << a[k][k] << ' ';
	}
	
}



int main(){
	int n;
	double **a;
	clock_t time;
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
		if ( n <= 10) {
			printmatr(n,a);
		}
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
	time = clock();
	answer(n, a, t, mode);
	time = - clock();
	cout << endl << setprecision(6) << 1000*((double) -time)/CLOCKS_PER_SEC << "ms" <<endl;
	freemem(a,n);
	freemem(t,n);
	return 0;
}