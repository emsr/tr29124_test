#include <complex.h>
#include <iostream.h>

int cbroyden(void (*f)(complex<double> *x,complex<double> *fv,int n),complex<double> *x,
    complex<double> *fv,int n,double *eps,int *iter);

void f(complex<double> *x,complex<double> *fv,int n)
{
    fv[0] = cos(2.0*x[0]) - cos(2.0*x[1]) -0.4;
    fv[1] = 2.0*(x[1]-x[0]) + sin(2.0*x[1]) - sin(2.0*x[0]) - 1.2;
}
void f2(complex<double> *x,complex<double> *fv,int n)
{
    fv[0] = x[0] + x[1] -3.0;
    fv[1] = x[0]*x[0] + x[1]*x[1] - 9.0;
}
int main()
{
    complex<double> *x,*r;
    double eps;
    int i,n,ecode,iter;

    n = 2;
    x = new complex<double>[n];
    r = new complex<double>[n];

    x[0] = 0.1;
    x[1] = 0.5;
    eps = 1e-8;
    iter = 10;
    ecode = cbroyden(f,x,r,n,&eps,&iter);
    cout << "Iterations: " << iter << endl;
    cout << "x[0] = " << x[0] << endl;
    cout << "x[1] = " << x[1] << endl;

    delete [] x;
    delete [] r;
    return 0;
}
