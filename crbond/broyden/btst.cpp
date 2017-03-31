#include <math.h>
#include <iostream.h>

int broyden(void (*f)(double *x,double *fv,int n),double *x,
    double *fv,int n,double *eps,int *iter);

void f(double *x,double *fv,int n)
{
    fv[0] = cos(2.0*x[0]) - cos(2.0*x[1]) -0.4;
    fv[1] = 2.0*(x[1]-x[0]) + sin(2.0*x[1]) - sin(2.0*x[0]) - 1.2;
}
void f2(double *x,double *fv,int n)
{
    fv[0] = x[0] + x[1] -3.0;
    fv[1] = x[0]*x[0] + x[1]*x[1] - 9.0;
}
int main()
 {
    double *x,*r,eps;
    int i,n,ecode,iter;

    n = 2;
    x = new double[n];
    r = new double[n];

    x[0] = 0.1;
    x[1] = 0.5;
    eps = 1e-12;
    iter = 20;
    ecode = broyden(f,x,r,n,&eps,&iter);
    cout << "Iterations: " << iter << endl;
    cout << "x[0] = " << x[0] << endl;
    cout << "x[1] = " << x[1] << endl;

    delete [] x;
    delete [] r;
    return 0;
}
