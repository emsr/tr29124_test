
     double McLachlan(char typ,int n, double h);

     double Estimatmcn(char typ,int n, double h);
     double Estimatmcnl(char typ,int n, double t);
     double nChain(char typ, int n, double t);


     double MCNRoot(double a, char typ, int n, double h, double accr, double &acc);
     double NewtonRaphsons(double a, char typ, int n, double h, double &acp);

     double Asympn(char typ, int n, double h);
     double Mcapn(int n, double h);
     double Approxn(char typ, int n, double  h);

     double Coefficients(char typ,int n, double h, double Bm[], int CDIM, double &mcnr, double &Norm);
     int CoefficientSec(char typ, int n, double h, double mcn, double Bm[], int CDIM, double Gm[], double &Norm);
