
   double MathuMZn(char typ,int n, double h, double u, short kind);
   double MathuMZDn(char typ,int n, double h, double u, short kind);

   double MathuQn(char typ,int n, double h, double v, short kind);
   double MathuQDn(char typ,int n, double h, double v, short kind);

   double  MathuIn(char typ, int n, double Am[], double Ims[], double Zm[], int CDIM);
   double  MathuKn(char typ, int n, double Am[], double Kms[], double Zm[], int CDIM);
   double  MathuIDn(char typ, int n, double h, double u, double Am[], double Ims[], double IDms[], double Zm[], double ZDm[], int CDIM);
   double  MathuKDn(char typ, int n, double h, double u, double Am[], double Ims[], double IDms[], double Zm[], double ZDm[], int CDIM);

   double  MathuQn(char typ, int n, double v, double Am[], int CDIM);
   double  MathuQDn(char typ, int n, double v, double Am[], int CDIM);
   double  IQn(char typ,int ni, double vi, double ti, double ui, double Am[], int CDIM);
   double  MathuEn(char typ, int n, double v, double Cm[], double Norm, double Am[], int CDIM);
   double  MathuEDn(char typ, int n, double v, double Cm[], double Norm, double Am[], int CDIM);


   double MCoefficients(char typ,int n, double h, double Am[], int CDIM, double *mcnr, double *Norm);
   int MCoefficientSec(char typ, int n, double h, double mcn, double Am[], int CDIM, double Cm[], double *Norm);
