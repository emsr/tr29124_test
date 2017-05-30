
   double MathuZn(char typ,int n, double h, double u, short kind);
   double MathuZDn(char typ,int n, double h, double u, short kind);

   double MathuSn(char typ,int n, double h, double v, short kind);
   double MathuSDn(char typ,int n, double h, double v, short kind);

   double  MathuJn(char typ, int n, double u, double Bm[], double Jm[], int CDIM);
   double  MathuJDn(char typ, int n, double h, double u, double Bm[], double Jm[], double JDm[], int CDIM);
   double  MathuYn(char typ, int n, double Bm[], double Jms[], double Ym[], int CDIM);
   double  MathuYDn(char typ, int n, double h, double u, double Bm[], double Jms[], double JDms[], double Ym[], double YDm[], int CDIM);

   double  MathuZn(char typ, int n, double Bm[], double Jms[], double Zm[], int CDIM);
   double  MathuZDn(char typ, int n, double h, double u, double Bm[], double Jms[], double JDms[], double Zm[], double ZDm[], int CDIM);

   double  MathuSn(char typ, int n, double v, double Bm[], int CDIM);
   double  MathuSDn(char typ, int n, double v, double Bm[], int CDIM);
   double  ISn(char typ,int ni, double vi, double ti, double ui, double Bm[], int CDIM);
   double  MathuFn(char typ, int n, double v, double Gm[], double Norm, double Bm[], int CDIM);
   double  MathuFDn(char typ, int n, double v, double Gm[], double Norm, double Bm[], int CDIM);

   double MmathuSn(char typ,int n, double h, double v);
   double  MmathuSn(char typ, int n, double v, double Bem[], double Bom[], int CDIM);

