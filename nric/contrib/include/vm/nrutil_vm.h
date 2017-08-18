#ifndef _NR_UTIL_H_
#define _NR_UTIL_H_

// This is a modified version of nrutil_nr.h for use with the VecMat software.

// This allows the NR recipes to use the VecMat classes.
#define NRVec Vector
#define NRMat Matrix
#include "vec_mat.h"
using namespace std;

typedef double DP;

template<class T>
inline const T SQR(const T a) {return a*a;}

template<class T>
inline const T MAX(const T &a, const T &b)
        {return b > a ? (b) : (a);}

inline float MAX(const double &a, const float &b)
        {return b > a ? (b) : float(a);}

inline float MAX(const float &a, const double &b)
        {return b > a ? float(b) : (a);}

template<class T>
inline const T MIN(const T &a, const T &b)
        {return b < a ? (b) : (a);}

inline float MIN(const double &a, const float &b)
        {return b < a ? (b) : float(a);}

inline float MIN(const float &a, const double &b)
        {return b < a ? float(b) : (a);}

template<class T>
inline const T SIGN(const T &a, const T &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const float &a, const double &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const double &a, const float &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

template<class T>
inline void SWAP(T &a, T &b)
	{T dum=a; a=b; b=dum;}

namespace NR {
	inline void nrerror(const string& error_text)
	// Numerical Recipes standard error handler
	{
		cerr << "Numerical Recipes run-time error..." << endl;
		cerr << error_text << endl;
		cerr << "...now exiting to system..." << endl;
		exit(1);
	}
}


template <class T>
class NRMat3d {
private:
	int nn;
	int mm;
	int kk;
	T ***v;
public:
	NRMat3d();
	NRMat3d(int n, int m, int k);
	inline T** operator[](const int i);	//subscripting: pointer to row i
	inline const T* const * operator[](const int i) const;
	inline int dim1() const;
	inline int dim2() const;
	inline int dim3() const;
	~NRMat3d();
};

template <class T>
NRMat3d<T>::NRMat3d(): nn(0), mm(0), kk(0), v(0) {}

template <class T>
NRMat3d<T>::NRMat3d(int n, int m, int k) : nn(n), mm(m), kk(k), v(new T**[n])
{
	int i,j;
	v[0] = new T*[n*m];
	v[0][0] = new T[n*m*k];
	for(j=1; j<m; j++)
		v[0][j] = v[0][j-1] + k;
	for(i=1; i<n; i++) {
		v[i] = v[i-1] + m;
		v[i][0] = v[i-1][0] + m*k;
		for(j=1; j<m; j++)
			v[i][j] = v[i][j-1] + k;
	}
}

template <class T>
inline T** NRMat3d<T>::operator[](const int i) //subscripting: pointer to row i
{
	return v[i];
}

template <class T>
inline const T* const * NRMat3d<T>::operator[](const int i) const
{
	return v[i];
}

template <class T>
inline int NRMat3d<T>::dim1() const
{
	return nn;
}

template <class T>
inline int NRMat3d<T>::dim2() const
{
	return mm;
}

template <class T>
inline int NRMat3d<T>::dim3() const
{
	return kk;
}

template <class T>
NRMat3d<T>::~NRMat3d()
{
	if (v != 0) {
		delete[] (v[0][0]);
		delete[] (v[0]);
		delete[] (v);
	}
}

//The next 3 classes are used in artihmetic coding, Huffman coding, and
//wavelet transforms respectively. This is as good a place as any to put them!

class arithcode {
private:
	NRVec<unsigned long> *ilob_p,*iupb_p,*ncumfq_p;
public:
	NRVec<unsigned long> &ilob,&iupb,&ncumfq;
	unsigned long jdif,nc,minint,nch,ncum,nrad;
	arithcode(unsigned long n1, unsigned long n2, unsigned long n3)
		: ilob_p(new NRVec<unsigned long>(n1)),
		iupb_p(new NRVec<unsigned long>(n2)),
		ncumfq_p(new NRVec<unsigned long>(n3)),
		ilob(*ilob_p),iupb(*iupb_p),ncumfq(*ncumfq_p) {}
	~arithcode() {
		if (ilob_p != 0) delete ilob_p;
		if (iupb_p != 0) delete iupb_p;
		if (ncumfq_p != 0) delete ncumfq_p;
	}
};

class huffcode {
private:
	NRVec<unsigned long> *icod_p,*ncod_p,*left_p,*right_p;
public:
	NRVec<unsigned long> &icod,&ncod,&left,&right;
	int nch,nodemax;
	huffcode(unsigned long n1, unsigned long n2, unsigned long n3,
		unsigned long n4) :
		icod_p(new NRVec<unsigned long>(n1)),
		ncod_p(new NRVec<unsigned long>(n2)),
		left_p(new NRVec<unsigned long>(n3)),
		right_p(new NRVec<unsigned long>(n4)),
		icod(*icod_p),ncod(*ncod_p),left(*left_p),right(*right_p) {}
	~huffcode() {
		if (icod_p != 0) delete icod_p;
		if (ncod_p != 0) delete ncod_p;
		if (left_p != 0) delete left_p;
		if (right_p != 0) delete right_p;
	}
};

class wavefilt {
private:
	NRVec<DP> *cc_p,*cr_p;
public:
	int ncof,ioff,joff;
	NRVec<DP> &cc,&cr;
	wavefilt() : cc(*cc_p),cr(*cr_p) {}
	wavefilt(const DP *a, const int n) :  //initialize to array
		cc_p(new NRVec<DP>(n)),cr_p(new NRVec<DP>(n)),
		ncof(n),ioff(-(n >> 1)),joff(-(n >> 1)),cc(*cc_p),cr(*cr_p) {
			int i;
			for (i=0; i<n; i++)
				cc[i] = *a++;
			DP sig = -1.0;
			for (i=0; i<n; i++) {
				cr[n-1-i]=sig*cc[i];
				sig = -sig;
			}
	}
	~wavefilt() {
		if (cc_p != 0) delete cc_p;
		if (cr_p != 0) delete cr_p;
	}
};

//Overloaded complex operations to handle mixed float and double
//This takes care of e.g. 1.0/z, z complex<float>

inline const complex<float> operator+(const double &a,
	const complex<float> &b) { return float(a)+b; }

inline const complex<float> operator+(const complex<float> &a,
	const double &b) { return a+float(b); }

inline const complex<float> operator-(const double &a,
	const complex<float> &b) { return float(a)-b; }

inline const complex<float> operator-(const complex<float> &a,
	const double &b) { return a-float(b); }

inline const complex<float> operator*(const double &a,
	const complex<float> &b) { return float(a)*b; }

inline const complex<float> operator*(const complex<float> &a,
	const double &b) { return a*float(b); }

inline const complex<float> operator/(const double &a,
	const complex<float> &b) { return float(a)/b; }

inline const complex<float> operator/(const complex<float> &a,
	const double &b) { return a/float(b); }

//some compilers choke on pow(float,double) in single precision. also atan2

inline float pow (float x, double y) {return pow(double(x),y);}
inline float pow (double x, float y) {return pow(x,double(y));}
inline float atan2 (float x, double y) {return atan2(double(x),y);}
inline float atan2 (double x, float y) {return atan2(x,double(y));}

#endif /* _NR_UTIL_H_ */
