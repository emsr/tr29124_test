#ifndef VM_TOOLS_H
#define VM_TOOLS_H

/******************************************************************************
 *
 *	vm_tools.h
 *	
 *	Header file for Vector and Matrix tools.
 *
 *	This header file is intended to be included by vec_mat.h, and it depends
 *	on the following other header files:
 *	
 *	vec_mat.h  vm_traits.h  vector.h  matrix.h
 *
 ******************************************************************************
 *
 *	VecMat Software, version 1.05
 *	By Kevin Dolan (2002)
 *	kdolan@mailaps.org
 *
 *****************************************************************************/

// Statistics.
template <class T>
void Stat(const Vector<T>& x, double& m, double& var);

class Histogram
{
protected:
	double low_, high_;
	size_t num_bins_;
	bool   use_range_, use_bins_;

public:
	explicit Histogram(double l = 0, double h = 0, size_t nb = 0) : low_(l), high_(h),
		num_bins_(nb), use_range_(false), use_bins_(false)
	{
		use_range_ = (low_ < high_);
		use_bins_  = (num_bins_ > 0);
	}
	Histogram(const Histogram& hist) : low_(hist.low_), high_(hist.high_), num_bins_(hist.num_bins_),
		use_range_(hist.use_range_), use_bins_(hist.use_bins_) {;}
	~Histogram() {;}

	Histogram& operator=(const Histogram& hist)
	{
		low_       = hist.low_;
		high_      = hist.low_;
		num_bins_  = hist.num_bins_;
		use_range_ = hist.use_range_;
		use_bins_  = hist.use_bins_;
		return *this;
	}
	void set_range(double l, double h)
	{
		if (low_ < high_)
		{
			low_  = l;
			high_ = h;
			use_range_ = true;
		}
		else
		{
			low_  = 0;
			high_ = 0;
			use_range_ = false;
		}
	}
	void set_bins(size_t nb)
	{
		if (nb > 0)
		{
			num_bins_ = nb;
			use_bins_ = true;
		}
		else
		{
			num_bins_ = 0;
			use_bins_ = false;
		}
	}
	double low() const {return low_;}
	double high() const {return high_;}
	size_t bins() const {return num_bins_;}

	template <class T>
	Vector<int> make_hist(const Vector<T>& x) const
	{
		size_t i, j, n, len = x.size();
		double l, h;
		if (!use_range_)
		{
			l = x.min();
			h = x.max();
		}
		else
		{
			l = low_;
			h = high_;
		}
		if (!use_bins_) 
		{
			size_t m = 0;
			for (i = 0; i < len; ++i)
			{
				if (x[i] <= h && l <= x[i]) m++;
			}
			if (m > 1) n = static_cast<size_t>(exp(0.626 + 0.4 * log(m - 1.0)));
			else n = 0;
		}
		else n = num_bins_;
		Vector<int> hist(0, n);

		double delta = (h - l) / n;
		for (i = 0; i < len; ++i)
		{
			j = static_cast<size_t>((x[i] - l) / delta);
			if (j >= 0 && j < n)
				++hist[j];
			if (x[i] == h) ++hist[n-1];
		}
		return hist;
	}
};


// Sorting algorithms.
template <class T>
class IndirectComp
{
protected:
	Vector<T>::const_iterator data_;

public:
	IndirectComp() : data_() {;}
	IndirectComp(const Vector<T>& x) : data_(x.begin()) {;}
	IndirectComp(const IndirectComp& comp) : data_(comp.data_) {;}
	~IndirectComp() {;}

	IndirectComp& operator=(const IndirectComp& comp)
	{
		data_ = comp.data_;
		return *this;
	}
	bool operator()(int a, int b) const
	{
		return (data_[a] < data_[b]);
	}
	void attach(const Vector<T>& x)
	{
		data_ = x.begin_();
	}
};

template <class T>
Vector<T> Sort(const Vector<T>& x);
template <class T>
Vector<int> Index(const Vector<T>& x);

// Matrix multiplication.
template <class T>
Matrix<T> MatMult(const Matrix<T>& x, const Matrix<T>& y);
template <class T>
Vector<T> MatMult(const Vector<T>& x, const Matrix<T>& y);
template <class T>
Vector<T> MatMult(const Matrix<T>& x, const Vector<T>& y);
template <class T>
Matrix<T> MatMult(const Vector<T>& x, const Vector<T>& y);
template <class T>
T DotProduct(const Vector<T>& x, const Vector<T>& y);


#include "vm/vm_tools.cpp"


#endif // VM_TOOLS_H

