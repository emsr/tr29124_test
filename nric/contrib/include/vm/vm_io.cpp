#ifndef VM_IO_CPP
#define VM_IO_CPP

/******************************************************************************
 *
 *	vm_io.cpp
 *	
 *	Source file containing streaming IO operators for Vectors and Matrices.
 *
 *	This source file is included by the header file vm_io.h, and it is
 *	dependent on the following header files:
 *
 *	vec_mat.h  vm_traits.h  vector.h  matrix.h  vm_io.h
 *
 ******************************************************************************
 *
 *	VecMat Software, version 1.05
 *	By Kevin Dolan (2002)
 *	kdolan@mailaps.org
 *
 *****************************************************************************/

// Stream operators for vectors.
template <class T>
inline std::ostream& operator<<(std::ostream& out, const Vector<T>& x)
{
	x.write(out);
	return out;
}

template <class T>
std::istream& operator>>(std::istream& in, Vector<T>& x)
{
	T tmp;
	size_t i = 0, n = x.size();
	while (in >> tmp && i < n)
	{
		x[i] = tmp;
		++i;
	}
	return in;
}


// Stream operators for matrices.
template <class T>
inline std::ostream& operator<<(std::ostream& out, const Matrix<T>& x)
{
	x.write(out);
	return out;
}

template <class T>
std::istream& operator>>(std::istream& in, Matrix<T>& x)
{
	T tmp;
	size_t i = 0, j = 0, n = x.nrows(), m = x.ncols();
	while (in >> tmp && i < n)
	{
		x(i, j) = tmp;
		++j;
		if (j == m)
		{
			j = 0;
			++i;
		}
	}
	return in;
}


// Binary file functions for vectors. 
template <class T>
void ReadBinary(const std::string filename, Vector<T>& array, size_t skip)
{
	size_t        i, len;
	struct stat	  file_status;
	std::ifstream file;

	i = stat(filename.c_str(), &file_status);

	if (i == -1)
	{
		array.free();
		return;
	}

	if ((file_status.st_size - skip) % sizeof(T))
	{
		array.free();
		return;
	}

	len = (file_status.st_size - skip) / sizeof(T);

	if (!len)
	{
		array.free();
		return;
	}

	array.reshape(len);

	file.open(filename.c_str(), std::ios::in | std::ios::binary);

	file.seekg(skip, std::ios::beg);
	file.read(reinterpret_cast<char *>(array.data()), len * sizeof(T));

	if (file.eof())
	{
		array.free();
		return;
	}

	file.close();

	return;
}

template <class T>
void WriteBinary(const std::string filename, const Vector<T>& array)
{
	std::ofstream file;

	file.open(filename.c_str(), std::ios::out | std::ios::binary);

	size_t    len = array.size();
	ptrdiff_t   s = array.stride();
	const T    *d = array.data();

	if (s == 1) file.write(reinterpret_cast<const char *>(d), len * sizeof(T));
	else
	{
		for (size_t i = 0; i < len; ++i)
		{
			file.write(reinterpret_cast<const char *>(d), sizeof(T));
			d += s;
		}
	}

	file.close();
}


// Binary file functions for matrices. 
template <class T>
void ReadBinary(const std::string filename, Matrix<T>& array, size_t skip)
{
	size_t        i, len;
	struct stat	  file_status;
	std::ifstream file;

	i = stat(filename.c_str(), &file_status);

	if (i == -1)
	{
		array.free();
		return;
	}

	if ((file_status.st_size - skip) % sizeof(T))
	{
		array.free();
		return;
	}

	len = (file_status.st_size - skip) / sizeof(T);

	if (!len)
	{
		array.free();
		return;
	}

	array.reshape(1, len);

	file.open(filename.c_str(), std::ios::in | std::ios::binary);

	file.seekg(skip, std::ios::beg);
	file.read(reinterpret_cast<char *>(array.data()), len * sizeof(T));

	if (file.eof())
	{
		array.free();
		return;
	}

	file.close();

	return;
}

template <class T>
void WriteBinary(const std::string filename, const Matrix<T>& array)
{
	std::ofstream file;

	file.open(filename.c_str(), std::ios::out | std::ios::binary);

	size_t    nr = array.nrows(), nc = array.ncols();
	ptrdiff_t rs = array.rstride(), cs = array.cstride();
	const T   *d = array.data();
	
	if (rs == nc && cs == 1) file.write(reinterpret_cast<const char *>(d), nr * nc * sizeof(T));
	else
	{
		size_t i, j;
		for (i = 0; i < nr; ++i)
		{
			d = array.data() + i * rs;
			for (j = 0; j < nc; ++j)
			{
				file.write(reinterpret_cast<const char *>(d), sizeof(T));
				d += cs;
			}
		}
	}

	file.close();
}


#endif // VM_IO_CPP
