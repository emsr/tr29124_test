

cp ../tr1_mathN.diff.bz2 .
svn checkout http://gcc.gnu.org/svn/gcc/trunk/libstdc++-v3/include
svn checkout http://gcc.gnu.org/svn/gcc/trunk/libstdc++-v3/testsuite
svn checkout http://gcc.gnu.org/svn/gcc/trunk/libstdc++-v3/docs
bzcat tr1_mathN.diff.bz2 | patch -p0

