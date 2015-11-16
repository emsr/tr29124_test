tool="$HOME/bin_specfun/bin/g++ -I/home/ed/gcc_specfun/libstdc++-v3/testsuite/util -D__STDCPP_WANT_MATH_SPEC_FUNCS__"
srcs=`find /home/ed/gcc_specfun/libstdc++-v3/testsuite/special_functions -name \*.cc`
srcs="$srcs `find /home/ed/gcc_specfun/libstdc++-v3/testsuite/ext/special_functions -name \*.cc`"
$tool $srcs 2> test_std.txt
