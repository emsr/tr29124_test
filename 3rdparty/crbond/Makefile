
INPUT_DIR = test_input
OUTPUT_DIR = test_output
BIN_DIR = .

all: test_rpoly test_bairstow test_broyden test_broyden_cmplx test_nlnewt test_cnlnewt

test_rpoly: rpoly.cpp
	g++ -g -o test_rpoly rpoly.cpp

test_bairstow: bairstow.c
	gcc -g -o test_bairstow bairstow.c

test_broyden: broyden/broyden.cpp broyden/btst.cpp broyden/ged.cpp
	g++ -g -o test_broyden broyden/broyden.cpp broyden/btst.cpp broyden/ged.cpp

test_broyden_cmplx: broyden_cmplx/cbroyden.cpp broyden_cmplx/cbtst.cpp broyden_cmplx/cged.cpp
	g++ -g -o test_broyden_cmplx broyden_cmplx/cbroyden.cpp broyden_cmplx/cbtst.cpp broyden_cmplx/cged.cpp

test_nlnewt: newton/newt.h newton/ged.cpp newton/nlnewt.cpp newton/test_nlnewt.cpp
	g++ -g -o test_nlnewt newton/ged.cpp newton/nlnewt.cpp newton/test_nlnewt.cpp

test_cnlnewt: newton_cmplx/newt.h newton_cmplx/cged.cpp newton_cmplx/cnlnewt.cpp newton_cmplx/test_cnlnewt.cpp
	g++ -g -o test_cnlnewt newton_cmplx/cged.cpp newton_cmplx/cnlnewt.cpp newton_cmplx/test_cnlnewt.cpp


test: $(OUTPUT_DIR)
	touch $(OUTPUT_DIR)/test_bairstow.prev && cp -f $(OUTPUT_DIR)/test_bairstow.prev $(OUTPUT_DIR)/test_bairstow.prev2
	touch $(OUTPUT_DIR)/test_bairstow.txt && cp -f $(OUTPUT_DIR)/test_bairstow.txt $(OUTPUT_DIR)/test_bairstow.prev
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_bairstow < $(INPUT_DIR)/test_bairstow.in > $(OUTPUT_DIR)/test_bairstow.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_bairstow < $(INPUT_DIR)/test_solver1.in >> $(OUTPUT_DIR)/test_bairstow.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_bairstow < $(INPUT_DIR)/test_solver2.in >> $(OUTPUT_DIR)/test_bairstow.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_bairstow < $(INPUT_DIR)/test_solver3.in >> $(OUTPUT_DIR)/test_bairstow.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_bairstow < $(INPUT_DIR)/test_solver4.in >> $(OUTPUT_DIR)/test_bairstow.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_bairstow < $(INPUT_DIR)/test_solver5.in >> $(OUTPUT_DIR)/test_bairstow.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_bairstow < $(INPUT_DIR)/test_solver6.in >> $(OUTPUT_DIR)/test_bairstow.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_bairstow < $(INPUT_DIR)/test_solver7.in >> $(OUTPUT_DIR)/test_bairstow.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_bairstow < $(INPUT_DIR)/test_solver8.in >> $(OUTPUT_DIR)/test_bairstow.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_rpoly < $(INPUT_DIR)/test_solver1.in > $(OUTPUT_DIR)/test_rpoly.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_rpoly < $(INPUT_DIR)/test_solver2.in >> $(OUTPUT_DIR)/test_rpoly.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_rpoly < $(INPUT_DIR)/test_solver3.in >> $(OUTPUT_DIR)/test_rpoly.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_rpoly < $(INPUT_DIR)/test_solver4.in >> $(OUTPUT_DIR)/test_rpoly.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_rpoly < $(INPUT_DIR)/test_solver5.in >> $(OUTPUT_DIR)/test_rpoly.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_rpoly < $(INPUT_DIR)/test_solver6.in >> $(OUTPUT_DIR)/test_rpoly.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_rpoly < $(INPUT_DIR)/test_solver7.in >> $(OUTPUT_DIR)/test_rpoly.txt
	LD_LIBRARY_PATH=$(CXX_LIB_DIR):$$LD_LIBRARY_PATH $(BIN_DIR)/test_rpoly < $(INPUT_DIR)/test_solver8.in >> $(OUTPUT_DIR)/test_rpoly.txt
