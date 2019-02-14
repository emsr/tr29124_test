
TEST_BIN_DIR = bin

TEST_OUT_DIR = test_output

$(TEST_BIN_DIR)/test_chebyshev: $(TEST_BIN_DIR) test_chebyshev.cpp include/ext/chebyshev.h include/ext/chebyshev.tcc
	g++ -std=gnu++17 -I./include -I../include -I../polynomial/include -g -Wall -Wextra -o $(TEST_BIN_DIR)/test_chebyshev test_chebyshev.cpp -lquadmath

clean:
	rm -rf $(TEST_BIN_DIR)/*

test: $(TEST_OUT_DIR) $(TEST_BIN_DIR)/test_chebyshev
	./$(TEST_BIN_DIR)/test_chebyshev > $(TEST_OUT_DIR)/test_chebyshev.txt
	fi

$(TEST_BIN_DIR):
	if test ! -d $(TEST_BIN_DIR); then \
	  mkdir $(TEST_BIN_DIR); \
	fi

$(TEST_OUT_DIR):
	if test ! -d $(TEST_OUT_DIR); then \
	  mkdir $(TEST_OUT_DIR); \
	fi
